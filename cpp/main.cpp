#include <iostream>
#include <map>
#include <algorithm>
#include <armadillo>
#include "sdp_branch_and_bound.h"
#include "config_params.h"

std::ofstream log_file;

std::map<std::string, std::string> read_params(std::string &config_file) {

    std::map<std::string, std::string> config_map = {};

    std::ifstream cFile (config_file);
    if (cFile.is_open()) {
        std::string line;
        while (getline(cFile, line)){
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find('=');
            auto key = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            config_map.insert(std::pair<std::string, std::string>(key, value));
        }

    }
    else {
        std::cerr << "Couldn't open config file for reading.\n";
    }

    return config_map;
}

// read data matrix
arma::mat read_data(const char *filename, int &n, int &d) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << strerror(errno) << "\n";
        exit(EXIT_FAILURE);
    }
    // read the header n, d
    file >> n >> d;
    arma::mat A(n, d);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < d; j++) {
            file >> A(i, j);
        }
    }
    return A;
}

void run(int argc, char **argv) {

    if (argc < 5) {
        std::cerr << "Usage: <config.txt> <DATA_PATH> <K> <LOG_PATH>" << std::endl;
        exit(EXIT_FAILURE);
    }

    auto *config = new Config();

    std::string config_file = argv[1];
    config->data_path = argv[2];
    int k = std::stoi(argv[3]);
    config->log_path = argv[4];
    log_file.open(config->log_path);

    std::map<std::string, std::string> config_map = read_params(config_file);

    int n, d;
    arma::mat W = read_data(config->data_path, n, d);
    auto *input_data = new InputData();
    input_data->A = W.t();
    input_data->k = k;

    //std::cout << A.n_rows << " " << A.n_cols << "\n";

    // Branch-and-bound settings
    config->bb_tol = std::stod(config_map["BB_TOL"]);
    config->bb_threads = std::stoi(config_map["BB_THREADS"]);
    config->bb_max_nodes = std::stoi(config_map["BB_MAX_NODES"]);
    config->bb_time_limit = std::stod(config_map["BB_TIME_LIMIT"]);
    config->bb_visiting_strategy = std::stoi(config_map["BB_VISITING_STRATEGY"]);

    // Matlab
    config->matlab_session_threads_root = std::stoi(config_map["MATLAB_SESSION_THREADS_ROOT"]);
    config->matlab_session_threads_child = std::stoi(config_map["MATLAB_SESSION_THREADS_CHILD"]);
    config->matlab_source_folder = config_map["MATLAB_SOURCE_FOLDER"].c_str();

    // SDP solver
    config->sdp_fr = std::stoi(config_map["SDP_FR"]);
    config->sdp_solver_tol = std::stod(config_map["SDP_SOLVER_TOL"]);
    config->sdp_solver_verbose = std::stoi(config_map["SDP_SOLVER_VERBOSE"]);

    // Cutting-plane settings
    config->cp_max_iter = std::stoi(config_map["CP_MAX_ITER"]);
    config->cp_tol = std::stod(config_map["CP_TOL"]);
    config->cp_max_ineq = std::stoi(config_map["CP_MAX_INEQ"]);
    config->cp_new_ineq = std::stoi(config_map["CP_NEW_INEQ"]);
    config->cp_eps_ineq = std::stod(config_map["CP_EPS_INEQ"]);
    config->cp_eps_active = std::stod(config_map["CP_EPS_ACTIVE"]);

    // Additional solvers
    config->gurobi_folder = config_map["GUROBI_FOLDER"].c_str();
    config->sdpnal_folder = config_map["SDPNAL_FOLDER"].c_str();

    log_file << "\n" << " DATA_PATH: " << config->data_path <<"\n\n";
    log_file << " N: " << n << "\n";
    log_file << " D: " << d << "\n";
    log_file << " K: " << k << "\n\n";

    log_file << " BB_TOL: " << config->bb_tol << "\n";
    log_file << " BB_THREADS: " << config->bb_threads << "\n";
    log_file << " BB_MAX_NODES: " << config->bb_max_nodes << "\n";
    log_file << " BB_TIME_LIMIT: " << config->bb_time_limit << "\n";
    log_file << " BB_VISITING_STRATEGY: " << config->bb_visiting_strategy << "\n\n";

    log_file << " MATLAB_SESSION_THREADS_ROOT: " << config->matlab_session_threads_root << "\n";
    log_file << " MATLAB_SESSION_THREADS_CHILD: " << config->matlab_session_threads_child << "\n";
    log_file << " MATLAB_SOURCE_FOLDER: " << config->matlab_source_folder << "\n\n";

    log_file << " SDP_FR: " << config->sdp_fr << "\n";
    log_file << " SDP_SOLVER_TOL: " << config->sdp_solver_tol << "\n";
    log_file << " SDP_SOLVER_VERBOSE: " << config->sdp_solver_verbose << "\n\n";

    log_file << " CP_MAX_ITER: " << config->cp_max_iter << "\n";
    log_file << " CP_TOL: " << config->cp_tol << "\n";
    log_file << " CP_MAX_INEQ: " << config->cp_max_ineq << "\n";
    log_file << " CP_NEW_INEQ: " << config->cp_new_ineq << "\n";
    log_file << " CP_EPS_INEQ: " << config->cp_eps_ineq << "\n";
    log_file << " CP_EPS_ACTIVE: " << config->cp_eps_active << "\n\n";

    log_file << " GUROBI_FOLDER: " << config->gurobi_folder << "\n";
    log_file << " SDPNAL_FOLDER: " << config->sdpnal_folder << "\n\n";

    sdp_branch_and_bound(config, input_data);

    log_file.close();

    delete (input_data);
    delete (config);

}

int main(int argc, char **argv) {

    // The B&B algorithm for solving the nonconvex quadratic programming problem:
    //
    // min_{x \in \{0,1\}^n} 1/k \sum_{i=1}^n x_i \|p_i - \mu(x)\|^2
    // s.t. sum(x) = k

    run(argc, argv);

    exit(EXIT_SUCCESS);
}
