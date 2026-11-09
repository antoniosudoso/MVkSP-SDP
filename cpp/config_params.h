#ifndef CLUSTERING_CONFIG_PARAMS_H
#define CLUSTERING_CONFIG_PARAMS_H

#define BEST_FIRST 0
#define DEPTH_FIRST 1
#define BREADTH_FIRST 2

#define ROOT -1
#define ZERO 0
#define ONE 1

extern std::ofstream log_file;

typedef struct Config {

    const char *data_path;
    const char *log_path;

    // Branch-and-bound settings
    double bb_tol;
    int bb_threads;
    int bb_max_nodes;
    double bb_time_limit;
    int bb_visiting_strategy;

    // Matlab settings
    int matlab_session_threads_root;
    int matlab_session_threads_child;
    const char *matlab_source_folder;

    // SDP settings
    int sdp_fr;
    double sdp_solver_tol;
    int sdp_solver_verbose;

    // Additional solvers
    const char *gurobi_folder;
    const char *sdpnal_folder;

    // cutting plane setting
    int cp_max_iter;
    double cp_tol;
    int cp_max_ineq;
    double cp_new_ineq;
    double cp_eps_ineq;
    double cp_eps_active;

} Config;


#endif //CLUSTERING_CONFIG_PARAMS_H
