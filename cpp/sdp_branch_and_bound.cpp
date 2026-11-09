#include <thread>
#include <iomanip>
#include "matlab_util.h"
#include "sdp_branch_and_bound.h"
#include "JobQueue.h"
#include "config_params.h"
#include "Node.h"
#include "ThreadPool.h"

void BBLog::header() {

    log_file << "\n"
              << "   --------------------------------------------------------------------"
                 "  SDP-BASED BRANCH-AND-BOUND ALGORITHM FOR MVkSP  "
                 "--------------------------------------------------------------------\n";

    log_file << "   "
              << std::setw(7)  << "N"        << "  "
              << std::setw(7) << "PAR ID"   << "  "
              << std::setw(7) << "ID"       << "  "
              << std::setw(6) << "TYPE"     << "  "
              << std::setw(7) << "BR IDX"   << "  "
              << std::setw(15) << "PAR LB"   << "  "
              << std::setw(15) << "LB"       << "  "
              << std::setw(6) << "ITER"     << "  "
              << std::setw(13) << "INEQ"     << "  "
              << std::setw(10) << "TIME (s)" << "  "
              << std::setw(15) << "UB"       << "  "
              << std::setw(15) << "GUB"      << "  "
              << std::setw(13) << "NGAP (%)" << "  "
              << std::setw(13) << "GGAP (%)" << "  "
              << std::setw(7) << "QUEUE"
              << "\n";

    log_file << "   "
              << std::string(7,  '-') << "  "
              << std::string(7, '-') << "  "
              << std::string(7, '-') << "  "
              << std::string(6, '-') << "  "
              << std::string(7, '-') << "  "
              << std::string(15, '-') << "  "
              << std::string(15, '-') << "  "
              << std::string(6, '-') << "  "
              << std::string(13, '-') << "  "
              << std::string(10, '-') << "  "
              << std::string(15, '-') << "  "
              << std::string(15, '-') << "  "
              << std::string(13, '-') << "  "
              << std::string(13, '-') << "  "
              << std::string(7, '-')
              << "\n";
}

void BBLog::footer() {
    log_file << "   --------------------------------------------------------------------"
                 "------------------------------------------------------"
                 "-----------------------------------------------------------------\n\n";
}

void BBLog::print() const {

    log_file << "   "
              << std::setw(7)  << this->n                                 << "  "
              << std::setw(7) << this->parent_id                          << "  "
              << std::setw(7) << this->id                                 << "  "
              << std::setw(6) << this->type                               << "  "
              << std::setw(7) << br_idx_i                                 << "  "
              << std::setw(15) << std::fixed << std::setprecision(4) << this->parent_lb << "  "
              << std::setw(15) << std::fixed << std::setprecision(4) << this->lb        << "  "
              << std::setw(6) << this->n_iter                             << "  "
              << std::setw(13) << this->n_ineq                            << "  "
              << std::setw(10) << std::fixed << std::setprecision(2) << this->time      << "  "
              << std::setw(15) << std::fixed << std::setprecision(4) << this->ub        << "  "
              << std::setw(15) << std::fixed << std::setprecision(4) << this->gub       << "  "
              << std::setw(13) << std::fixed << std::setprecision(4) << (this->gap  * 100.0) << "  "
              << std::setw(13) << std::fixed << std::setprecision(4) << (this->ggap * 100.0) << "  "
              << std::setw(7) << this->open
              << "\n";
    log_file.flush();
}

// Root node
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory, InputData *input_data, Config *config) {

    matlab::data::TypedArray<double> A_matlab = arma_to_matlab_matrix(factory, input_data->A);
    matlab::data::TypedArray<double> k_matlab = factory.createScalar<double>(input_data->k);

    // Create StructArray
    std::vector<std::string> f = {"n_threads", "bb_tol", "sdp_fr", "sdp_verbose", "sdp_tol",
                                  "cp_maxiter", "cp_tol", "cp_maxineq", "cp_newineq", "cp_epsineq", "cp_activeineq"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1},f);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(config->matlab_session_threads_root);
    struct_matlab[0]["bb_tol"] = factory.createScalar<double>(config->bb_tol);
    struct_matlab[0]["sdp_fr"] = factory.createScalar<int>(config->sdp_fr);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(config->sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(config->sdp_solver_tol);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(config->cp_max_iter);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(config->cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(config->cp_max_ineq);
    struct_matlab[0]["cp_newineq"] = factory.createScalar<double>(config->cp_new_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(config->cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(config->cp_eps_active);

    std::vector<matlab::data::Array> args({A_matlab, k_matlab, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    std::vector<matlab::data::Array> result;
    result = matlabPtr->feval(u"call_solve_MVkSP_root",n_return, args);
    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::TypedArray<double> field_best_x_gub = structArray[0]["best_x_ub"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::TypedArray<double> field_ineq_list = structArray[0]["ineq_list"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_Bcell"];
    matlab::data::TypedArray<double> field_best_l = structArray[0]["best_l"];
    matlab::data::TypedArray<double> field_br_idx = structArray[0]["br_idx"];
    matlab::data::TypedArray<double> field_n = structArray[0]["n"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::vec x_heuristic = matlab_to_arma_vector(field_best_x_gub);
    arma::vec l_vector = matlab_to_arma_vector(field_best_l);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    int n_ineq = (int) field_ineq_list[cp_iter];
    std::vector<arma::sp_mat> B_vector = matlab_to_arma_vector_sp_mat(field_best_B_cell);
    int br_idx = (int) field_br_idx[0];
    int n = (int) field_n[0];

    return SDPResult{n, lb, ub, x_heuristic, n_ineq, cp_iter, cp_flag, br_idx, B_vector, l_vector};
}

// Branching nodes
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory,
        std::vector<arma::sp_mat> &parent_B_vector, arma::vec &parent_l_vector,
        double global_ub, arma::vec &global_x, std::vector<double> &global_F0, std::vector<double> &global_F1,
        InputData *input_data, Config *config) {

    // convert data
    matlab::data::TypedArray<double> A_matlab = arma_to_matlab_matrix(factory, input_data->A);
    matlab::data::TypedArray<double> k_matlab = factory.createScalar<double>(input_data->k);
    matlab::data::TypedArray<double> matlab_F0 = vector_to_matlab_vector(factory, reinterpret_cast<std::vector<double> &>(global_F0));
    matlab::data::TypedArray<double> matlab_F1 = vector_to_matlab_vector(factory, reinterpret_cast<std::vector<double> &>(global_F1));
    matlab::data::CellArray matlab_B_cell = arma_to_matlab_cell(factory, parent_B_vector);
    matlab::data::TypedArray<double> matlab_parent_l = arma_to_matlab_vector(factory, parent_l_vector);
    matlab::data::TypedArray<double> matlab_global_x = arma_to_matlab_vector(factory, global_x);
    matlab::data::TypedArray<double> matlab_gub = factory.createScalar<double>(global_ub);

    // Create StructArray
    std::vector<std::string> f = {"n_threads", "bb_tol", "sdp_fr", "sdp_verbose", "sdp_tol",
                                  "cp_maxiter", "cp_tol", "cp_maxineq", "cp_newineq", "cp_epsineq", "cp_activeineq"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1},f);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(config->matlab_session_threads_root);
    struct_matlab[0]["bb_tol"] = factory.createScalar<double>(config->bb_tol);
    struct_matlab[0]["sdp_fr"] = factory.createScalar<int>(config->sdp_fr);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(config->sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(config->sdp_solver_tol);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(config->cp_max_iter);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(config->cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(config->cp_max_ineq);
    struct_matlab[0]["cp_newineq"] = factory.createScalar<double>(config->cp_new_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(config->cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(config->cp_eps_active);


    std::vector<matlab::data::Array> args({A_matlab, k_matlab, matlab_F0, matlab_F1, matlab_gub, matlab_global_x,
                                           matlab_B_cell, matlab_parent_l, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    std::vector<matlab::data::Array> result;
    result = matlabPtr->feval(u"call_solve_MVkSP_child", n_return, args);

    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::TypedArray<double> field_best_x_gub = structArray[0]["best_x_ub"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::TypedArray<double> field_ineq_list = structArray[0]["ineq_list"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_Bcell"];
    matlab::data::TypedArray<double> field_best_l = structArray[0]["best_l"];
    matlab::data::TypedArray<double> field_br_idx = structArray[0]["br_idx"];
    matlab::data::TypedArray<double> field_n = structArray[0]["n"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::vec x_heuristic = matlab_to_arma_vector(field_best_x_gub);
    arma::vec l_vector = matlab_to_arma_vector(field_best_l);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    int n_ineq = (int) field_ineq_list[cp_iter];
    std::vector<arma::sp_mat> B_vector = matlab_to_arma_vector_sp_mat(field_best_B_cell);
    int br_idx = (int) field_br_idx[0];
    int n = (int) field_n[0];

    return SDPResult{n, lb, ub, x_heuristic, n_ineq, cp_iter, cp_flag, br_idx, B_vector, l_vector};

}



std::vector<JobData *> create_children_jobs(double node_gap, SDPNode *node, int branching_index,
                                            NodeData *parent, SharedData *shared_data, Config *config) {

    if (std::isinf(node->lb) || node_gap <= config->bb_tol) {
        delete(node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return {};
    }

    if (branching_index == -1) {

        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        log_file << "\t PRUNING BY OPTIMALITY " << node->id << "\n";
        delete (node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return {};

        // mutex is automatically released when lock goes out of scope
    }

    std::vector<JobData *> branch;

    auto *child_data_zero = new NodeData();
    child_data_zero->node = new SDPNode(*node);
    child_data_zero->i = branching_index;

    auto *child_data_one = new NodeData();
    child_data_one->node = new SDPNode(*node);
    child_data_one->i = branching_index;

    auto *zero_job_data = new JobData();
    zero_job_data->type = ZERO;
    zero_job_data->node_data = child_data_zero;

    auto *one_job_data = new JobData();
    one_job_data->type = ONE;
    one_job_data->node_data = child_data_one;

    branch.emplace_back(zero_job_data);
    branch.emplace_back(one_job_data);

    if (parent != nullptr) {
        delete (parent->node);
        delete (parent);
    }

    delete (node);

    return branch;

}


std::vector<JobData *> solve_child_problem(MatlabStruct *matlab, int type, NodeData *node_data, InputData *input_data, SharedData  *shared_data, Config *config) {

    auto child_node = new SDPNode();
	auto parent = node_data->node;

	double parent_gap = (shared_data->global_ub - parent->lb) / std::max(1.0, shared_data->global_ub);
	if (parent_gap <= config->bb_tol)
        return {};

    if (type == ZERO) {
        child_node->F0 = parent->F0;
        child_node->F1 = parent->F1;
        child_node->F0.push_back(node_data->i);

        /*
        std::cout << "\n\t C++ printing F0: ";
        for (auto &elem : child_node->F0) std::cout << elem << " ";
        std::cout << "\n\n";
        */

    }

    if (type == ONE) {
        child_node->F0 = parent->F0;
        child_node->F1 = parent->F1;
        child_node->F1.push_back(node_data->i);

        /*
        std::cout << "\n\t C++ printing F1: ";
        for (auto &elem : child_node->F1) std::cout << elem << " ";
        std::cout << "\n\n";
        */

    }

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult  sdp_result = solve_sdp(matlab->matlabPtr, matlab->factory,
              parent->B_vector, parent->l_vector, shared_data->global_ub,
              shared_data->global_x, child_node->F0, child_node->F1, input_data, config);

    int n = sdp_result.n;
    int cp_iter = sdp_result.cp_iter;
    // int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    child_node->lb = std::max(sdp_result.lb, parent->lb);
    child_node->ub = sdp_result.ub;
    child_node->B_vector = sdp_result.B_vector;
    child_node->l_vector = sdp_result.l_vector;
    arma::vec x_heuristic = sdp_result.x_ub;
    int br_idx = sdp_result.branching_index;

    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        if (child_node->ub - shared_data->global_ub <= -config->bb_tol) {
            // update global upper bound
            shared_data->global_ub = child_node->ub;
            shared_data->global_x = x_heuristic;
        }

        child_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;
        int open = shared_data->queue->getSize();

        node_gap = (shared_data->global_ub - child_node->lb) / std::max(1.0, shared_data->global_ub);
        Node *min_lb_node = shared_data->queue->getMinLB();
        if (min_lb_node != nullptr && min_lb_node->lb > shared_data->best_lb) {
            shared_data->global_gap = (shared_data->global_ub - min_lb_node->lb) / std::max(1.0, shared_data->global_ub);
            shared_data->best_lb = min_lb_node->lb;
        }

        auto end_all = std::chrono::high_resolution_clock::now();
        auto elapsed = (std::chrono::duration_cast<std::chrono::microseconds> (end_all - start_time).count());
        auto time = (double) elapsed / 1000000;

        BBLog log{n, parent->id, child_node->id, type, br_idx,
                  parent->lb, child_node->lb, cp_iter, n_ineq, time, child_node->ub,
                  shared_data->global_ub, node_gap, shared_data->global_gap, open};
        log.print();

    }

    return create_children_jobs(node_gap, child_node, br_idx, node_data, shared_data, config);
}


std::vector<JobData *> solve_root_problem(MatlabStruct *matlab, InputData *input_data, SharedData *shared_data, Config *config) {

    // init root
    SDPNode *root;
    root = new SDPNode();
    root->F0 = {};
    root->F1 = {};
    root->id = shared_data->n_nodes;

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult sdp_result = solve_sdp(matlab->matlabPtr, matlab->factory, input_data, config);

    int n = sdp_result.n;
    int cp_iter = sdp_result.cp_iter;
    // int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    root->lb = sdp_result.lb;
    root->ub = sdp_result.ub;
    root->B_vector = sdp_result.B_vector;
    root->l_vector = sdp_result.l_vector;
    arma::vec x_ub = sdp_result.x_ub;
    int br_idx = sdp_result.branching_index;

    // Update shared data
    shared_data->best_lb = root->lb;
    shared_data->global_ub = root->ub;
    shared_data->global_x = x_ub;
    shared_data->n_nodes++;

    int open = (int) shared_data->queue->getSize();
    double inf = std::numeric_limits<double>::infinity();

    double node_gap = (shared_data->global_ub - root->lb) / std::max(1.0, shared_data->global_ub);
    shared_data->global_gap = node_gap;

    auto end_all = std::chrono::high_resolution_clock::now();
    auto elapsed = (std::chrono::duration_cast<std::chrono::microseconds> (end_all - start_time).count());
    auto time = (double) elapsed / 1000000;

    shared_data->root_ub = shared_data->global_ub;
    shared_data->root_lb = shared_data->best_lb;
    shared_data->root_gap = shared_data->global_gap;
    shared_data->root_time = time;
    shared_data->root_cp = cp_iter;
    shared_data->root_cuts = n_ineq;

    BBLog log{n, -1, root->id, ROOT, br_idx,
              -inf, root->lb, cp_iter, n_ineq, time, root->ub,
              shared_data->global_ub, node_gap, shared_data->global_gap, open};
    log.print();

    return create_children_jobs(node_gap, root, br_idx, nullptr, shared_data, config);
}

bool is_thread_pool_working(std::vector<bool> &thread_state) {
    int count = 0;
    for (auto && i : thread_state) {
        if (i) count++;
    }
    if (count == 0)
        return false;
    return true;
}


void sdp_branch_and_bound(Config *config, InputData *input_data) {

    int n_thread = config->bb_threads;

    log_file << "\n";
    auto **matlab = new MatlabStruct *[n_thread];
    for (int i = 0; i < n_thread; i++) {
        auto *matlab_struct = new MatlabStruct();
        log_file << " Starting MATLAB session..." << "\n";
        matlab_struct->matlabPtr = start_matlab(config);
        log_file << " MATLAB started!" << "\n";
        matlab[i] = matlab_struct;
    }
    log_file << "\n";

    JobAbstractQueue *queue;
    switch (config->bb_visiting_strategy) {
        case DEPTH_FIRST:
            queue = new JobStack();
            break;
        case BEST_FIRST:
            queue = new JobPriorityQueue();
            break;
        case BREADTH_FIRST:
            queue = new JobQueue();
            break;
        default:
            queue = new JobPriorityQueue();
    }

    auto *shared_data = new SharedData();
    shared_data->best_lb = -std::numeric_limits<double>::infinity();
    shared_data->global_gap = std::numeric_limits<double>::infinity();
    shared_data->global_ub = std::numeric_limits<double>::infinity();
    shared_data->n_nodes = 0;
    shared_data->queue = queue;

    shared_data->root_time = 0;
    shared_data->root_lb = - std::numeric_limits<double>::infinity();
    shared_data->root_ub = std::numeric_limits<double>::infinity();
    shared_data->root_gap = std::numeric_limits<double>::infinity();
    shared_data->root_cp = 0;

    shared_data->threadStates.reserve(n_thread);
    for (int i = 0; i < n_thread; i++) {
        shared_data->threadStates.push_back(false);
    }

    ThreadPool pool(matlab, config, shared_data, input_data, n_thread);

    BBLog::header();
    log_file.flush();

    auto start_all = std::chrono::high_resolution_clock::now();

    std::vector<JobData *> jobs = solve_root_problem(matlab[0], input_data, shared_data, config);

    // submit jobs to the thread pool
    for (auto & job : jobs) {
        pool.addJob(job);
    }

    while (true) {

        {
            std::unique_lock<std::mutex> l(shared_data->queueMutex);
            while (is_thread_pool_working(shared_data->threadStates) && shared_data->n_nodes < config->bb_max_nodes &&
                   std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_all).count() <
                   config->bb_time_limit) {
                shared_data->mainConditionVariable.wait(l);
            }

            if (shared_data->queue->empty() || shared_data->n_nodes >= config->bb_max_nodes ||
                std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_all).count() >=
                config->bb_time_limit) {
                break;
            }

        }

    }

    auto end_all = std::chrono::high_resolution_clock::now();
    auto elapsed = (std::chrono::duration_cast<std::chrono::microseconds> (end_all - start_all).count());

    pool.quitPool();

    if (queue->empty()) {
        shared_data->global_gap = 0.0;
    } else {
        double valid_lb = queue->getMinLB()->lb;
        shared_data->global_gap = (shared_data->global_ub - valid_lb) / std::max(1.0, shared_data->global_ub);
    }

    BBLog::footer();
    log_file.flush();

    log_file << "\n";
    log_file << " Root time (s): " << shared_data->root_time << "\n";
    log_file << " Root LB: " << shared_data->root_lb << "\n";
    log_file << " Root UB: " << shared_data->root_ub << "\n";
    log_file << " Root Gap (%): " << shared_data->root_gap * 100 << "\n";
    log_file << " Root iter: " << shared_data->root_cp << "\n";
    log_file << " Root cuts: " << shared_data->root_cuts << "\n";

    log_file << "\n";
    log_file << " Objective: " << shared_data->global_ub << "\n";
    log_file << " Gap (%): " << shared_data->global_gap * 100 << "\n";
    log_file << " Nodes: " << shared_data->n_nodes << "\n";
    log_file << " Time (s): " << (double) elapsed / 1000000 << "\n";
    log_file << " Solution: ";
    int n = (int) shared_data->global_x.size();
    for (int i = 0; i < n; i++) {
        if (shared_data->global_x(i) == 1.0)
            log_file << i << " ";
    }
    log_file << "\n\n";

    // Free memory
    delete (queue);
    delete (shared_data);
    for (int i = 0; i < n_thread; i++) {
        delete (matlab[i]);
    }
    delete[] (matlab);

}
