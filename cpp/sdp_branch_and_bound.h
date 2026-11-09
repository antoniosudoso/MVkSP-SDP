#ifndef CLUSTERING_SDP_BRANCH_AND_BOUND_H
#define CLUSTERING_SDP_BRANCH_AND_BOUND_H

#include <armadillo>
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include "JobQueue.h"
#include "config_params.h"

typedef struct MatlabStruct {

    std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;
    matlab::data::ArrayFactory factory;

} MatlabStruct;


typedef struct SharedData {

    // Between workers and main
    std::condition_variable mainConditionVariable;
    std::vector<bool> threadStates;

    // Queue of requests waiting to be processed
    JobAbstractQueue *queue;
    // This condition variable is used for the threads to wait until there is work to do
    std::condition_variable queueConditionVariable;
    // Mutex to protect queue
    std::mutex queueMutex;

    double best_lb;
    double global_ub;
    arma::vec global_x;
    double global_gap;
    int n_nodes;

    double root_gap;
    double root_lb;
    double root_ub;
    int root_cp;
    double root_time;
    int root_cuts;

} SharedData;

typedef struct InputData {

    arma::mat A;
    int k;

} InputData;

class BBLog {

public:

    int n; // problem size
    int parent_id; // parent node id
    int id; // node id
    int type; // node type
    int br_idx_i; // branching index i
    double parent_lb; // parent lb
    double lb; // node lower bound
    int n_iter; // number of cutting-plane iterations
    int n_ineq; // number of inequalities
    double time; // elapsed time in seconds
    double ub; // upper bound at the node
    double gub; // global upper bound
    double gap; // gap at the node
    double ggap; // global gap
    int open; // number of open nodes

public:

    static void header();
    static void footer();
    void print() const;

};


void sdp_branch_and_bound(Config *config, InputData *input_data);
std::vector<JobData *> solve_root_problem(MatlabStruct *matlab, InputData *input_data, SharedData *shared_data, Config *config);
std::vector<JobData *> solve_child_problem(MatlabStruct *matlab, int type, NodeData *node_data, InputData *input_data, SharedData *shared_data, Config *config);

#endif //CLUSTERING_SDP_BRANCH_AND_BOUND_H
