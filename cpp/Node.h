#ifndef CLUSTERING_NODE_H
#define CLUSTERING_NODE_H

#include <armadillo>
#include <map>
#include <set>
#include <vector>

class Node {

public:

    // Variables indices to fix to 0
    std::vector<double> F0;
    // Variable indices to fix to 1
    std::vector<double> F1;

    // lower bound
    double lb;
    // upper bound
    double ub;
    // node id
    int id;

};


class SDPNode : public Node {


public:

    std::vector<arma::sp_mat> B_vector;
    arma::vec l_vector;

};

typedef struct NodeData {

    SDPNode *node;
    int i;

} NodeData;

typedef struct JobData {

    int type;
    NodeData *node_data;

} JobData;

typedef struct SDPResult {

    int n;
    double lb;
    double ub;
    arma::vec x_ub;
    int n_ineq;
    int cp_iter;
    int cp_flag;
    int branching_index;
    std::vector<arma::sp_mat> B_vector;
    arma::vec l_vector;

} SDPResult;


#endif //CLUSTERING_NODE_H
