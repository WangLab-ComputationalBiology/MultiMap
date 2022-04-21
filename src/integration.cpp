#include <RcppEigen.h>
#include <progress.hpp>
#include <unordered_map>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

typedef Eigen::Triplet<double> T;

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> FindWeightsC(NumericVector cells2,
                                         Eigen::MatrixXd distances,
                                         std::vector<std::string> anchor_cells2,
                                         std::vector<std::string> integration_matrix_rownames,
                                         Eigen::MatrixXd cell_index,
                                         Eigen::VectorXd anchor_score,
                                         double min_dist,
                                         double sd,
                                         bool display_progress) {
  std::vector<T> tripletList;
  tripletList.reserve(anchor_cells2.size() * 10);
  std::unordered_map<int, std::vector<int>> cell_map;
  Progress p(anchor_cells2.size() + cells2.size() , display_progress);
  // build map from anchor_cells2 to integration_matrix rows
  for(int i=0; i<anchor_cells2.size(); ++i){
    std::vector<int> matches;
    std::vector<std::string>::iterator iter = integration_matrix_rownames.begin();

    // what is cell map ?  ////////////////////////////////////////////////////
    // multiple same anchor ?? ////////////////////////////////////////////////
    // 1 ref anchor <=> N query anchor, thus, query anchors could occcur multiple times

    while ((iter = std::find(iter, integration_matrix_rownames.end(), anchor_cells2[i])) != integration_matrix_rownames.end()) {
      int idx = std::distance(integration_matrix_rownames.begin(), iter);
      matches.push_back(idx);
      iter++;
    }

    // anchor i : all its matches?? ///////////////////////////////////////////
    cell_map[i] = matches;
    p.increment();
  }

  // Construct dist_weights matrix
  for(auto const &cell : cells2){

    Eigen::VectorXd dist = distances.row(cell);
    Eigen::VectorXd indices = cell_index.row(cell);
    int k=0; //number of anchors used so far; a cell in the neighbor list may contribute to multiple anchors
    for(int i=0; i<indices.size() && k<indices.size(); ++i){ //index in neighbor list
      std::vector<int> mnn_idx = cell_map[indices[i]-1];
      for(int j=0; j<mnn_idx.size() && k<indices.size(); ++j){

        // anchor score may have multiple value for the same query anchor cell
        double to_add = 1 - exp(-1 * dist[i] * anchor_score[mnn_idx[j]]/ pow(2/sd, 2));

        tripletList.push_back(T(mnn_idx[j], cell, to_add));
        k++;
      }
    }
    p.increment();

  }

  Eigen::SparseMatrix<double> return_mat;

  if(min_dist == 0){

    Eigen::SparseMatrix<double> dist_weights(integration_matrix_rownames.size(), cells2.size());

    dist_weights.setFromTriplets(tripletList.begin(), tripletList.end(), [] (const double&, const double &b) { return b; });

    Eigen::VectorXd colSums = dist_weights.transpose() * Eigen::VectorXd::Ones(dist_weights.rows());

    for (int k=0; k < dist_weights.outerSize(); ++k){

      for (Eigen::SparseMatrix<double>::InnerIterator it(dist_weights, k); it; ++it){
        it.valueRef() = it.value()/colSums[k];
      }

    }

    return_mat = dist_weights;

  } else {

    Eigen::MatrixXd dist_weights = Eigen::MatrixXd::Constant(integration_matrix_rownames.size(), cells2.size(), min_dist);

    for(int i = 0; i < dist_weights.cols(); ++i){
      for(int j = 0; j < dist_weights.rows(); ++j){
        dist_weights(j, i) = 1 - exp(-1 * dist_weights(j, i) * anchor_score[j]/ pow(2/sd, 2) );
      }
    }

    for(auto const &weight : tripletList){
      dist_weights(weight.row(), weight.col()) = weight.value();
    }

    Eigen::VectorXd colSums = dist_weights.colwise().sum();

    for(int i = 0; i < dist_weights.cols(); ++i){
      for(int j = 0; j < dist_weights.rows(); ++j){
        dist_weights(j, i) = dist_weights(j, i) / colSums[i];
      }
    }

    return_mat = dist_weights.sparseView();
  }

  return(return_mat);
}
