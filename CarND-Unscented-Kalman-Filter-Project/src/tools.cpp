#include <iostream>

#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  // create a container for the rmse calc
  VectorXd rmse(4);

  rmse << 0., 0., 0., 0.;

  // if estimation vectors don't have the right size or the number of estimations and 
  // ground truth points don't match
  if (   estimations[0].size() != 4 
  	  || estimations.size() != ground_truth.size()){
  	cerr << "Tools::CalculateRMSE exited due to improper input..." << endl;
    return rmse;
  }

  // define temp variable for storing rmse of each channel
  VectorXd res(4);

  // accummulate error in each channel
  for (size_t i=0; i < ground_truth.size(); i++){
  	// take vector difference
  	res = ground_truth[i] - estimations[i];
  	
  	// perform by-element multiplication
  	res = res.array() * res.array();

  	// do the accummulation
  	rmse += res;
  }

  // take the mean
  rmse /= ground_truth.size();

  // take the square root and return it
  return rmse.array().sqrt();
}
