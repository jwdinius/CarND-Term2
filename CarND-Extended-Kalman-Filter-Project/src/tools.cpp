#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

using namespace std;

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

  // at this point, it's the mean squared error (it's enough to debug)
  if (false){
    cout << rmse << endl;
    cout << "" << endl;
  }
  // take the square root and return it
  return rmse.array().sqrt();

}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  // calculate a jacobian matrix for the measurement process h
  MatrixXd H(3,4);

  // get state
  float x  = x_state(0);
  float y  = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // avoid divide-by-zero
  if (x == 0. && y == 0.){
    cerr << "position is (0,0) giving bad range info" << endl;
    H << 0., 0., 0., 0.,
         0., 0., 0., 0.,
         0., 0., 0., 0.;

    return H;
  }

  // define helper variables
  float r2 = x*x + y*y;
  float r  = sqrt(r2);
  float r3 = r2*r;
  float xv = x*vy - y*vx;

  H <<  x / r    ,   y/r     ,    0.,    0.,
       -y / r2   ,   x/r2    ,    0.,    0.,
       -y*xv / r3,  x*xv / r3,   x/r,   y/r;

  return H;

}
