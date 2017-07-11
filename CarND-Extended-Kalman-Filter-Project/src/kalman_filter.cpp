#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace::std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */

  // state prediction
  x_ = F_ * x_;

  // covariance prediction
  P_ = F_ * P_ * F_.transpose() + Q_;

  return;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  // compute measurement residual
  VectorXd res = z - H_*x_;

  // DEBUG
  if (false){
    cout << "" << endl;
    cout << z << endl;
    cout << "" << endl;
    cout << H_ << endl;
    cout << "" << endl;
    cout << x_ << endl;
    cout << "" << endl;
  }

  // compute measurement covariance
  MatrixXd S   = H_*P_*H_.transpose() + R_;

  // compute Kalman gain matrix
  MatrixXd K   = P_*H_.transpose()*S.inverse();

  // state update
  x_ += K * res;

  // covariance update
  MatrixXd I = MatrixXd::Identity(4,4);
  P_  = (I - K*H_) * P_;

  return;

}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  // compute measurement residual
  VectorXd res = z - h_(x_);

  // normalize phi difference to be between -pi and +pi
  if (res(1) < -M_PI){
    res(1) += 2.*M_PI;
  }
  else if (res(1) > M_PI){
    res(1) -= 2.*M_PI;
  }

  // compute measurement Jacobian
  Tools tools;
  MatrixXd H = tools.CalculateJacobian(x_);

  // compute measurement covariance
  MatrixXd S   = H*P_*H.transpose() + R_;

  // compute Kalman gain matrix
  MatrixXd K   = P_*H.transpose()*S.inverse();

  // state update
  x_ += K * res;

  // covariance update
  MatrixXd I = MatrixXd::Identity(4,4);
  P_  = (I - K*H) * P_;

  return;

}


// construct pseudomeasurement from state vector
VectorXd KalmanFilter::h_(const VectorXd &state){
  
  float  x = state(0);
  float  y = state(1);
  float vx = state(2);
  float vy = state(3);

  float r2 = x*x + y*y;
  float r  = sqrt(r2);
  float rv = x*vx + y*vy;

  VectorXd out(3);

  // output is radial distance, angle, and range rate
  out << r, atan2(y,x), rv / r;

  return out;
}