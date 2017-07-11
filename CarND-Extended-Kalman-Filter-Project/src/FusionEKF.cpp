#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  Q_ = MatrixXd::Identity(4,4);

  VectorXd x(4);
  x << 0., 0., 0., 0.;
  MatrixXd P = MatrixXd::Identity(4,4);
  P *= 1000.;
  MatrixXd F = MatrixXd::Identity(4,4);
  F(0,2) = 1.;
  F(1,3) = 1.;

  H_laser_ << 1. , 0. ,  0. ,  0.,
              0. , 1. ,  0. ,  0.;

  // just do initialization on laser matrices, we can overwrite later if needed
  ekf_.Init(x, P, F, H_laser_, R_laser_, Q_);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_ = hinv_(measurement_pack.raw_measurements_);
      ekf_.R_ = R_radar_;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1),
                 0., 0.;
      ekf_.R_ = R_laser_;
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  if (false){
    cout << ekf_.F_ << endl;
  }
  
  float dt2 = dt*dt;
  float dt3 = dt2*dt;
  float dt4 = dt3*dt;
  
  float noise_ax = 9.;
  float noise_ay = noise_ax;
  
  ekf_.Q_(0,0) = noise_ax * dt4 / 4.;
  ekf_.Q_(0,2) = noise_ax * dt3 / 2.;
  ekf_.Q_(1,1) = noise_ay * dt4 / 4.;
  ekf_.Q_(1,3) = noise_ay * dt3 / 2.;
  ekf_.Q_(2,0) = ekf_.Q_(0,2);
  ekf_.Q_(2,2) = noise_ax * dt2;
  ekf_.Q_(3,1) = ekf_.Q_(1,3);
  ekf_.Q_(3,3) = noise_ay * dt2;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output when debugging
  if (false){
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
  }

  // update the previous timestamp for the next pass
  previous_timestamp_ = measurement_pack.timestamp_;

  return;
}

// compute measurement to state function (inverse of pseudomeasurement function)
VectorXd FusionEKF::hinv_(const VectorXd &z){
  float r    = z(1);
  float phi  = z(2);
  float rdot = z(3);

  VectorXd out(4);

  float x  = r*cos(phi);
  float y  = r*sin(phi);
  float vx = rdot*cos(phi);
  float vy = rdot*sin(phi);

  out << x, y, vx, vy;

  return out;

}
