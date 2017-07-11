#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  time_us_        = 0;
  n_aug_          = 7;
  n_x_            = 5;
  nis_laser_      = 0.;
  nis_radar_      = 0.;
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2 * n_aug_ + 1);
  lambda_ = 3 - n_aug_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF: " << endl;


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double r   = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rd  = meas_package.raw_measurements_(2);

      x_(0) = r * cos(phi);
      x_(1) = r * sin(phi);
      x_(2) = rd;
      x_(3) = phi;
      x_(4) = 0.;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double x = meas_package.raw_measurements_(0);
      double y = meas_package.raw_measurements_(1);
      
      x_(0) = x;
      x_(1) = y;
      x_(2) = 0.;
      x_(3) = atan2(y,x);
      x_(4) = 0.;

    }

    P_  = MatrixXd::Identity(n_x_,n_x_);

    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    weights_.fill(0.);
    Xsig_pred_.fill(0.);
  
    for (unsigned i = 0; i < 2*n_aug_+1; i++){
      if (i == 0){
        weights_(i) = lambda_ / (lambda_ + double(n_aug_));
      }
      else{
        weights_(i) = 1. / (2. * (lambda_ + double(n_aug_)));
      }
    }
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.;
  
  Prediction(delta_t);
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

  time_us_ = meas_package.timestamp_;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.;
  x_aug(6) = 0.;
  
  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_ * std_a_,            0.          ,
              0.      , std_yawdd_ * std_yawdd_;
  P_aug.bottomRightCorner(2,2) = Q;
  
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  //set sigma points as columns of matrix Xsig
  // first point is the mean
  Xsig_aug.col(0) = x_aug;
  
  // now, for the rest, i and i+n_x can be filled via a loop:
  for (unsigned i=1; i <= n_aug_; i++){
    Xsig_aug.col(i)        = x_aug + sqrt(lambda_+double(n_aug_))*A.col(i-1);
    Xsig_aug.col(i+n_aug_) = x_aug - sqrt(lambda_+double(n_aug_))*A.col(i-1);
  }

  double tol = 1.e-3;
  //predict sigma points
  for (unsigned i = 0; i < Xsig_pred_.cols(); i++){
    VectorXd y    = Xsig_aug.col(i);
    VectorXd x    = y.head(n_x_);
    double yawd   = y(4);
    double dt2    = delta_t*delta_t;
    VectorXd x_up = VectorXd(n_x_);
    
    //avoid division by zero
    if (fabs(yawd) < tol){
      x_up(0) = y(2) * cos(y(3)) * delta_t + .5 * dt2 * cos(y(3)) * y(5);
      x_up(1) = y(2) * sin(y(3)) * delta_t + .5 * dt2 * sin(y(3)) * y(5);
      x_up(2) = delta_t * y(5);
      x_up(3) = .5 * dt2 * y(6);
      x_up(4) = delta_t * y(6);
      
    }
    else{
      double t = y(3) + y(4) * delta_t;
      x_up(0)  = (y(2) / y(4)) * ( sin(t) - sin(y(3)) ) + .5 * dt2 * y(5) * cos(y(3));
      x_up(1)  = (y(2) / y(4)) * (-cos(t) + cos(y(3)) ) + .5 * dt2 * y(5) * sin(y(3));
      x_up(2)  = delta_t * y(5);
      x_up(3)  = y(4) * delta_t + .5 *dt2 * y(6);
      x_up(4)  = delta_t * y(6);
    }
    //write predicted sigma points into right column
    Xsig_pred_.col(i) = x + x_up;
  }

  //predict state mean
  x_.fill(0.);
  for (unsigned i = 0; i < 2*n_aug_+1; i++){
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.);
  VectorXd res = VectorXd(n_x_);
  for (unsigned i = 0; i < 2*n_aug_+1; i++){
      res = Xsig_pred_.col(i) - x_;  
      P_ += weights_(i) * res * res.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.);

  //transform sigma points into measurement space
  VectorXd tmp = VectorXd(n_x_);
  double x, y;
  for (unsigned i = 0; i < 2*n_aug_+1; i++){
    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);
  }
  
  //calculate measurement covariance matrix S
  MatrixXd R   = MatrixXd(n_z,n_z);
  R << std_laspx_ * std_laspx_, 0.                     ,
       0.                     , std_laspy_ * std_laspy_;

  // call common filter Update routine
  Update(Zsig, R, meas_package);
  
}

/**
 * Updates the state and the state covariance matrix using a lidar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.);
  
  //transform sigma points into measurement space
  VectorXd tmp = VectorXd(n_x_);
  double x, y, x2, y2, phi, v, phi_do, r;
  for (unsigned i = 0; i < Xsig_pred_.cols(); i++){
    x   = Xsig_pred_(0,i);
    y   = Xsig_pred_(1,i);
    v   = Xsig_pred_(2,i);
    phi = Xsig_pred_(3,i);
    x2  = x * x;
    y2  = y * y;
    r   = sqrt(x2 + y2);

    Zsig(0,i) = r;
    Zsig(1,i) = atan2(y,x);
    Zsig(2,i) = (x * cos(phi) * v + y * sin(phi) * v) / r; 
  }

  //calculate measurement covariance matrix S
  MatrixXd R   = MatrixXd(n_z,n_z);
  R << std_radr_ * std_radr_,              0.          ,             0.         ,
                  0.        , std_radphi_ * std_radphi_,             0.         ,
                  0.        ,              0.          , std_radrd_ * std_radrd_;

  // call common filter Update routine
  Update(Zsig, R, meas_package);
}

void UKF::Update(MatrixXd Zsig, MatrixXd R, MeasurementPackage meas_package){
  
  // measurement size = n_z
  int n_z = Zsig.rows();

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.);

  //calculate mean predicted measurement
  for(int i = 0; i < Zsig.cols(); i++){
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  // calculate measurement covariance
  for(int i = 0; i < Zsig.cols(); i++){
    VectorXd z_diff = VectorXd(n_z);
    z_diff = Zsig.col(i) - z_pred;
    
    NormalizeAngle(z_diff(1));
    
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  S += R;
  
  //measurement
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  Tc.fill(0.);
  //calculate cross correlation matrix
  for(int i = 0; i < Zsig.cols(); i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3));
    
    VectorXd z_diff = Zsig.col(i) - z_pred;
    NormalizeAngle(z_diff(1));
    
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  //update state mean and covariance matrix
  
  VectorXd z_diff = z - z_pred;
  NormalizeAngle(z_diff(1));
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  double nis = z_diff.transpose() * S.inverse() * z_diff;
  
  //calculate the Normalized Innovation Square
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    nis_radar_= nis;
  }else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    nis_laser_= nis;
  }

}

void UKF::NormalizeAngle(double &angle){
  while (angle> M_PI) angle-=2.*M_PI;
  while (angle<-M_PI) angle+=2.*M_PI;
}
