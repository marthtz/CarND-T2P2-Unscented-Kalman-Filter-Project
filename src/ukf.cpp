#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define PI_X2 (2.*M_PI)
#define DIV0_LIMIT_CHECK (0.0000001)

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  ///* State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  /// Number of sigma points
  n_sigma_points_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;
  std_a_square_ = std_a_ * std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 4;
  std_yawdd_square_ = std_yawdd_ * std_yawdd_;

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

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_points_);

  // Weights of sigma points
  weights_ = VectorXd(n_sigma_points_);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<n_sigma_points_; i++)
  {
    //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // Measurement dimensions
  n_z_radar_ = 3;
  n_z_lidar_ = 2;

  // Measurement covariance matrices
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(n_z_lidar_, n_z_lidar_);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

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

  /*****************************************************************************
   *  Initialize
   ****************************************************************************/
  if (!is_initialized_)
  {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    // Initialize in case there's no RADAR or LIDAR measurement 
    x_.fill(1.0);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Get all 3 RADAR measurements
      double rho_measured = meas_package.raw_measurements_[0];
      double phi_measured = meas_package.raw_measurements_[1];
      double rhodot_measured = meas_package.raw_measurements_[2];

      // Convert from polar to cartesian coordinates
      double x = rho_measured * cos(phi_measured);
      double y = rho_measured * sin(phi_measured);

      // According to tips and tricks, there's not enough information to set
      // vx and vy with a RADAR initialization,. Therefore, set to 0.
      double vx = 0;// rhodot_measured * cos(phi_measured);
      double vy = 0;// rhodot_measured * sin(phi_measured);

      // Set state vector with first measurement
      x_ << x,  y, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity

      // Get both LIDAR measurements - set vx and vy to 0
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // Set first timestamp
    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // Compute the time elapsed between the current and previous measurements
  // dt is expressed in seconds
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Start a new prediction
  Prediction(dt);


  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates with current measurement
    UpdateRadar(meas_package);
  }
  else
  {
    // Laser updates with current measurement
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Loop counter
  int i;

  /*****************************************************************************
   * Augment Sigma Points
   ****************************************************************************/
  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  // Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_points_);

  // Create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_square_;
  P_aug(6,6) = std_yawdd_square_;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (i=0; i<n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  /*****************************************************************************
   * Predict Sigma Points
   ****************************************************************************/
  for (i=0; i<n_sigma_points_; i++)
  {
    // Extract values for better readability
    double p_x  = Xsig_aug(0,i);
    double p_y  = Xsig_aug(1,i);
    double v    = Xsig_aug(2,i);
    double yaw  = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // Predicted state values
    double px_p, py_p;

    // Avoid div 0
    if (fabs(yawd) > DIV0_LIMIT_CHECK)
    {
      px_p = p_x + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw+yawd*delta_t));
    }
    else
    {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    // Add noise
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // Write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  /*****************************************************************************
   * Predict Mean and Covariance
   ****************************************************************************/
  // Predicted state mean
  x_.fill(0.0);

  // Iterate over sigma points
  for (i=0; i<n_sigma_points_; i++)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // Predicted state covariance matrix
  P_.fill(0.0);

  // Iterate over sigma points
  for (i=0; i<n_sigma_points_; i++)
  {
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= PI_X2;
    while (x_diff(3)<-M_PI) x_diff(3) += PI_X2;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  // Loop counter
  int i;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_lidar_, n_sigma_points_);

  // Transform sigma points into measurement space
  for (i=0; i<n_sigma_points_; i++)
  {
    // Measurement model
    Zsig(0, i) = Xsig_pred_(0,i);     // x
    Zsig(1, i) = Xsig_pred_(1,i);     // y
  }

  // Call common update code
  UpdateUKF(meas_package, R_lidar_, n_z_lidar_, Zsig);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  // Loop counter
  int i;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, n_sigma_points_);

  // Transform sigma points into measurement space
  for (i=0; i<n_sigma_points_; i++)
  {
    // Extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // Measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                       //r
    Zsig(1, i) = atan2(p_y, p_x);                               //phi
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  // Call common update code
  UpdateUKF(meas_package, R_radar_, n_z_radar_, Zsig);
}



void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd R, int n_z, MatrixXd Zsig)
{
  // Loop counter
  int i;

  /*****************************************************************************
   * Predict Measurement
   ****************************************************************************/
  // Incoming measurement
  VectorXd z = meas_package.raw_measurements_;

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (i=0; i<n_sigma_points_; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (i=0; i<n_sigma_points_; i++)
  {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= PI_X2;
    while (z_diff(1)<-M_PI) z_diff(1) += PI_X2;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  S = S + R;


  /*****************************************************************************
   * Update State
   ****************************************************************************/
  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (i=0; i<n_sigma_points_; i++)
  {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= PI_X2;
    while (z_diff(1)<-M_PI) z_diff(1) += PI_X2;
    while (x_diff(3)> M_PI) x_diff(3) -= PI_X2;
    while (x_diff(3)<-M_PI) x_diff(3) += PI_X2;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // Residual
  VectorXd z_diff = z - z_pred;

  // Angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= PI_X2;
  while (z_diff(1)<-M_PI) z_diff(1) += PI_X2;

  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

//CALCULATE_NIS
#if 0
  float NIS;
  string out_file_name;

  if (n_z == n_z_radar_)
  {
    out_file_name = "NIS_Radar.txt";
  }
  else
  {
    out_file_name = "NIS_Lidar.txt";
  }

  NIS = z_diff.transpose() * S.inverse() * z_diff;
  ofstream out_file_(out_file_name.c_str(), ofstream::out | ofstream::app);
  out_file_ << NIS << endl;
#endif
}
