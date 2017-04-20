#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  //VectorXd x_aug;

  ///* state covariance matrix
  MatrixXd P_;
  //MatrixXd P_aug;

  // sigma points
  //MatrixXd Xsig;
  MatrixXd Xsig_aug;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  double previous_timestamp_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  float NIS_radar_;
  float Num_radar_meas;
  float Num_radar_meas_outlines;

  ///* the current NIS for laser
  float NIS_laser_;
  float Num_laser_meas;
  float Num_laser_meas_outlines;

  UKF(); // Constructor
  virtual ~UKF(); // Destructor

  void AugmentedSigmaPoints();
  void SigmaPointPrediction(double delta_t);
  void Prediction();
  void UpdateLidar(const VectorXd &z);
  void UpdateRadar(const VectorXd &z);
  void ProcessMeasurement(MeasurementPackage meas_package);

};

#endif /* UKF_H */
