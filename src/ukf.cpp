#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

 // Initializes Unscented Kalman filter
UKF::UKF() {
  // number of (aug) states
  n_x_ = 5; // num of states
  n_aug_ = 7; // number of states plus independent noises
  use_laser_ = true; // false - laser measurements will be ignored (except during init)
  use_radar_ = true; // false - radar measurements will be ignored (except during init)

  x_ = VectorXd(n_x_); // initial state vector, x = VectorXd(5)
  P_ = MatrixXd(n_x_, n_x_);//   // initial covariance matrix, P_ = MatrixXd(5, 5);

  std_a_ = 2; // Process noise standard deviation longitudinal acceleration in m/s^2
  std_yawdd_ = 2; // Process noise standard deviation yaw acceleration in rad/s^2
  std_laspx_ = 0.15;   // Laser measurement noise standard deviation position1 in m
  std_laspy_ = 0.15;   // Laser measurement noise standard deviation position2 in m
  std_radr_ = 0.3;   // Radar measurement noise standard deviation radius in m
  std_radphi_ = 0.03; // Radar measurement noise standard deviation angle in rad
  std_radrd_ = 0.3;   // Radar measurement noise standard deviation radius change in m/s

  P_ << 1, 0, 0,  0,   0,
        0, 1, 0,  0,   0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0,  0, 1;

  //Xsig       = MatrixXd(n_x_, 2 * n_x_ + 1); // sigma point of x
  Xsig_aug   = MatrixXd(n_aug_, 2 * n_aug_ + 1); //sigma point of x augment
  Xsig_pred_ = MatrixXd(n_x_,2* n_aug_+1); // predicted sigma points x = f(x, nu)
  weights_   = VectorXd(2*n_aug_+1);

  lambda_    = 3 - n_aug_;
  NIS_radar_ = 0;
  NIS_laser_ = 0;
  Num_radar_meas=0;
  Num_radar_meas_outlines=0;
  Num_laser_meas=0;
  Num_radar_meas_outlines=0;
}

UKF::~UKF() {}

// Generate Augment Sigma Points
void UKF::AugmentedSigmaPoints() {

    VectorXd x_aug = VectorXd(n_aug_); // augment vector to x
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); // augment matrix to P

    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_,n_x_) = P_;
    P_aug(n_x_,n_x_) = std_a_*std_a_;
    P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

    x_aug.head(n_x_) = x_;

    for (int i=n_x_; i<n_aug_; i++)
    {
      x_aug(i) = 0;
    }
    MatrixXd L = P_aug.llt().matrixL();     //create square root matrix
    Xsig_aug.col(0)  = x_aug;     //set first column of sigma point matrix


  //set remaining sigma points
    for (int i = 0; i < n_aug_; i++)
      {
        Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
      }
}

void UKF::SigmaPointPrediction(double delta_t) {

  //predict sigma points
  //cout<<"\ndt="<<delta_t<<endl;
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
      //extract values for better readability
      double p_x  = Xsig_aug(0,i);
      double p_y  = Xsig_aug(1,i);
      double v    = Xsig_aug(2,i);
      double yaw  = Xsig_aug(3,i);
      double yawd = Xsig_aug(4,i);
      double nu_a = Xsig_aug(5,i);
      double nu_yawdd = Xsig_aug(6,i);

      //predicted state values
      double px_p, py_p;

      //avoid division by zero
      if (fabs(yawd) > 0.001) {
          px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
          py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
      }
      else {
          px_p = p_x + v*delta_t*cos(yaw);
          py_p = p_y + v*delta_t*sin(yaw);
      }

      double v_p    = v;
      double yaw_p  = yaw + yawd*delta_t;
      double yawd_p = yawd;

      //add noise
      px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
      py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
      v_p = v_p + nu_a*delta_t;

      yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
      yawd_p = yawd_p + nu_yawdd*delta_t;

      //write predicted sigma point into right column
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
    }
}


void UKF::Prediction(){

  // set weights
   weights_(0) = lambda_/(lambda_+n_aug_);
   for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
     weights_(i) = 0.5/(n_aug_+lambda_);
   }

  //predicted state mean
  x_.fill(0.0); // Reset x_ to 0;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0); // Reset P_ to 0;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    VectorXd x_diff = Xsig_pred_.col(i) - x_;     // state difference
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI; // state psi
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI; // state psi

    while (x_(3)> M_PI) x_(3)-=2.*M_PI; // state psi
    while (x_(3)<-M_PI) x_(3)+=2.*M_PI; // state psi

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

  }
}

// Updates the state and the state covariance matrix using a laser measurement.
//void UKF::UpdateLidar(MeasurementPackage meas_package)
void UKF::UpdateLidar(const VectorXd &z) {
      //MatrixXd z = meas_package.raw_measurements_;
      MatrixXd H_laser =  MatrixXd(2, n_x_);
      H_laser << 1,0,0,0,0,
                  0,1,0,0,0;

      MatrixXd R_laser = MatrixXd(2,2);
      R_laser  << 0.0225, 0,
                 0, 0.0225;
      /*R_laser << 0, 0,
                 0, 0;*/

      MatrixXd S = H_laser*P_*H_laser.transpose() + R_laser;

      VectorXd z_pred = VectorXd(2);
      z_pred = H_laser*x_;

      VectorXd z_diff = z - z_pred;

      MatrixXd K = P_*H_laser.transpose()*S.inverse();
      x_ = x_ + K * z_diff;

      while (x_(3)> M_PI) x_(3)-=2.*M_PI;
      while (x_(3)<-M_PI) x_(3)+=2.*M_PI;

      //P_ = P_ - K*S*K.transpose();
      P_ = P_ - K*H_laser*P_;
      cout<<"Laser P=\n"<<P_<<endl;

      float error = z_diff.transpose()*S.inverse()*z_diff;
      Num_laser_meas += 1;
      if (Num_laser_meas>=10 && error>5.99) {
        Num_laser_meas_outlines += 1;
        NIS_laser_ = Num_laser_meas_outlines/Num_laser_meas;
      }
}

// Updates the state and the state covariance matrix using a radar measurement.
void UKF::UpdateRadar(const VectorXd &z) {
     int n_z = z.size(); //=3
     MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
     //transform sigma points into measurement space
     for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

       // extract values for better readibility
       double p_x  = Xsig_pred_(0,i);
       double p_y  = Xsig_pred_(1,i);
       double v    = Xsig_pred_(2,i);
       double yaw  = Xsig_pred_(3,i);
       double yawd = Xsig_pred_(4,i);
       double v1 = cos(yaw)*v;
       double v2 = sin(yaw)*v;

       // measurement model
       Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y); //r
       Zsig(1,i) = atan2(p_y, p_x);
       Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);
       //Zsig(1,i) = yaw; //psi

       float r = Zsig(0,i);
       if (r<0.01) {
         r=0.01;
         Zsig(1,i) = atan2(p_y, r)*p_x/(fabs(p_x)+0.00001);                              //phi
         Zsig(2,i) = (p_x*v1 + p_y*v2 ) / r; }  //r_dot

       /*if (r < 0.01) {
          r=0.01;
          Zsig(1,i) = atan(p_y/r)*p_x/(fabs(p_x)+0.0001);
          Zsig(2,i)  = (p_x*v1+p_y*v2)/r;}*/
     }

     //mean predicted measurement
     VectorXd z_pred = VectorXd(n_z);
     z_pred.fill(0.0);
     for (int i=0; i < 2*n_aug_+1; i++) {
         z_pred = z_pred + weights_(i) * Zsig.col(i);
     }
     //measurement covariance matrix S
     MatrixXd S = MatrixXd(n_z,n_z);
     S.fill(0.0);
     for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
       //residual
       VectorXd z_diff = Zsig.col(i) - z_pred;

       //angle normalization
       while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
       while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

       S = S + weights_(i) * z_diff * z_diff.transpose();
     }

     //add measurement noise covariance matrix
      MatrixXd R_radar = MatrixXd(n_z,n_z);
      R_radar <<    std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;
      S = S + R_radar;

      MatrixXd Tc = MatrixXd(n_x_, n_z);
      //calculate cross correlation matrix
      Tc.fill(0.0);
      for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
         //residual
         VectorXd z_diff = Zsig.col(i) - z_pred;
         //angle normalization
         while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
         while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

         // state difference
         VectorXd x_diff = Xsig_pred_.col(i) - x_;
         //angle normalization
         while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
         while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

         Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
       }

       //Kalman gain K;
       MatrixXd K = Tc * S.inverse();

       //residual
       VectorXd z_diff = z - z_pred;

       //angle normalization
       while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
       while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

       //update state mean and covariance matrix
       x_ = x_ + K * z_diff;

       while (x_(3)> M_PI) x_(3)-=2.*M_PI;
       while (x_(3)<-M_PI) x_(3)+=2.*M_PI;

       P_ = P_ - K*S*K.transpose();
       cout<<"Radar P_=\n"<<P_<<endl;

       // NIS: normalization innovation Squared
       float error = z_diff.transpose()*S.inverse()*z_diff;
       Num_radar_meas += 1;
       if (Num_radar_meas>=10 && error>7.8) {
         Num_radar_meas_outlines += 1;
         NIS_radar_ = Num_radar_meas_outlines/(Num_radar_meas-9);
       }

}

void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {

    if (!is_initialized_) {
      x_ = VectorXd(n_x_);
      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
         //cout << "Sensor_Type: RADAR" << '\n'<<endl;
        float ro = measurement_pack.raw_measurements_[0];
        float phi = measurement_pack.raw_measurements_[1];
        //float dot_psi = measurement_pack.raw_measurements_[2];
        x_ << ro*cos(phi), ro*sin(phi), 0, 0, 0;
        previous_timestamp_ = measurement_pack.timestamp_;

      }
      else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        //cout << "Sensor_Type: LASER" << '\n'<<endl;
        x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
        previous_timestamp_ = measurement_pack.timestamp_;
      }

      // done initializing, no need to predict or update
      is_initialized_ = true;
      return;
    }

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    //dt=0.05;
    cout<<"\ndt="<<dt<<endl;
    cout<< "previous time="<<previous_timestamp_/1000000.0<<"current time="<<measurement_pack.timestamp_ /1000000.0<<"\n"<<endl;
    previous_timestamp_ = measurement_pack.timestamp_;

    UKF::AugmentedSigmaPoints();
    UKF::SigmaPointPrediction(dt);
    UKF::Prediction();

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        previous_timestamp_ = measurement_pack.timestamp_;
        UKF::UpdateRadar(measurement_pack.raw_measurements_);
        }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        previous_timestamp_ = measurement_pack.timestamp_;
        UKF::UpdateLidar(measurement_pack.raw_measurements_);
      }
}
