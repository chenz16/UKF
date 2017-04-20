# Unscented Kalman Filter Project Submission
Self-Driving Car Engineer Nanodegree Program

---

## Dependencies

  * cmake >= 3.5
  * make >= 4.1
  * gcc/g++ >= 5.4

## Basic Build Instructions

  1. Clone this repo.
  
  2. Make a build directory: `mkdir build && cd build`
  
  3. Compile: `cmake .. && make` 
    * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
    
  4. Run it: `./ExtendedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
   
      - eg. `./ExtendedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../data/output3.txt`


## Go trhough project rubric 

### compile without errors 

Code compiled without errors with cmake and make.


### Accuacy 
run `./ExtendedKF ../data/obj_pose-laser-radar-synthetic-input.txt  
../data/output3.txt` to get the data accuacy

    Accuracy - RMSE:
    0.0747166
    0.088438
    0.382311
    0.308456

which is quite close to the criteria:
`[.09, .10, .40, .30].`

This result assumed the following process noise factor levels

     std_a_ = 3; // Process noise standard deviation longitudinal acceleration in m/s^2
     std_yawdd_ = 3; // Process noise standard deviation yaw acceleration in rad/s^2

To figure out if measurement and prediction of measurement from process model statistically aligns, the NIS is calculated as :

       Number of Laser Measurement =240 NIS Laser =0.0201613. 
       Number of Radar Measurement =241 NIS Radar =0.039604
       
Here the first 9 laser/radar measurement is not accounted to make sure the kalman filter prediction is stable. The NIS number is close to 5%. 

To ensure the error is exactly within the criteria: `[.09, .10, .40, .30].`. The process noises are modified as

     std_a_ = 2; // Process noise standard deviation longitudinal acceleration in m/s^2
     std_yawdd_ = 2; // Process noise standard deviation yaw acceleration in rad/s^2
     
 The NIS is calculated as:
 
     Number of Laser Measurement =240 NIS Laser =0.0282258
     Number of Radar Measurement =241 NIS Radar =0.0343348
which are still within 5%. 

The accuracy for the prediction is:

    Accuracy - RMSE:
    0.0700293
    0.085612
    0.355426
    0.253763

which falls into the criteria `[.09, .10, .40, .30].`

### Follows the Correct Algorithm

`Your Sensor Fusion algorithm follows the general processing flow as taught in the preceding lessons.
Your Kalman Filter algorithm handles the first measurements appropriately.
Your Kalman Filter algorithm first predicts then updates.
Your Kalman Filter can handle radar and lidar measurements.`

The algorithm treats the first measurement as initialization of the states (x). Normal steps of Unscented kalman filter calculation are applied from the second measurement.  

To predict the state in terms of mean and covariance in next time step, sigma points of current step states x (i.e. sampling points of state x within the defined uncertain range through covariance matrix P) are first generated. Next, the new states based on these sigma points are predicted; then the mean and covariance of the state at next time step are calculated based on these predicted points.

In the update step, we handle radar and lidar data separately given lidar measurement is linear function of state x. The update of lidar measurement follows the standard kalman filter method, which i skip explanation here. 

The radar measurement is nonlinear function of the predicted state x. we use the sigma points of predicted states obtained in last step to calculate the mean and covariance of the error between the measumrenet and prediction of measurement. The kalman gain is calculated accordingly. The the state and its covariance is updated based on predicted state, kalman gain, error between measurement and measurement prediction. 

All the key step of UKF are implemented in ukf.cpp. 

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

To check if the process model align with the measurement, i calculated the NIS of radar and lidar. I assumed the sensing standard deviation is accurate and the process noise is something i need to tune. To tune those process noise standard deviation, i am able to make the RMSE within the criteria and also make NIS within 5%. 

### Code Efficiency

`Your algorithm should avoid unnecessary calculations.`

The code follows the general approach of the sample code from the course material. 

