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
../data/output3.txt` to get the first data accuacy

Accuracy - RMSE:
0.0747166
0.088438
0.382311
0.308456

which is quite close to the criteria:
`[.09, .10, .40, .30].`

This result assumed the following process noise factor 

     std_a_ = 3; // Process noise standard deviation longitudinal acceleration in m/s^2
     std_yawdd_ = 3; // Process noise standard deviation yaw acceleration in rad/s^2

To figure out if measurement and prediction of measurement from process model statistically makes sense, the NIS is calcuated:

       Number of Laser Measurement =240 NIS Laser =0.0201613. 
       Number of Radar Measurement =241 NIS Radar =0.039604
       
Here the first 10 laser/radar measurement is not accounted to make sure the kalman filter prediction is stable. The NIS number is close to 5%. 

To ensure the error is exactly within the criteria: `[.09, .10, .40, .30].`. The process noises are modifed as

     std_a_ = 2; // Process noise standard deviation longitudinal acceleration in m/s^2
     std_yawdd_ = 2; // Process noise standard deviation yaw acceleration in rad/s^2
     
 The NIS is calcuated as:
 
     Number of Laser Measurement =240 NIS Laser =0.0282258
     Number of Radar Measurement =241 NIS Radar =0.0343348
which are still within 5%. 

The accuacy for the prediction is:

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

The algorithm treats the first measurement as initlization of the states (x). Normal steps of kalman filter calcuation are applied from the second measurement.  

The algorithm predicted through calling the prediction function in the instance of `ekf_` defined from the class   `KalmanFilter`;  the algorithm then updated the states and covariance matrix through calling the update functions in the instance of `ekf_` defined from the class `KalmanFilter`. The update functions are different for laser and radar data. When laser data was dected, the update function is a standard update function of Kalman filter; when radar data was dected, the update function used extended kalman filter method. In the extended kalman filter, the prediction of current measurement is calcated through non-linear measurement functions h(x); the Jacobian matrix was obtained from h(x) and used for calculating the kalman gain K.    

     
    ekf_.Predict();

    /*****************************************************************************
    *  Update
    ****************************************************************************/
    // Renew H, R before update
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Update measurement related matrix H and R
      cout << "Sensor_Type: RADAR" << '\n'<<endl;
      ekf_.H_ = tool.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      // Radar updates
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)  {
      // Update measurement related matrix H and R
      cout << "Sensor_Type: LASER" << '\n'<<endl;
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;

      // Laser updates
      ekf_.Update(measurement_pack.raw_measurements_);
    }



### Code Efficiency

`Your algorithm should avoid unnecessary calculations.`

The code follows the general approach of the sample code from the course material. 


