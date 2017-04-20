#include <iostream>
#include "Eigen/Dense"
#include <vector>
#include "ukf.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

int main() {

	UKF ukf;
  ukf.x_ <<5.7441,
         1.3800,
         2.2049,
         0.5015,
         0.3528;

  VectorXd z = VectorXd(3);
  z << 10, 0.5, 0.05;

  ukf.P_ << 0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
cout<<"P=" << ukf.P_<<endl;
 ukf.AugmentedSigmaPoints();
 ukf.SigmaPointPrediction(0.01);
 ukf.Prediction();
 //ukf.UpdateLidar(z);
 ukf.UpdateRadar(z);

 cout<< "New P="<<ukf.P_ <<endl;

	return 0;
}
