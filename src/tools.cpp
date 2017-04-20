#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

        int n_measure = estimations[0].size();
        VectorXd rmse = VectorXd(n_measure);
      	rmse.fill(0);

      	if(estimations.size() != ground_truth.size()
      			|| estimations.size() == 0){
      		std::cout << "Invalid estimation or ground_truth data" << std::endl;
      		return rmse;
      	}

      	//accumulate squared residuals
      	for(unsigned int i=0; i < estimations.size(); ++i){
      		VectorXd residual = estimations[i] - ground_truth[i];

      		//coefficient-wise multiplication
      		residual = residual.array()*residual.array();
      		rmse += residual;
      	}

      	//calculate the mean
      	rmse = rmse/estimations.size();

      	//calculate the squared root
      	rmse = rmse.array().sqrt();

      	//return the result
      	return rmse;
}
