#ifndef PERFORMANCE_HPP
#define PERFORMANCE_HPP

#include <Eigen/Dense>

// Autocorrelation Time.
Eigen::VectorXd autocorrelationTime(const std::vector<Eigen::MatrixXd> & x, std::string method = "bulk", int maxLag = -1);

// Effective Sample Size.
double ess(const std::vector<Eigen::VectorXd> & x, std::string method = "bulk", int maxLag = -1);
Eigen::VectorXd ess(const std::vector<Eigen::MatrixXd> &chains, std::string method = "bulk", int maxLag = -1);

// Markov Chain Standard Error.
double mcse(const std::vector<Eigen::VectorXd> & x, bool var = false, std::string method = "bulk", int maxLag = -1);
Eigen::VectorXd mcse(const std::vector<Eigen::MatrixXd> &chains, bool var = false, std::string method = "bulk", int maxLag = -1);

// R-hat.
Eigen::VectorXd rhat(const std::vector<Eigen::MatrixXd> &chains);

#endif 
