#ifndef PERFORMANCE_HPP
#define PERFORMANCE_HPP

#include <Eigen/Dense>

// Basic Statistics
Eigen::VectorXd columnMean(const Eigen::MatrixXd &);
Eigen::VectorXd columnVariance(const Eigen::MatrixXd &);

// Autocorrelation
Eigen::VectorXd autocorrelation(const Eigen::VectorXd &x);
double autocorrelationTimeMean(const Eigen::VectorXd & x, int maxLag = -1);
double autocorrelationTimeBulk(const Eigen::VectorXd & x, int maxLag = -1);

// Effective Sample Size
Eigen::VectorXd effectiveSampleSizeMean(const Eigen::MatrixXd &x);
Eigen::VectorXd effectiveSampleSizeBulk(const Eigen::MatrixXd &x);

double multivariateESSMean(const Eigen::MatrixXd &);
double multivariateESSBulk(const Eigen::MatrixXd &);

// Markov Chain Standard Error
Eigen::VectorXd mcsePerDimension(const Eigen::MatrixXd &x, int maxLag = 100);
Eigen::VectorXd mcseMultivariate(const Eigen::MatrixXd &);


Eigen::VectorXd autocorrelationNew(const Eigen::VectorXd &x, int maxLag);
double computeTauHat(const std::vector<Eigen::VectorXd> &autocorrs, const Eigen::VectorXd &withinVar, double varEst);
double effectiveSampleSize(const std::vector<Eigen::MatrixXd> &chains, int dim, int maxLag = -1);
Eigen::VectorXd effectiveSampleSize(const std::vector<Eigen::MatrixXd> &chains, int maxLag = -1);

#endif 
