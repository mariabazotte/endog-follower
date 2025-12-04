#include "performance.hpp"

// Mean of each sample column.
Eigen::VectorXd columnMean(const Eigen::MatrixXd &samples) {
    return samples.colwise().mean();
}

// Variance of each sample column.
Eigen::VectorXd columnVariance(const Eigen::MatrixXd &samples) {
    Eigen::VectorXd mean = columnMean(samples);
    Eigen::MatrixXd centered = samples.rowwise() - mean.transpose();
    return (centered.array().square().colwise().sum() / (samples.rows() - 1)).transpose();
}

// Autocorrelation for a vector up to maxLag.
Eigen::VectorXd autocorrelation(const Eigen::VectorXd &x, int maxLag) {
    int N = x.size();
    Eigen::VectorXd ac(maxLag + 1);
    double mean = x.mean();
    double var = (x.array() - mean).square().sum() / (N - 1); // unbiased

    for (int lag = 0; lag <= maxLag; lag++) {
        double c = 0.0;
        for (int i = 0; i < N - lag; i++)
            c += (x[i] - mean) * (x[i + lag] - mean);
        ac(lag) = c / ((N - 1) * var); // unbiased normalization
    }
    return ac;
}

// Autocorrelation time (sum until first negative)
double autocorrelationTimeMean(const Eigen::VectorXd &x, int maxLag) {
    int N = x.size();
    if (maxLag <= 0) maxLag = N - 1;

    Eigen::VectorXd ac = autocorrelation(x, maxLag);
    double tau = 0.0;
    for (int k = 1; k <= maxLag; k++) {
        if (ac[k] <= 0) break;
        tau += ac[k];
    }
    return 1.0 + 2.0 * tau;
}

// Autocorrelation time (Geyer IPS / bulk)
double autocorrelationTimeBulk(const Eigen::VectorXd &x, int maxLag) {
    int N = x.size();
    if (maxLag <= 0) maxLag = N - 1;

    Eigen::VectorXd ac = autocorrelation(x, maxLag);

    double tau = 1.0;
    for (int k = 1; k + 1 <= maxLag; k += 2) {
        double pair_sum = ac[k] + ac[k + 1];
        if (pair_sum <= 0) break;
        tau += 2.0 * pair_sum;
    }

    return tau;
}

// Univariate ESS per sample dimension.
Eigen::VectorXd effectiveSampleSizeMean(const Eigen::MatrixXd &samples) {
    int d = samples.cols();
    int N = samples.rows();
    Eigen::VectorXd ess(d);

    for (int j = 0; j < d; j++) {
        double tau = autocorrelationTimeMean(samples.col(j));
        ess(j) = N / tau;
    }
    return ess;
}

// Univariate ESS per column (bulk)
Eigen::VectorXd effectiveSampleSizeBulk(const Eigen::MatrixXd &samples) {
    int d = samples.cols();
    int N = samples.rows();
    Eigen::VectorXd ess(d);

    for (int j = 0; j < d; j++) {
        double tau = autocorrelationTimeBulk(samples.col(j));
        ess(j) = N / tau;
    }
    return ess;
}

// Multivariate ESS 
double multivariateESSMean(const Eigen::MatrixXd &samples) {
    int N = samples.rows();
    int d = samples.cols();

    double tau_sum = 0.0;
    for (int j = 0; j < d; j++) {
        double tau_j = autocorrelationTimeMean(samples.col(j)); // old ESS per dimension
        tau_sum += tau_j;
    }

    double tau_avg = tau_sum / d;

    return N / tau_avg; // approximate multivariate ESS
}

// Multivariate ESS (bulk)
double multivariateESSBulk(const Eigen::MatrixXd &samples) {
    int N = samples.rows();
    int d = samples.cols();

    // Compute univariate bulk autocorrelation time per dimension
    double tau_sum = 0.0;
    for (int j = 0; j < d; j++) {
        double tau_j = autocorrelationTimeBulk(samples.col(j));
        tau_sum += tau_j;
    }

    double tau_avg = tau_sum / d;

    return N / tau_avg; // approximate multivariate ESS
}

// Per-dimension MCSE
Eigen::VectorXd mcsePerDimension(const Eigen::MatrixXd &samples, int maxLag) {
    Eigen::VectorXd essVec = effectiveSampleSizeMean(samples);
    Eigen::VectorXd stdDev = columnVariance(samples).array().sqrt();
    return stdDev.array() / essVec.array().sqrt();
}

// Multivariate MCSE (covariance-based)
Eigen::VectorXd mcseMultivariate(const Eigen::MatrixXd &samples) {
    double essMulti = multivariateESSMean(samples);
    Eigen::MatrixXd cov = (samples.rowwise() - columnMean(samples).transpose()).transpose() *
                   (samples.rowwise() - columnMean(samples).transpose()) / (samples.rows() - 1);
    Eigen::MatrixXd covMean = cov / essMulti;
    return covMean.diagonal().array().sqrt();
}






// Autocorrelation of a vector up to maxLag
Eigen::VectorXd autocorrelationNew(const Eigen::VectorXd &x, int maxLag) {
    int N = x.size();
    if (maxLag < 0) maxLag = N - 1;
    Eigen::VectorXd ac(maxLag + 1);
    double mean = x.mean();
    double var = (x.array() - mean).square().sum() / N;

    for (int lag = 0; lag <= maxLag; lag++) {
        double c = 0.0;
        for (int i = 0; i < N - lag; i++) {
            c += (x[i] - mean) * (x[i + lag] - mean);
        }
        ac(lag) = c / ((N - lag) * var);
    }
    return ac;
}


// Compute tauHat using Geyer’s initial monotone sequence
double computeTauHat(const std::vector<Eigen::VectorXd> &autocorrs, const Eigen::VectorXd &withinVar, double varEst) {
    int numChains = autocorrs.size();
    int numDraws = autocorrs[0].size();
    std::vector<double> rhoHat;

    double rhoEven = 0, rhoOdd = 0;

    for (int t = 0; t < numDraws / 2; t++) {
        double acov = 0.0;
        for (int m = 0; m < numChains; m++) {
            acov += (numDraws - 1) * withinVar(m) * autocorrs[m](2 * t);
        }
        acov /= numChains * numDraws;
        rhoEven = (t == 0 ? 1.0 : 1.0 - (withinVar.mean() - acov) / varEst);

        acov = 0.0;
        for (int m = 0; m < numChains; m++) {
            acov += (numDraws - 1) * withinVar(m) * autocorrs[m](2 * t + 1);
        }
        acov /= numChains * numDraws;
        rhoOdd = 1.0 - (withinVar.mean() - acov) / varEst;

        if (rhoEven + rhoOdd <= 0) break;

        rhoHat.push_back(rhoEven);
        rhoHat.push_back(rhoOdd);
    }

    if (rhoEven > 0) rhoHat.push_back(rhoEven);

    // Apply initial monotone sequence
    for (size_t t = 1; t < (rhoHat.size() - 2) / 2; t++) {
        if (rhoHat[2 * t] + rhoHat[2 * t + 1] > rhoHat[2 * t - 2] + rhoHat[2 * t - 1]) {
            rhoHat[2 * t] = (rhoHat[2 * t - 2] + rhoHat[2 * t - 1]) / 2;
            rhoHat[2 * t + 1] = rhoHat[2 * t];
        }
    }

    double tauHat = -1.0;
    for (size_t t = 0; t < rhoHat.size() / 2; t++) {
        tauHat += 2.0 * (rhoHat[2 * t] + rhoHat[2 * t + 1]);
    }
    if (rhoHat.size() % 2) tauHat += rhoHat.back();

    return tauHat;
}

// Univariate ESS per dimension with multiple chains
double effectiveSampleSize(const std::vector<Eigen::MatrixXd> &chains, int dim, int maxLag) {
    int numChains = chains.size();
    if (numChains == 0) throw std::invalid_argument("No chains for ESS.");

    int numDraws = chains[0].rows();
    if (numDraws == 0) throw std::invalid_argument("Chains have zero draws.");

    // 1. Compute intra-chain means and variances
    Eigen::VectorXd intraMean(numChains);
    Eigen::VectorXd intraVar(numChains);
    std::vector<Eigen::VectorXd> autocorrs(numChains);

    double interChainMean = 0.0;

    for (int m = 0; m < numChains; m++) {
        intraMean(m) = chains[m].col(dim).mean();
        interChainMean += intraMean(m);
        intraVar(m) = ((chains[m].col(dim).array() - intraMean(m)).square().sum()) / (numDraws - 1);
        autocorrs[m] = autocorrelationNew(chains[m].col(dim), maxLag);
    }
    interChainMean /= numChains;

    // 2. Between-chain variance
    double betweenVar = 0.0;
    if (numChains > 1) {
        for (int m = 0; m < numChains; m++)
            betweenVar += std::pow(intraMean(m) - interChainMean, 2);
        betweenVar *= numDraws / (numChains - 1);
    }

    // 3. Within-chain variance
    double withinVar = intraVar.mean();

    // 4. Variance estimate
    double varEst = ((numDraws - 1) * withinVar + betweenVar) / numDraws;

    // 5. Compute tauHat
    double tauHat = computeTauHat(autocorrs, intraVar, varEst);

    // 6. ESS
    return std::min(numChains * numDraws / tauHat, numChains * numDraws * std::log10(numChains * numDraws));
}

// Compute ESS for all dimensions
Eigen::VectorXd effectiveSampleSize(const std::vector<Eigen::MatrixXd> &chains, int maxLag) {
    int dims = chains[0].cols();
    Eigen::VectorXd ess(dims);
    for (int d = 0; d < dims; d++)
        ess(d) = effectiveSampleSize(chains, d, maxLag);
    return ess;
}