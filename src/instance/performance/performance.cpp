#include "performance.hpp"

// Compute inter and intra mean and variance.
void computeMeanVars(const std::vector<Eigen::VectorXd> &x, 
                    Eigen::VectorXd & intra_mean, Eigen::VectorXd & intra_var, 
                    double & inter_mean, double & betw_var, double & with_var, double & est_var) {
    int numChains = x.size();
    int numDraws = x[0].rows();

    intra_mean.resize(numChains);
    intra_var.resize(numChains);

    for (int m = 0; m < numChains; m++) {
        intra_mean(m) = x[m].mean();
        intra_var(m) = ((x[m].array() - intra_mean(m)).square().sum()) / (numDraws - 1);
    }
                    
    // Compute between inter-chain mean.
    inter_mean = intra_mean.mean();

    // Compute between inter-chain variance.
    betw_var = 0.0;
    if (numChains > 1) {
        for (int m = 0; m < numChains; m++)
            betw_var += std::pow(intra_mean(m) - inter_mean, 2);
        betw_var *= numDraws / (numChains - 1);
    }

    // Compute within intra-chain variance.
    with_var = intra_var.mean();

    // Variance estimate.
    est_var = ((numDraws - 1) * with_var + betw_var) / numDraws;
}

// Autocorrelation of a vector up to maxLag. (multi-chain pooled)
Eigen::VectorXd autocorrelation(const std::vector<Eigen::VectorXd> & x, int maxLag) {
    int numChains = x.size();
    int N = x[0].size();
    
    Eigen::VectorXd ac(maxLag + 1);

    // Compute mean and variance.
    Eigen::VectorXd intra_mean, intra_var;
    double inter_mean, betw_var, with_var, est_var;
    computeMeanVars(x, intra_mean, intra_var, inter_mean, betw_var, with_var, est_var);

    // Case zero variance.
    if (est_var <= 0.0) {
        ac.setZero();
        ac(0) = 1.0;
        return ac;
    }

    // Compute autocorrelation.
    for (int lag = 0; lag <= maxLag; lag++) {
        double avg_c = 0.0;
        
        for (int m = 0; m < numChains; ++m) {
            const Eigen::VectorXd & chain = x[m];
            double mean_m = intra_mean(m);
            
            double c = 0.0;
            for (int i = 0; i < N - lag; i++) {
                c += (chain(i) - mean_m) * (chain(i + lag) - mean_m);
            }
            c /= (N - lag);
            avg_c += c;
        }
        avg_c /= (double)numChains; 

        ac(lag) = 1.0 - (with_var - avg_c) / est_var;
    }
    ac(0) = 1.0;

    return ac;
}

// Autocorrelation time. (sum until first negative)
double autocorrelationTimeMean(const std::vector<Eigen::VectorXd> & x, int maxLag) {
    Eigen::VectorXd ac = autocorrelation(x, maxLag);
    double tau = 0.0;
    for (int k = 1; k <= maxLag; k++) {
        if (ac[k] <= 0) break;
        tau += ac[k];
    }
    return 1.0 + 2.0 * tau;
}

// Autocorrelation time. (Geyer IPS / bulk)
double autocorrelationTimeBulk(const std::vector<Eigen::VectorXd> & x, int maxLag) {
    Eigen::VectorXd ac = autocorrelation(x, maxLag);
    double tau = 1.0;
    int k = 1;
    while (k + 1 <= maxLag) {
        double pair_sum = ac[k] + ac[k + 1];
        if (pair_sum <= 0) break;
        tau += 2.0 * pair_sum;
        k += 2;
    }
    if (k <= maxLag && ac[k] > 0) {
        tau += ac[k];  
    }

    return tau;
}

///////////////////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd autocorrelationTime(const std::vector<Eigen::MatrixXd> & chains, std::string method, int maxLag) {
    int numChains = (int)chains.size();
    int numDraws = chains[0].rows();
    int numDim = chains[0].cols();

    if(numChains == 0) throw std::invalid_argument("autocorrelationTime: No chains for computation.");
    if(numDraws == 0) throw std::invalid_argument("autocorrelationTime: Chains have zero draws.");
    if(numDim == 0) throw std::invalid_argument("autocorrelationTime: Chains have zero dimension.");
    for (const auto & chain : chains) {
        if (chain.rows() != numDraws) throw std::invalid_argument("autocorrelationTime: Inconsistent parameter draws across chains.");
        if (chain.cols() != numDim) throw std::invalid_argument("autocorrelationTime: Inconsistent parameter dimension across chains.");
    }

    if(method != "mean" && method != "bulk") 
        throw std::invalid_argument("autocorrelationTime: Method must be \"mean\" or \"bulk\".");

    if(maxLag < 0 || maxLag > numDraws - 1) maxLag = numDraws - 1;

    Eigen::VectorXd acc_time(numDim);
    for (int d = 0; d < numDim; d++){
        std::vector<Eigen::VectorXd> x;
        for(int m = 0; m < numChains; ++m) x.push_back(chains[m].col(d));

        if(method == "mean") acc_time(d) = autocorrelationTimeMean(x, maxLag);
        else if(method == "bulk") acc_time(d) = autocorrelationTimeBulk(x, maxLag);
    }

    return acc_time;
}

// Compute efficient sample size per dimension with multiple chains.
double ess(const std::vector<Eigen::VectorXd> & x, std::string method, int maxLag) {
    int numChains = x.size();
    int numDraws = x[0].size();

    if(numChains == 0) throw std::invalid_argument("ess: No chains for computation.");
    if(numDraws == 0) throw std::invalid_argument("ess: Chains have zero draws.");
    for (const auto & chain : x) {
        if (chain.size() != numDraws) throw std::invalid_argument("ess: Inconsistent parameter draws across chains.");
    }

    if(method != "mean" && method != "bulk") 
        throw std::invalid_argument("ess: Method must be \"mean\" or \"bulk\".");

    if (maxLag < 0 || maxLag > numDraws - 1) maxLag = numDraws - 1;

    // Computing tau.
    double tau = 1.0;
    if(method == "mean") tau = autocorrelationTimeMean(x, maxLag);
    else if(method == "bulk") tau = autocorrelationTimeBulk(x, maxLag);
    else throw std::invalid_argument("ess: method must be \"mean\" or \"bulk\"");
    if (!(tau > 0)) tau = 1.0;

    return (numChains * numDraws / tau);
}

// Compute efficient sample size for all dimensions.
Eigen::VectorXd ess(const std::vector<Eigen::MatrixXd> &chains, std::string method, int maxLag) {
    int numChains = chains.size();
    int numDraws = chains[0].rows();
    int numDim = chains[0].cols();

    if(numChains == 0) throw std::invalid_argument("ess: No chains for computation.");
    if(numDraws == 0) throw std::invalid_argument("ess: Chains have zero draws.");
    if(numDim == 0) throw std::invalid_argument("ess: Chains have zero dimension.");
    for (const auto & chain : chains) {
        if (chain.rows() != numDraws) throw std::invalid_argument("ess: Inconsistent parameter draws across chains.");
        if (chain.cols() != numDim) throw std::invalid_argument("ess: Inconsistent parameter dimension across chains.");
    }

    if(method != "mean" && method != "bulk") 
        throw std::invalid_argument("ess: Method must be \"mean\" or \"bulk\".");

    if (maxLag < 0 || maxLag > numDraws - 1) maxLag = numDraws - 1;

    Eigen::VectorXd vector_ess(numDim);
    for (int d = 0; d < numDim; d++){
        std::vector<Eigen::VectorXd> chains_dim;
        for(int m = 0; m < numChains; ++m) chains_dim.push_back(chains[m].col(d));

        vector_ess(d) = ess(chains_dim, method, maxLag);
    }
    return vector_ess;
}

// Compute markov chain standard error for one dimensions.
double mcse(const std::vector<Eigen::VectorXd> & x, bool var, std::string method, int maxLag) {
    int numChains = x.size();
    int numDraws = x[0].size();

    if(numChains == 0) throw std::invalid_argument("mcse: No chains for computation.");
    if(numDraws == 0) throw std::invalid_argument("mcse: Chains have zero draws.");
    for (const auto & chain_x : x) {
        if (chain_x.size() != numDraws) throw std::invalid_argument("mcse: Inconsistent parameter draws across chains.");
    }

    if(method != "mean" && method != "bulk") 
        throw std::invalid_argument("mcse: Method must be \"mean\" or \"bulk\".");

    if (maxLag < 0 || maxLag > numDraws - 1) maxLag = numDraws - 1;

    // Compute variance estimate for this dimension.
    Eigen::VectorXd intra_mean, intra_var;
    double inter_mean, betw_var, with_var, est_var;
    computeMeanVars(x, intra_mean, intra_var, inter_mean, betw_var, with_var, est_var);

    // Compute ESS using your existing function.
    double value_ess = ess(x, method, maxLag);

    // Define MCSE: if ESS <= 0 or est_var <=0, set MCSE to 0.
    if(var == true) return (value_ess >= 1e-12 && est_var >= 1e-12) ? (est_var / value_ess) : 0.0;
    else return (value_ess >= 1e-12 && est_var >= 1e-12) ? std::sqrt(est_var / value_ess) : 0.0;
}

// Compute markov chain standard error for all dimensions.
Eigen::VectorXd mcse(const std::vector<Eigen::MatrixXd> & chains, bool var, std::string method, int maxLag) {
    int numChains = (int)chains.size();
    int numDraws = chains[0].rows();
    int numDim = chains[0].cols();

    if(numChains == 0) throw std::invalid_argument("mcse: No chains for computation.");
    if(numDraws == 0) throw std::invalid_argument("mcse: Chains have zero draws.");
    if(numDim == 0) throw std::invalid_argument("mcse: Chains have zero dimension.");
    for (const auto & chain : chains) {
        if (chain.rows() != numDraws) throw std::invalid_argument("mcse: Inconsistent parameter draws across chains.");
        if (chain.cols() != numDim) throw std::invalid_argument("mcse: Inconsistent parameter dimension across chains.");
    }

    if(method != "mean" && method != "bulk") 
        throw std::invalid_argument("mcse: Method must be \"mean\" or \"bulk\".");

    if (maxLag < 0 || maxLag > numDraws - 1) maxLag = numDraws - 1;

    Eigen::VectorXd vector_mcse(numDim);
    for (int d = 0; d < numDim; d++){
        // Define auxiliar vector with column.
        std::vector<Eigen::VectorXd> chains_dim;
        for(int m = 0; m < numChains; ++m) chains_dim.push_back(chains[m].col(d));

        // Define MCSE
        vector_mcse(d) = mcse(chains_dim, var, method, maxLag);
    }
    return vector_mcse;
}

// Gelman-Rubin R-hat (per-dimension)
Eigen::VectorXd rhat(const std::vector<Eigen::MatrixXd> &chains) {
    int numChains = (int)chains.size();
    int numDraws = chains[0].rows();
    int numDim = chains[0].cols();

    if(numChains <= 1) throw std::invalid_argument("rhat Needs at least two chains for computation.");
    if(numDraws <= 1) throw std::invalid_argument("rhat: Need at least 2 draws per chain for computation.");
    if(numDim == 0) throw std::invalid_argument("rhat: Chains have zero dimension.");
    for (const auto &chain : chains) {
        if (chain.rows() != numDraws) throw std::invalid_argument("rhat: Inconsistent parameter draws across chains.");
        if (chain.cols() != numDim) throw std::invalid_argument("rhat: Inconsistent parameter dimension across chains.");
    }

    Eigen::VectorXd out(numDim);
    for (int d = 0; d < numDim; ++d) {

        std::vector<Eigen::VectorXd> chains_dim;
        for(int m = 0; m < numChains; ++m) chains_dim.push_back(chains[m].col(d));

        // Compute variance estimate for this dimension.
        Eigen::VectorXd intra_mean, intra_var;
        double inter_mean, betw_var, with_var, est_var;
        computeMeanVars(chains_dim, intra_mean, intra_var, inter_mean, betw_var, with_var, est_var);

        // If with_var is zero (all chains constant), define rhat = 1.
        if (with_var <= 0.0) {
            out(d) = 1.0;
            continue;
        }

        // Prevent negative or zero.
        if (est_var <= 0.0) {
            out(d) = 1.0;
            continue;
        }

        out(d) = std::sqrt(est_var/with_var);
    }
    return out;
}