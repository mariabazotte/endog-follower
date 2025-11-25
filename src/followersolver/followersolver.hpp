#ifndef ABSTRACT_FOLLOWER_SOLVER_HPP
#define ABSTRACT_FOLLOWER_SOLVER_HPP

#include "gurobi_c++.h"
#include "../input/input.hpp"
#include "../instance/instance.hpp"

class LeaderSolver;

class AbstractFollowerSolver {
    protected:
        const Input & input;
        const Instance & instance;

        LeaderSolver *leader;        // Leader problem.

        GRBVar Fs;                   // Optimal leader optimistic value. Adapted for optimality or near-optimality cases.
        GRBVar Fw;                   // Optimal leader pessimistic value. Adapted for optimality or near-optimality cases.
        
        GRBVar fs;                   // Optimal follower value. 
        GRBVar fs_eps;               // Near optimal optimistic follower value. Corresponds to near-optimal solution chosen by the cooperative follower.
        GRBVar fw_eps;               // Near optimal pessimistic follower value. Corresponds to near-optimal solution chosen by the adversarial follower.
        GRBVar fi_eps;               // Near optimal interior follower value. Corresponds to near-optimal solution that maximizes the min slack from constraints. 

        GRBVar si;                   // Maximum min slack obtained to find an interior optimal or near-optimal solution.

        GRBVar *ys = NULL;           // Optimal follower optimistic solution.
        GRBVar *yw = NULL;           // Optimal or near-optimal follower pessimistic solution. (Formulation allows the definition of only yw). 
        GRBVar *yi = NULL;           // Optimal or near-optimal follower interior solution. (Formulation allows the definition of only yi).
        GRBVar *ys_eps = NULL;       // Near-optimal follower optimistic solution. (Formulation needs definition of both ys and ys_eps, as we need a variable computing the optimal follower value).

        double Fs_;                  // Auxiliar parameters taking the value of the variables when problem is solved.
        double Fw_;

        double fs_;
        double fs_eps_;
        double fw_eps_;
        double fi_eps_;

        double *ys_ = NULL;
        double *yw_ = NULL;
        double *yi_ = NULL;
        double *ys_eps_ = NULL;

        double eval_avg;                 // Parameters for optimal leader objective (part corresponding to follower variables) at evaluation problem.
        std::vector<double> eval_sce;    

        // Define primal feasibility.
        void definePrimalVars(GRBVar *&, std::string);
        void definePrimalConstrs(GRBVar *&, GRBVar *&, GRBVar *&, GRBVar *&, bool, std::string);
        void definePrimalConstrs(GRBVar *&, std::string);
        
        // Define dual feasibility.
        void defineDualVars(GRBVar *&, GRBVar *&, GRBVar *&, bool, bool);
        void defineDualConstrs(GRBVar *&, GRBVar *&, GRBVar *&, bool, bool);
        
        // Define complementarity constraints.
        void defineCompConstrs(GRBVar *&, GRBVar *&, GRBVar *&, GRBVar *&, GRBVar *&, GRBVar *&);
        
        // Define follower and leader objectives.
        bool fs_not_init;
        void defineFollowerOptimalObj(GRBVar *&, std::string);
        void defineFollowerNearOptimalObj(GRBVar &, GRBVar *&, GRBVar &, std::string);
        void defineLeaderObj(GRBVar &, GRBVar *&, std::string);

        // Define optimal follower.
        void defineOptimalFollower();

        // Define optimistic, pessimistic, and interior followers solutions.
        void defineOptimisticFollower();
        void definePessimisticFollower();
        void defineInteriorFollower();

        // Abstract methods according to the specific follower behavior.
        virtual void create() = 0;
        virtual void defineLeaderObj() = 0;

        // Compute optimistic (strong), pessimistic (weak), and interior solutions 
        // given the optimal leader decision obtained after solving the problem.
        void computeStrongWeakInteriorSolutions(bool, bool, bool);
        
    public:
        AbstractFollowerSolver(const Input & input,const Instance & instance, LeaderSolver *leader) :
                                input(input), instance(instance), leader(leader), fs_not_init(true) {}

        virtual ~AbstractFollowerSolver() { 
            if(ys_) delete[] ys_;
            if(yw_) delete[] yw_;
            if(yi_) delete[] yi_;
            if(ys_eps_) delete[] ys_eps_;
        }

        // Getters for follower optimal objective values.
        double getFollowerOptimalObj() const { return fs_; }
        double getFollowerNearOptimalOptObj() const { return fs_eps_; }
        double getFollowerNearOptimalPesObj() const { return fw_eps_; }

        // Getters for leader optimal objectuve values.
        double getLeaderOptObj() const { return Fs_; }
        double getLeaderPesObj() const { return Fw_; }

        // Getter for follower optimistic/strong solution.
        double getYs_(int i) const { if(input.isFollowerNearOpt() == true) 
                                        return ys_eps_[i];
                                     else return ys_[i]; }
        const double * getYs_() const { if(input.isFollowerNearOpt() == true) 
                                            return ys_eps_;
                                        else return ys_; }

        // Getter for foolower pessimistic/weak solution.
        double getYw_(int i) const { return yw_[i]; }
        const double * getYw_() const { return yw_; }
        
        // Getter for follower interior solution.  
        double getYi_(int i) const { return yi_[i]; }                         
        const double * getYi_() const { return yi_; }
        
        // Getter for difference between scenario value and average value of the leader objective corresponding
        // to follower variables for the evaluation problem. Used to compute variance of the estimator.
        double getDiffEvalScenario(int s) const { return (eval_sce[s] - eval_avg); }
        
        // Function to update the optimal solution after solving problem.
        void upd_solution();
        
        // Function to evaluate the optimal leader solution obtained after solving
        // the problem. It either computes the true expected value (enumerates all scenarios)
        //  or uses the evaluation SAA problem depending on the follower behavior.
        virtual double evaluate() = 0;

        // Function to compute the 
        virtual void computeStrongWeakInteriorSolutions() = 0;
};

#endif