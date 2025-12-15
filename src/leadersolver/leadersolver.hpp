#ifndef LEADER_SOLVER_HPP
#define LEADER_SOLVER_HPP

#include "gurobi_c++.h"
#include "status/status.hpp"
#include "../input/input.hpp"
#include "../instance/instance.hpp"
#include "../followersolver/followersolverfactory.hpp"

class AbstractFollowerSolver;

class AbstractLeaderSolver {
    protected:
        const Input & input;
        Instance & instance;

        std::string name;

        GRBEnv *env = NULL;

        double lb;
        double ub;
        double gap;
        double time_;
        double nb_bnb;
        int nb_sol;
        Status status;

        void writeCompFixStrongWeak(std::string &, double, double, double, double) const;
        void writeCompDepStrongWeak(std::string &, Input::TypesDepStrongWeak, double, double, double, double) const;
        void writeCompDepGeneral(std::string &, Input::TypesDepGeneral, int, double, double, double, double, double, double, double, double) const;

    public:
        AbstractLeaderSolver(const Input & input, Instance & instance,std::string name, GRBEnv * env) : input(input), instance(instance), name(name), env(env) {
            lb = -std::numeric_limits<double>::infinity();
            ub = std::numeric_limits<double>::infinity();
            gap = 1.0;
            time_ = 0.0;
            nb_bnb = 0.0;
            nb_sol = 0;
            status = Status::Not_Solved;
        } 

        virtual ~AbstractLeaderSolver() { delete env; }

        double getLB() const { return lb; }
        double getUB() const { return ub; }
        double getGAP() const { return gap; }
        double getTime() const { return time_; }
        double getNbBnB() const { return nb_bnb; }
        double getNbSol() const { return nb_sol; }
        Status getStatus() const { return status; }  

        // Method to solve problem.
        virtual void solve() = 0;

        // Getters for leader decision and follower.
        virtual double getX_(int) const = 0;
        virtual const double * getX_() const = 0;
        virtual AbstractFollowerSolver * getFollower() const = 0;

        // Methods to write solution.
        virtual std::string write() const = 0;
        virtual std::string writeHead() const = 0; 
        std::string writeComp() const;
        std::string writeHeadComp() const;
};

class LeaderSolver : public AbstractLeaderSolver {
    private:
        int pr;

        GRBModel *model = NULL;
        GRBVar *x = NULL;
        double *x_ = NULL;
    
        AbstractFollowerSolver *follower = NULL;

        void params();
        void write_inf();
        void write_sol();
        void print();
        void verify_stat();
        void upd_solution();
        void create();

    public:
        LeaderSolver(const Input &, Instance &, std::string);
        LeaderSolver(const Input &, Instance &, std::string, int, GRBEnv *);

        ~LeaderSolver();

        int getProblem() const { return pr; }

        GRBEnv * getGRBEnv() { return env; }
        GRBModel * getGRBModel() { return model; }
        GRBVar & getX(int i) { return x[i]; }

        double getX_(int i) const { return x_[i]; }
        const double * getX_() const { return x_; }
        AbstractFollowerSolver * getFollower() const { return follower; }

        void solve();  
        std::string write() const;
        std::string writeHead() const; 
};

class SAALeaderSolver : public AbstractLeaderSolver {
    private:
        std::vector<LeaderSolver*> solvers;

        double var_lb;
        double var_ub;
        double var_gap;
        double var_nb_bnb;

        double stat_lb;
        double stat_ub;
        double stat_gap;
        int nb_notopt;
        int best_problem;

        double time_lb;
        double time_eval;

        double f_mean;
        double f_var;

        std::vector<double> lb_vector;
        std::vector<double> ub_vector;
        std::vector<double> gap_vector;
        std::vector<double> nb_bnb_vector;

        void estimators();
        void testProblem(int);

    public:
        SAALeaderSolver(const Input &, Instance &, std::string);

        ~SAALeaderSolver();

        double getX_(int i) const { return solvers[best_problem]->getX_(i); }
        const double * getX_() const { return solvers[best_problem]->getX_(); }
        AbstractFollowerSolver * getFollower() const { return solvers[best_problem]->getFollower(); }

        void solve();
        std::string write() const;
        std::string writeHead() const; 
};

#endif

