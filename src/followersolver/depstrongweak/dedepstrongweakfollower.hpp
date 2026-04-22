#ifndef DE_STRONG_WEAK_FOLLOWER_SOLVER_HPP
#define DE_STRONG_WEAK_FOLLOWER_SOLVER_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../instance/instance.hpp"
#include "../followersolver.hpp"
#include <julia.h>

class DEDepStrWkFollowerSolver : public AbstractFollowerSolver {
    protected:
        GRBVar beta;
        GRBVar obj;
        double time_pwl;

        void create(){
            defineOptimalFollower(); // Define follower's optimality.
            defineOptimisticFollower();
            definePessimisticFollower();
            defineLeaderObj();
        }

        void defineLeaderObj();
        void computePointsPWL();
        void fixedPointsPWL();

    public:
        DEDepStrWkFollowerSolver(const Input & input, Instance & instance, LeaderSolver *leader) :
                                AbstractFollowerSolver(input,instance,leader), time_pwl(0.0) { create(); }
        
        ~DEDepStrWkFollowerSolver() {}

        void evaluate(double &, double &, double &, double &);

        void computeStrongWeakInteriorSolutions();

        void solve();
};

#endif