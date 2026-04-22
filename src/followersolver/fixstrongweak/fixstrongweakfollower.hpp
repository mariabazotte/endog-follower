#ifndef FIX_STRONG_WEAK_FOLLOWER_SOLVER_FACTORY_HPP
#define FIX_STRONG_WEAK_FOLLOWER_SOLVER_FACTORY_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../instance/instance.hpp"
#include "../followersolver.hpp"

class FixStrongWeakFollowerSolver : public AbstractFollowerSolver {
    protected:
        void create(){
            defineOptimalFollower(); // Define follower's optimality.
            if(input.getFixCoopLevel() >= 0.0001) defineOptimisticFollower();
            if(input.getFixCoopLevel() <= 0.9999) definePessimisticFollower();
            defineLeaderObj();
        }

        void defineLeaderObj();

    public:
        FixStrongWeakFollowerSolver(const Input & input, Instance & instance, LeaderSolver *leader) : 
                                    AbstractFollowerSolver(input,instance,leader) { create(); }

        ~FixStrongWeakFollowerSolver() {}

        void evaluate(double &, double &, double &, double &);

        void computeStrongWeakInteriorSolutions();
};

#endif
