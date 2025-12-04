#ifndef DE_STRONG_WEAK_FOLLOWER_SOLVER_HPP
#define DE_STRONG_WEAK_FOLLOWER_SOLVER_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../instance/instance.hpp"
#include "../followersolver.hpp"

class DEDepStrWkFollowerSolver : public AbstractFollowerSolver {
    protected:
        GRBVar beta;
        GRBVar obj;

        void create(){
            defineOptimisticFollower();
            definePessimisticFollower();
            defineLeaderObj();
        }

        void defineLeaderObj();

    public:
        DEDepStrWkFollowerSolver(const Input & input, Instance & instance, LeaderSolver *leader) :
                                AbstractFollowerSolver(input,instance,leader) { create(); }
        
        ~DEDepStrWkFollowerSolver() {}

        double evaluate();

        void computeStrongWeakInteriorSolutions();
};

#endif