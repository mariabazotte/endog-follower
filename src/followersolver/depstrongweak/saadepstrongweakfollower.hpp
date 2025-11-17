#ifndef SAA_STRONG_WEAK_FOLLOWER_SOLVER_HPP
#define SAA_STRONG_WEAK_FOLLOWER_SOLVER_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../instance/instance.hpp"
#include "../followersolver.hpp"

class SAADepStrWkFollowerSolver : public AbstractFollowerSolver {
    protected:
        int pr;

        GRBVar * z = NULL;

        void create(){
            defineOptFollower();
            definePesFollower();
            defineLeaderObj();
        }

        void defineLeaderObj();

    public:
        SAADepStrWkFollowerSolver(const Input & input,const Instance & instance, LeaderSolver *leader, int pr) :
                                AbstractFollowerSolver(input,instance,leader), pr(pr) {}
        
        ~SAADepStrWkFollowerSolver() { delete[] z; }

        double evaluate();

        void computeStrongWeakSolutions();
};

#endif