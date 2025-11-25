#ifndef GENERAL_FOLLOWER_SOLVER_HPP
#define GENERAL_FOLLOWER_SOLVER_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../instance/instance.hpp"
#include "../followersolver.hpp"

class DepGeneralFollowerSolver : public AbstractFollowerSolver {
    protected:
        int pr; 

        std::vector<GRBVar*> y;
        GRBVar * alpha = NULL;
        GRBVar * alpha_min = NULL;
        GRBVar * alpha_max = NULL;
        GRBVar * q = NULL;
        GRBVar * w = NULL;

        void create(){
            defineInteriorFollower();
            defineLeaderObj();
        }

        void defineLeaderObj();

        // Define Markov Chain Monte Carlo variables, constraints and objective.
        void defineMCMCVars();
        void defineMCMCConstrs();
        void defineMCMCObj();

        // Define step (alpha) variables explicitly.
        void defineExplicitStep();

        // Define step (alpha) variables with KKT.
        void defineKKTStep();
        void defineAlphaPrimalConstrs(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);
        void defineAlphaDualVars(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);
        void defineAlphaDualConstrs(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);
        void defineAlphaCompConstrs(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,
                                    std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);

    public:
        DepGeneralFollowerSolver(const Input &,const Instance &, LeaderSolver *, int);
        
        ~DepGeneralFollowerSolver() {}

        double evaluate();

        void computeStrongWeakInteriorSolutions();
};

#endif