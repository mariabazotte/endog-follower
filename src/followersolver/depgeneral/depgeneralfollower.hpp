#ifndef GENERAL_FOLLOWER_SOLVER_HPP
#define GENERAL_FOLLOWER_SOLVER_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../instance/instance.hpp"
#include "../followersolver.hpp"

class DepGeneralFollowerSolver : public AbstractFollowerSolver, public GRBCallback  {
    protected:
        int pr;
        bool metropolis_hastings; 

        std::vector<GRBVar*> y;
        std::vector<GRBVar*> b_alpha_min;
        std::vector<GRBVar*> b_alpha_max;
        std::vector<GRBVar*> bound_alpha_min;
        std::vector<GRBVar*> bound_alpha_max;
        
        GRBVar * alpha = NULL;
        GRBVar * alpha_min = NULL;
        GRBVar * alpha_max = NULL;
        GRBVar * q = NULL;
        GRBVar * w = NULL;
        GRBVar * w_delta = NULL;

        double * x_;
        double * alpha_ = NULL;
        double * alpha_min_ = NULL;
        double * alpha_max_ = NULL;

        std::vector<double*> y_;
        std::vector<double*> bound_alpha_min_;
        std::vector<double*> bound_alpha_max_;
        int * b_alpha_min_ = NULL;
        int * b_alpha_max_ = NULL;
        
        double numeric = 1;

        bool relaxation = false;
        bool lazy_callback = true;
        bool bigm = true;

        void create(){
            defineInteriorFollower();
            if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
                defineAuxiliarOptimisticFollower();
                definePessimisticFollower();
            }
            defineLeaderObj();
        }

        void defineLeaderObj();

        // Define Markov Chain Monte Carlo variables, constraints and objective.
        void defineMCMCVars();
        void defineMCMCConstrs();
        void defineMCMCObj();

        // Define interval coefficients.
        void defineIntCoeffConstrs();

        // Define sequence.
        void defineSequence();
        void defineSequenceWithMetropolisHastings();
        void defineSequenceWithSliceSampling();

        // Define step (alpha) variables explicitly.
        void defineExplicitStep();
        void defineNewExplicitStep();
        void defineExplicitStepMaxMinConstr();

        // Define step (alpha) variables with KKT.
        void defineKKTStep();
        void defineAlphaPrimalConstrs(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);
        void defineAlphaDualVars(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);
        void defineAlphaDualConstrs(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);
        void defineAlphaCompConstrs(std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,
                                    std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &,std::vector<GRBVar*> &);
        
        void getExpression(GRBLinExpr & expr, int s, int min_or_max, int it);

        int getNbConstrsMin(int s) const {
            int nb_constrs_min = instance.getGeneralNbConstrsAlpha(pr,s,0) + 
                                 instance.getGeneralNbBdlbsAlpha(pr,s,0) + 
                                 instance.getGeneralNbBdubsAlpha(pr,s,0) + 1;
            if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
                if(metropolis_hastings == false) nb_constrs_min += 1;
            }
            return nb_constrs_min;
        }

        int getNbConstrsMax(int s) const {
            int nb_constrs_max = instance.getGeneralNbConstrsAlpha(pr,s,1) + 
                                 instance.getGeneralNbBdlbsAlpha(pr,s,1) + 
                                 instance.getGeneralNbBdubsAlpha(pr,s,1) + 1;
            if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
                if(metropolis_hastings == false) nb_constrs_max += 1;
            }
            return nb_constrs_max;
        }

    public:
        DepGeneralFollowerSolver(const Input &, Instance &, LeaderSolver *, int);
        
        ~DepGeneralFollowerSolver() {
            for(int s = 0; s < instance.getNbScenarios(); ++s) {
                if(y[s]) delete[] y[s];
                if(bound_alpha_min[s]) delete[] bound_alpha_min[s];
                if(bound_alpha_max[s]) delete[] bound_alpha_max[s];
            }
            if(alpha) delete[] alpha;
            if(alpha_min) delete[] alpha_min;
            if(alpha_max) delete[] alpha_max;
            if(q) delete[] q;
            if(w) delete[] w;
            if(w_delta) delete[] w_delta;

            if(x_) delete[] x_;
            if(alpha_) delete[] alpha_;
            if(alpha_min_) delete[] alpha_min_;
            if(alpha_max_) delete[] alpha_max_;
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                if(y_[s]) delete[] y_[s];
                if(bound_alpha_min_[s]) delete[] bound_alpha_min_[s];
                if(bound_alpha_max_[s]) delete[] bound_alpha_max_[s];
            }
            if(b_alpha_min_) delete[] b_alpha_min_;
            if(b_alpha_max_) delete[] b_alpha_max_;
        }

        void evaluate(double &, double &, double &, double &);

        void computeStrongWeakInteriorSolutions();

        void upd_solution();

        void callback();
};

#endif