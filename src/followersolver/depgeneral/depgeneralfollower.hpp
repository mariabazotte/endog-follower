#ifndef GENERAL_FOLLOWER_SOLVER_HPP
#define GENERAL_FOLLOWER_SOLVER_HPP

#include "gurobi_c++.h"
#include "../../input/input.hpp"
#include "../../instance/instance.hpp"
#include "../followersolver.hpp"

class DepGeneralFollowerSolver : public AbstractFollowerSolver, public GRBCallback {
    protected:
        int pr;
        bool metropolis_hastings; 

        double rho_init = 1.0;
        double rho_mult = 2.0;
        double tol_bi_feas = 1e-4;
        double tol_part_min = 1e-4;
        double tol_part_min_obj = 1e-3;
        int rho_max_it = 30;

        int nb_max_it = 10000;
        int max_delta = 10;
        int num_pairs = 10;

        std::vector<GRBVar*> y;
        std::vector<GRBVar*> b_alpha_min;
        std::vector<GRBVar*> b_alpha_max;
        std::vector<GRBVar*> bound_alpha_min;
        std::vector<GRBVar*> bound_alpha_max;
        
        // GRBVar * alpha = NULL;
        GRBVar * alpha_min = NULL;
        GRBVar * alpha_max = NULL;
        
        GRBVar * q = NULL;
        GRBVar * w = NULL;
        GRBVar * w_delta = NULL;

        // double * alpha_ = NULL;
        double * alpha_min_ = NULL;
        double * alpha_max_ = NULL;

        std::vector<double*> y_;
        int * b_alpha_min_ = NULL;
        int * b_alpha_max_ = NULL;
        double * bound_alpha_min_ = NULL;
        double * bound_alpha_max_ = NULL;
        
        bool grb_minmax;
        bool big_m;
        double numeric;

        GRBModel * model_follower = NULL;
        GRBVar * y_follower = NULL;
        std::vector<GRBConstr> follower_constrs;

        GRBVar theta;
        bool cuts_mipnode;

        void create(){
            if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel || 
               input.getMethodDepGeneral() == Input::MethodDepGeneral::Relax || 
               input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty){
                defineOptimalFollower(); // Define follower's optimality.
                defineInteriorFollower();
                if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
                    defineAuxiliarOptimisticFollower();
                    definePessimisticFollower();
                }
                defineLeaderObj();
            }else if(input.getMethodDepGeneral() == Input::MethodDepGeneral::LocalSearch){
                defineOptimalFollower(); // Define follower's optimality.
                defineOptimisticFollower();
                defineInteriorFollower();
                definePessimisticFollower();
                Fs.set(GRB_DoubleAttr_Obj,0.5); // Used for initial solution.
                Fw.set(GRB_DoubleAttr_Obj,0.5);
                defineFollowerProblem();
            }else if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Callback){
                defineOptimalFollower(); // Define follower's optimality.
                defineOptimisticFollower();
                // definePessimisticFollower();
                defineCallback();
                if(cuts_mipnode == true) defineFollowerProblem(); // Cuts at mipnode (partial fractional solutions) 
            }
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
        
        void solve_relax();
        void solve_penalty();
        void solve_sub_penalty(double);
        void get_variables(double &);
        double computeDualityError();

        void solve_local_search();
        bool verify_leader_feas(double *);
        bool upd_follower_sol(double *);
        void defineFollowerProblem();
        bool solveFollowerProblem();
        void evaluateSAA(double &, double &, double *, double *);
        
        void solve_callback();
        void defineCallback();

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
                if(b_alpha_min[s]) delete[] b_alpha_min[s];
                if(b_alpha_max[s]) delete[] b_alpha_max[s];
                if(bound_alpha_min[s]) delete[] bound_alpha_min[s];
                if(bound_alpha_max[s]) delete[] bound_alpha_max[s];
            }
            // if(alpha) delete[] alpha;
            if(alpha_min) delete[] alpha_min;
            if(alpha_max) delete[] alpha_max;
            if(q) delete[] q;
            if(w) delete[] w;
            if(w_delta) delete[] w_delta;

            // if(alpha_) delete[] alpha_;
            if(alpha_min_) delete[] alpha_min_;
            if(alpha_max_) delete[] alpha_max_;
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                if(y_[s]) delete[] y_[s];
            }
            if(b_alpha_min_) delete[] b_alpha_min_;
            if(b_alpha_max_) delete[] b_alpha_max_;
            if(bound_alpha_min_) delete[] bound_alpha_min_;
            if(bound_alpha_max_) delete[] bound_alpha_max_;

            if(y_follower) delete[] y_follower;
            if(model_follower) delete model_follower;
        }

        void evaluate(double &, double &, double &, double &);

        void computeStrongWeakInteriorSolutions();

        void upd_solution();

        void solve();

        void callback();
};

#endif