#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include "../input/input.hpp"
#include "bilevelmodel/bilevelmodel.hpp"
#include "performance/performance.hpp"

#include <list>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
namespace py = pybind11;
using namespace pybind11::literals; 

class Instance {
    private:
        const Input & input;                           /* Input. */

        BilevelModel * bilevelmodel = NULL;            /* Bilevel model with variables and constraints information. */
        
        std::string instance_file;                     /* File with instance information. */
        std::string name;                              /* Name of the instance. */ 

        int nb_saa_problems;                           /* Number of SAA problems. */
        int nb_saa_scenarios;                          /* Number of scenarios for each SAA problem. */
        int nb_saa_thinning;                           /* Number of thinning scenarios for the SAA problem. */
        
        int nb_val_problems;                           /* Number of problems (chains) for the validation test. */
        int nb_val_scenarios;                          /* Number of scenarios for the validation test. */
        int nb_val_thinning;                           /* Number of thinning scenarios for the validation test. */

        double critical_tstudent;                      /* Critical value for t-student distribution. */
        double critical_normal;                        /* Critical value for normal distribution. */

        double seed;                                   /* Seed for random data generation (scenarios). */
        std::default_random_engine rd_engine;          /* Random engine for data generation (scenarios).*/ 
        std::uniform_real_distribution<> uniform_eps;  /* Uniform distribution for sampling scenarios. (decision-dependent strong weak) */
        std::uniform_real_distribution<> uniform;      /* Uniform distribution for sampling scenarios. (decision-dependent general step and accept random variables) */
        std::normal_distribution<> normal;             /* Normal distribution for sampling scenarios. (decision-dependent general direction random vector for full hit and run) */
        std::uniform_int_distribution<> uniform_int;   /* Uniform int distribution for sampling scenarios. (decision-dependent general direction random vector for coordinate hit and run) */

        bool coordinate_mcmc;                          /* 0-> Do not use coordinate (use full instead) markov chain monte carlo hit and run, 1-> Use coordinate markov chain monte carlo hit-and-run. */
        bool constr_normal_transform;                  /* 0-> Do not use constraint normal transform for hit and run, 1-> Use constraint normal transform for hit-and-run. */
        bool polyround_transform;                      /* 0-> Do not use polyround transform for hit-and-run, 1-> Use polyround transform for hit-and-run. */
        bool polyround_transform_specific;             /* 0-> Do not use polyround transform specific (for a fixed leader decision at evaluation step only) for hit and run, 1-> Use polyround transform specific for hit-and-run. */
        bool transform_specific_executed;              /* 0-> Polyround transform specific has not been executed yet, 1-> Polyround transform specific has already been executed. */
        bool metropolis_acceptance;                    /* */

        Eigen::MatrixXd A_follower;                    /* Inequality constraints matrix on follower problem. (only follower variables)*/
        Eigen::VectorXd b_follower;                    /* RHS value inequality constraints on follower problem. */
        Eigen::MatrixXd Aeq_follower;                  /* Equality constraints matrix on follower problem. (only follower variables) */

        Eigen::MatrixXd nullspace;                     /* Null space matrix for equality constraints on follower problem. */
        Eigen::MatrixXd constrnorm_tr;                 /* Transformation matrix when constr_normal_transform = true. */
        Eigen::MatrixXd polyround_tr;                  /* Transformation matrix when polyround_transform = true, or polyround_transform_specific = true. */

        // Parameters for General decision-dependent Proportional and Strong-Fragile.
        double gen_min_probab;                         /* pi_min: minimum value of the function proportional to the decision-dependent probability. */
        double gen_max_probab;                         /* pi_max: maximum value of the function proportional to the decision-dependent probability. */
        double gen_ref_pt_probab;                      /* pi_c: reference value point of the function proportional to the decision-dependent probability. */
       
        // Scenarios for Decision-dependent strong-weak case.
        std::vector<std::vector<double>> scenarios_dep_strwk;                   /* scenarios_dep_strwk[pr][s]: inverse value of beta or cooperation level sampled from U(0,1) for transformed scenario s of the SAA problem pr. */

        // Scenarios for Decision-dependent general case.
        std::vector<std::vector<std::vector<double>>> dir_scenarios_dep_gen;    /* dir_scenarios_dep_gen[pr][s]: (upsilon) direction vector for general case transformed scenario s of the SAA problem pr. */
        std::vector<std::vector<double>> step_scenarios_dep_gen;                /* step_scenarios_dep_gen[pr][s]: (tau) sample U(0,1) to define step for general case transformed scenario s of the SAA problem pr.  */
        std::vector<std::vector<double>> accep_scenarios_dep_gen;               /* accep_scenarios_dep_gen[pr][s]: (phi) sample U(0,1) to define acceptance rate for general case transformed scenario s of the SAA problem pr. */

        std::vector<std::vector<double>> coeff_obj_alpha;                       /* coeff_obj_alpha[pr][s]: coefficient of the alpha variable at the follower near optimal objective value constraint at problem pr and scenario s. */
        std::vector<std::vector<std::vector<double>>> coeff_alpha;              /* coeff_alpha[pr][s][c]: coefficient of the alpha variable at the follower constraint c at problem pr and scenario s. */
        std::vector<std::vector<std::vector<std::vector<int>>>> constrs_alpha;  /* constrs_alpha[pr][s][0/1]: list of follower constraint defining alpha min (0), and alpha max (1) in problem pr and scenario s. */
        std::vector<std::vector<std::vector<std::vector<int>>>> bdlbs_alpha;    /* bdlbs_alpha[pr][s][0/1]: list of lower bounds on follower variables defining alpha min (0), and alpha max (1) in problem pr and scenario s. */
        std::vector<std::vector<std::vector<std::vector<int>>>> bdubs_alpha;    /* bdubs_alpha[pr][s][0/1]: list of upper bounds on follower variables defining alpha min (0), and alpha max (1) in problem pr and scenario s. */

        // Evaluation scenarios for Decision-dependent general case.
        std::vector<std::vector<double>> dir_eval_dep_gen;                      /* dir_eval_dep_gen[s]: (upsilon) direction for general case for evaluation step.  */
        std::vector<double> step_eval_dep_gen;                                  /* step_eval_dep_gen[s]: (tau) sample U(0,1) to define step for general case for evaluation step. */
        std::vector<double> accep_eval_dep_gen;                                 /* accep_eval_dep_gen[s]: (phi) sample U(0,1) to define acceptance rate for general case for evaluation step. */
        
        std::vector<double> eval_coeff_obj_alpha;                               /* eval_coeff_obj_alpha[s]: */
        std::vector<std::vector<double>> eval_coeff_alpha;                      /* eval_coeff_alpha[s][c]: */
        std::vector<std::vector<std::vector<int>>> eval_constrs_alpha;          /* eval_constrs_alpha[s][0/1]: */
        std::vector<std::vector<std::vector<int>>> eval_bdlbs_alpha;            /* eval_bdlbs_alpha[s][0/1]: list of lower bounds on follower variables defining alpha min (0), and alpha max (1) in scenario s of evaluation step.*/
        std::vector<std::vector<std::vector<int>>> eval_bdubs_alpha;            /* eval_bdubs_alpha[s][0/1]: list of upper bounds on follower variables defining alpha min (0), and alpha max (1) in scenario s of evaluation step. */

        // Problem configurations for comparison.
        std::vector<double> fix_strwk_cg = {0.0, 0.3, 0.5, 0.7, 1.0};
        std::vector<double> dep_strwk_thr_cg = {0.5, 2.0, 5.0, 10.0};
        std::vector<double> dep_strwk_str_cg = {0.5, 2.0, 5.0, 10.0};
        std::vector<double> dep_strwk_frg_cg = {0.5, 2.0, 5.0, 10.0};
        std::vector<int> dep_gen_int_cg = {4, 5};
        std::vector<double> dep_gen_scal_cg = {0.5, 2.0, 5.0, 10.0};

        void defineName();                                              /* Define instance name. */
        void defineCriticalValues();                                    /* Define critical values. */
        void defineDepGenParams();                                      /* Define auxiliar parameters for General Decision-dependent case. */

        void defineFollowerProblemEigenMatrix();                        /* Define Eigen matrix and vector representing the follower problem with only follower variables. */
        void computeDepGenConstrNormTransform();                        /* Function to compute transformation for follower problem based on constraint normal for hit-and-run used on General decision-dependent cases. */
        void computeDepGenPolyRoundTransform(const Eigen::VectorXd &);  /* Function to compute transformation for follower problem based on the polyround library for hit-and-run used on General decision-dependent cases. */
        void computeDepGenPolyRoundTransform(const double *, double);   /* Function to compute transformation for follower problem considering a fixed leader decision based on the polyround library for hit-and-run. */
        void nullSpaceFollowerEqConstrs(bool);                          /* Function to define null space for equality constraints of the follower problem. */
        void updInfoEvaluateDepGeneral(const double *, double);         /* Function to update sampling information when polyround_transform_specific = true. */ 
         
        void generateScenarios();                                       /* Generate scenarios. */
        double inverseDepStrongWeakScenario();                          /* Function to sample scenarios for decision-dependent strong weak case. */
        void directionDepGeneralScenario(std::vector<double> & dir, bool eval = false);  /* Function to sample direction scenarios for decision-dependent general case. */
        void directionDepGeneralScenarioEvalModify();                   /* Function to update direction samples when polyround_transform_specific = true. */ 
        void defineGeneralStepInfo(int,                                 /* Function to define scenario step information for decision-dependent general case. */
            std::vector<double> &, std::vector<std::vector<double>> &, 
            std::vector<std::vector<std::vector<int>>> &, 
            std::vector<std::vector<std::vector<int>>> &, 
            std::vector<std::vector<std::vector<int>>> &);  
        
        double computeFollowerObj(const double * y_) const;
        double computeLeaderObjLeaderVars(const double *) const;
        double computeLeaderObjFollowerVars(const double *) const;

        // Getter for current parameter values (F_c, pi_c, m) for evaluation general decision-dependent.
        void computeParamsDepGeneral(double,Input::TypesDepGeneral, int, double, double &, double &, double &) const;
        
        // Getters for current step/alpha min and max for evaluation general decision-dependent.
        double defineEvalMinStep(int, const double *, const double *, double) const;
        double defineEvalMaxStep(int, const double *, const double *, double) const;
        
        // Getters for sampling current step/alpha final value for evaluation general decision-dependent.
        double sampleEvalAlpha(int, double, double, double, double, double, double, Input::TypesDepGeneral) const;
        double sampleEvalUniformAlpha(int, double, double) const;
        double sampleEvalLinearAlpha(int, double, double, double, double) const;
        double sampleEvalExponentialAlpha(int, double, double, double) const;

        // Getter for general probability for current follower optimal objective value.
        double getEvalDepGeneralProb(double, double, double, double, Input::TypesDepGeneral) const;

        void performEvalMCMC(const Eigen::MatrixXd &) const;
        void evaluateDepGeneralLibrary(const double *, double, double, double, double, double, Input::TypesDepGeneral);

    public:
        Instance(const Input &);                 
        ~Instance(){ delete bilevelmodel; }

        BilevelModel * getModel() const { return bilevelmodel; }
        
        std::string getName() const { return name; }

        int getNbProblems() const { return nb_saa_problems; }
        int getNbScenarios() const { return nb_saa_scenarios; }
        int getNbValidateScenarios() const { return nb_val_scenarios; }
        int getNbScenarios(int pr) { if(pr == nb_saa_problems) return nb_val_scenarios; 
                                        else return nb_saa_scenarios; }

        double getCriticalNormal() const { return critical_normal;}
        double getCriticalTStudent() const { return critical_tstudent;}

        double getStepSizeInterval() const { return 100; }
        double getGeneralMaxProbab() const { return gen_max_probab; }
        double getGeneralMinProbab() const { return gen_min_probab; }
        double getGeneralRefProbab() const { return gen_ref_pt_probab; }
        double getRefLeaderObj() const { return (bilevelmodel->leader_lb + bilevelmodel->leader_ub)/2.0; }

        // Getters for scenario information for strong-weak decision-dependent case.
        double getStrongWeakScenario(int pr, int s) const { return scenarios_dep_strwk[pr][s]; }
        
        // Getters for scenario information for general decision-dependent case.
        double getGeneralDirection(int pr, int s, int i) const { return dir_scenarios_dep_gen[pr][s][i]; }
        double getGeneralStep(int pr, int s) const { return step_scenarios_dep_gen[pr][s]; }
        double getGeneralAccep(int pr, int s) const { return accep_scenarios_dep_gen[pr][s]; }   
        double getGeneralCoeffAlphaObj(int pr, int s) const { return coeff_obj_alpha[pr][s]; }
        double getGeneralCoeffAlpha(int pr, int s, int c) const { return coeff_alpha[pr][s][c]; }
        int getGeneralNbConstrsAlpha(int pr, int s, int m) const { return constrs_alpha[pr][s][m].size(); }
        int getGeneralConstrsAlpha(int pr, int s, int m, int j) const { return constrs_alpha[pr][s][m][j]; }
        int getGeneralNbBdlbsAlpha(int pr, int s, int m) const { return bdlbs_alpha[pr][s][m].size(); }
        int getGeneralBdlbsAlpha(int pr, int s, int m, int j) const { return bdlbs_alpha[pr][s][m][j]; }
        int getGeneralNbBdubsAlpha(int pr, int s, int m) const { return bdubs_alpha[pr][s][m].size(); }
        int getGeneralBdubsAlpha(int pr, int s, int m, int j) const { return bdubs_alpha[pr][s][m][j]; }

        // Getters for evaluation scenarios information for general decision-dependent case.
        double getGeneralDirectionEval(int s, int i) const { return dir_eval_dep_gen[s][i]; }
        double getGeneralStepEval(int s) const { return step_eval_dep_gen[s]; }
        double getGeneralAccepEval(int s) const { return accep_eval_dep_gen[s]; }

        // Getters for configurations to compare current model optimal decision.
        const std::vector<double> & getFixStrongWeakCoopConfigs() const { return fix_strwk_cg; }
        const std::vector<double> & getDepStrongWeakThrConfigs() const { return dep_strwk_thr_cg; }
        const std::vector<double> & getDepStrongWeakStrConfigs() const { return dep_strwk_str_cg; }
        const std::vector<double> & getDepStrongWeakFrgConfigs() const { return dep_strwk_frg_cg; }
        const std::vector<int> & getDepGeneralIntConfigs() const { return dep_gen_int_cg; }
        const std::vector<double> & getDepGeneralScalConfigs() const { return dep_gen_scal_cg; }
        
        // Getters for strong-weak probability for current follower optimal objective value.
        double getStrongWeakProbab(double) const;                                            // Return probability optimistic case for current decision-dependent strong-weak follower behavior and obtained optimal leader solution and optimal follower value.
        double getEvalDepStrongWeakProbab(double, Input::TypesDepStrongWeak, double) const;  // Return probability optimistic case for given decision-dependent strong-weak follower behavior and obtained optimal leader solution and optimal follower value.

        // Evaluation of general decision-dependent case for current leader decision.
        void evaluateDepGeneral(double &, double &, double &, double &, const double *, const double *, double, Input::TypesDepGeneral, int, double);
        
        // Display instance information.
        void display() const;
};

#endif
