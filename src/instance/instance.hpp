#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include "../input/input.hpp"
#include "bilevelmodel/bilevelmodel.hpp"

#include <list>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>

class Instance {
    private:
        const Input & input;                       /* Input. */

        BilevelModel * bilevelmodel = NULL;        /* Bilevel model with variables and constraints information. */
        
        std::string instance_file;                 /* File with instance information. */
        std::string name;                          /* Name of the instance. */ 

        double seed;                               /* Seed for random data generation (scenarios). */
        std::default_random_engine rd_engine;      /* Random engine for data generation (scenarios).*/ 

        std::uniform_real_distribution<> uniform_eps;  /* Uniform distribution for sampling scenarios. (decision-dependent strong weak) */
        std::uniform_real_distribution<> uniform;      /* Uniform distribution for sampling scenarios. (decision-dependent general) */
        std::normal_distribution<> normal;             /* Normal distribution for sampling scenarios. (decision-dependent general) */

        double critical_tstudent;                  /* Critical value for t-student distribution. */
        double critical_normal;                    /* Critical value for normal distribution. */

        // Parameters for General Decision Dependent Proportion and Strong Fragile 
        double gen_min_probab;                     /* pi_min: minimum value of the function proportional to the decision-dependent probability. */
        double gen_max_probab;                     /* pi_max: maximum value of the function proportional to the decision-dependent probability. */
        double gen_ref_pt_probab;                  /* pi_c: reference value point of the function proportional to the decision-dependent probability. */
        int gen_nb_intervals;                      /* Number of intervals to define cooperation coefficient. */
        std::vector<double> gen_intervals;         /* Interval values to define cooperation coefficient. */
        std::vector<double> gen_coeff_intervals;   /* Cooperation coefficient value for each interval. */

        int nb_saa_problems;                       /* Number of SAA problems. */
        int nb_saa_scenarios;                      /* Number of scenarios for each SAA problems. */
        int nb_val_scenarios;                      /* Number of scenarios for the validation test. */

        // Scenarios for Decision-dependent strong-weak case.
        std::vector<std::vector<double>> scenarios_dep_strwk;   /* scenarios_dep_strwk[pr][s]: inverse value of beta or cooperation level sampled from U(0,1) for transformed scenario s of the SAA problem pr. */

        // Scenarios for Decision-dependent general case.
        std::vector<std::vector<std::vector<double>>> dir_scenarios_dep_gen;   /* dir_scenarios_dep_gen[pr][s]: (upsilon) direction vector for general case transformed scenario s of the SAA problem pr. */
        std::vector<std::vector<double>> step_scenarios_dep_gen;               /* step_scenarios_dep_gen[pr][s]: (tau) sample U(0,1) to define step for general case transformed scenario s of the SAA problem pr.  */
        std::vector<std::vector<double>> accep_scenarios_dep_gen;              /* accep_scenarios_dep_gen[pr][s]: (phi) sample U(0,1) to define acceptance rate for general case transformed scenario s of the SAA problem pr. */

        std::vector<std::vector<double>> coeff_obj_alpha;                       /* coeff_obj_alpha[pr][s]: Coefficient of the alpha variable at the follower near optimal objective value constraint at problem pr and scenario s. */
        std::vector<std::vector<std::vector<double>>> coeff_alpha;              /* coeff_alpha[pr][s][c]: Coefficient of the alpha variable at the follower constraint c at problem pr and scenario s. */
        std::vector<std::vector<std::vector<std::vector<int>>>> constrs_alpha;  /* constrs_alpha[pr][s][0/1]: List of follower constraint defining alpha min (0), and alpha max (1) in problem pr and scenario s. */
        std::vector<std::vector<std::vector<std::vector<int>>>> bdlbs_alpha;    /* bdlbs_alpha[pr][s][0/1]: List of lower bounds on follower variables defining alpha min (0), and alpha max (1) in problem pr and scenario s. */
        std::vector<std::vector<std::vector<std::vector<int>>>> bdubs_alpha;    /* bdubs_alpha[pr][s][0/1]: List of upper bounds on follower variables defining alpha min (0), and alpha max (1) in problem pr and scenario s. */

        // Evaluation scenarios for Decision-dependent general case.
        std::vector<std::vector<double>> dir_eval_dep_gen;   /* dir_eval_dep_gen[s]: (upsilon) direction for general case for evaluation step.  */
        std::vector<double> step_eval_dep_gen;               /* step_eval_dep_gen[s]: (tau) sample U(0,1) to define step for general case for evaluation step. */
        std::vector<double> accep_eval_dep_gen;              /* accep_eval_dep_gen[s]: (phi) sample U(0,1) to define acceptance rate for general case for evaluation step. */
        
        std::vector<double> eval_coeff_obj_alpha;                 /* eval_coeff_obj_alpha[s]: */
        std::vector<std::vector<double>> eval_coeff_alpha;        /* eval_coeff_alpha[s][c]: */
        std::vector<std::vector<std::vector<int>>> eval_constrs_alpha;  /* eval_constrs_alpha[s][0/1]: */
        std::vector<std::vector<std::vector<int>>> eval_bdlbs_alpha;    /* eval_bdlbs_alpha[s][0/1]: */
        std::vector<std::vector<std::vector<int>>> eval_bdubs_alpha;    /* eval_bdubs_alpha[s][0/1]: */

        // Problem configurations for comparison.
        std::vector<double> fix_strwk_cg = {0.0, 0.3, 0.5, 0.7, 1.0};
        std::vector<double> dep_strwk_thr_cg = {0.5, 2.0, 5.0, 10.0};
        std::vector<double> dep_strwk_str_cg = {0.5, 2.0, 5.0, 10.0};
        std::vector<double> dep_strwk_frg_cg = {0.5, 2.0, 5.0, 10.0};
        std::vector<int> dep_gen_int_cg = {4, 5};
        std::vector<double> dep_gen_scal_cg = {0.5, 2.0, 5.0, 10.0};

        void defineName();                       /* Define instance name. */
        void defineCriticalValues();             /* Define critical values. */
        void defineDepGenParams();               /* Define auxiliar parameters for General Decision-dependent case. */
        
        void generateScenarios();                                              /* Generate scenarios. */
        double inverseDepStrongWeakScenario();                                 /* Function to sample scenarios for decision-dependent strong weak case. */
        void nullSpaceFollowerEqConstrs(Eigen::MatrixXd &);                    /* Function to define null space for equality constraints. */
        void dirDepGeneralScenario(Eigen::MatrixXd &, std::vector<double> &);  /* Function to sample direction scenarios for decision-dependent general case. */
        void defineGeneralStepInfo(int, std::vector<double> &, 
                                std::vector<std::vector<double>> &, 
                                std::vector<std::vector<std::vector<int>>> &, 
                                std::vector<std::vector<std::vector<int>>> &, 
                                std::vector<std::vector<std::vector<int>>> &);  /* Function to define scenario step information for decision-dependent general case. */

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

        double getStepSizeInterval() const { return 10; }
        double getGeneralMaxProbab() const { return gen_max_probab; }
        double getGeneralMinProbab() const { return gen_min_probab; }
        double getGeneralRefProbab() const { return gen_ref_pt_probab; }

        double getGeneralNbIntervals() const { return gen_nb_intervals; }
        double getGeneralIntervalValue(int i) const { return gen_intervals[i]; }
        double getGeneralIntervalCoeff(int i) const { return gen_coeff_intervals[i]; }
        double getRefLeaderObj() const { return (bilevelmodel->leader_lb + bilevelmodel->leader_ub)/2.0; }

        // Getters for 
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

        double getStrongWeakProbab(double) const;                                            // Return probability optimistic case for current decision-dependent strong-weak follower behavior and obtained optimal leader solution and optimal follower value.
        double getEvalDepStrongWeakProbab(double, Input::TypesDepStrongWeak, double) const;  // Return probability optimistic case for given decision-dependent strong-weak follower behavior and obtained optimal leader solution and optimal follower value.

        double getDepGeneralProbab(double, const double *, std::vector<double> &) const;
        double getEvalDepGeneralProb(double, const double *, std::vector<double> &, Input::TypesDepGeneral, int, double) const;
        
        double defineEvalMinStep(int, const double *, std::vector<double> &, double) const;
        double defineEvalMaxStep(int, const double *, std::vector<double> &, double) const;

        void display() const;
};

#endif
