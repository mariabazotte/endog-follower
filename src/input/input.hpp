#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>

class Input {
    public:

        enum FollowerBehavior {     /* Types of Follower Behavior. */
            FixStrongWeak = 0,      /* Strong Weak problem with Fixed Cooperation Level. */
            DepStrongWeak = 1,      /* Strong Weak problem with Decision-dependent Cooperation Level. */
            DepGeneral = 2,         /* General problem with Decision-dependent Cooperation. */
        };

        enum TypesDepStrongWeak {   /* Types of Strong Weak problem  with Decision-dependent Cooperation Level. */
            Proportional = 0,       /* Proportional. */
            Threshold = 1,          /* Threshold. */
            Strong = 2,             /* Strong or Flexible. */
            Fragile = 3             /* Fragile or Inflexible. */
        };

        enum TypesDepGeneral {      /* Types of General problem Decision-dependent Cooperation. */
            Neutral = 0,            /* Neutral. */
            GenProportional = 1,    /* Proportional. */
            StrongFragile = 2       /* Strong Fragile. */
        };

        enum SolverApproach {       /* Approaches for solving the stochastic program with decision-dependent uncertainty. */
            TR_SAA = 0,             /* Transformed SAA program. */
            DEP = 1                 /* Deterministic equivalent program. */
        };

    private:
        
        std::string instance_file;      /* File with instance dataset. */

        bool is_follower_near_optimal;  /* 0 -> Optimal follower S(x). 1 -> Near-optimal follower S(x,varepsilon). */
        double eps_near_optimal;        /* Varepsilon parameter defining near-optimal follower. */

        int int_follower_behavior;
        FollowerBehavior follower_behavior;      /* Type of Follower Behavior. */

        int int_type_dep_strongweak;
        TypesDepStrongWeak type_dep_strongweak;  /* Type of Decision-dependent strong-weak cooperation. */

        int int_type_dep_general;
        TypesDepGeneral type_dep_general;        /* Type of Decision-dependent general cooperation. */

        double fix_cooperation_level;   /* Fixed cooperation level for fixed strong weak case (beta). */

        double threshold_param;         /* Parameter > 0 defining slope/steepness if Decision-dependent strong-weak Threshold type.  */
        double strong_param;            /* Parameter > 0 defining decay rate if Decision-dependent strong-weak Strong type.  */
        double fragile_param;           /* Parameter > 0 defining growth/decay rate if Decision-dependent strong-weak Fragile type.  */
        int strwk_nb_pwl_intervals;     /* Number of intervals used at the PWL for the Decidion-dependent strong weak case. */

        int gen_nb_intervals;           /* Number of intervals to define cooperation cooeficient for Decision-Dependent General Proportional and Strong Fragile cases. */
        double gen_scaling_param;       /* Parameter > 0 defining scaling for cooperation cooeficient if Decision-dependent general Proportional and Strong Fragile cases. */
        bool gen_step_explicit;         /* 0 -> Do not use explicit formulation for step variables. 1-> Use explicit formulation for step variables. */
        
        int int_solver_approach;
        SolverApproach solver_approach;  /* Approach to solve the stochastic program with decision dependent uncertainty. */

        std::string solution_file;     /* File to write results. */
        std::ofstream solFile;         /* Ofstream to write the results in solution_file. */
        std::string comparison_file;   /* File to write comparison results. */
        std::ofstream compFile;        /* Ofstream to write the comparison results in comparison_file. */

        double seed;                   /* Seed for scenario generation. */
        double time_limit;             /* Time limit for optimization. */
        double nb_threads;             /* Number of threads for gurobi. */
        double verbose;                /* Verbose (print configuration). */
        
        double eps_bigm;               /* Small value epsilon used at bigm constraints. */

        int nbproblemsSAA;             /* Number of problems used in the SAA method. */
        int nbscenariosSAA;            /* Number of scenarios for each SAA problem. */
        int nbvalidatescenarios;       /* Number of scenarios for the validation problem. */

        void defaultParams();          /* Define default values for parameters. */
        void testParameters();         /* Test parameters value. */
        void defineSolutionFile();     /* Define name solution file. */
        
    public:
        Input() { defaultParams(); }
        Input(int argc, char* argv[]);
        ~Input(){ solFile.close(); compFile.close(); }

        // Instance
        std::string getInstanceFile() const { return instance_file; }
        double getEpsBigM() const { return eps_bigm; }

        // Problem definition
        Input::FollowerBehavior getFollowerBehavior() const { return follower_behavior; }
        Input::TypesDepStrongWeak getTypeDepStrongWeak() const { return type_dep_strongweak; }
        Input::TypesDepGeneral getTypeDepGeneral() const { return type_dep_general; }
        
        double getFixCoopLevel() const { return fix_cooperation_level; } 
        
        double getParamThreshold() const { return threshold_param; }
        double getParamStrong() const { return strong_param; }
        double getParamFragile() const { return fragile_param; }
        int getNbPWLIntervalsStrongWeak() const { return strwk_nb_pwl_intervals; }
        
        int getNbIntervalsGeneral() const { return gen_nb_intervals; }
        double getScalingGeneral() const { return gen_scaling_param; }
        bool useStepExplicitGeneral() const { return gen_step_explicit; }

        bool isFollowerNearOpt() const { return is_follower_near_optimal; }  
        double getEpsNearOpt() const { return eps_near_optimal; }     
        
        // Problem resolution
        Input::SolverApproach getSolverApproach() const { return solver_approach; }

        // General parameters
        double getTimeLimit() const { return time_limit; }
        int getNbThreads() const { return nb_threads; }
        double getSeed() const { return seed; }
        int getVerbose() const { return verbose; }

        // Scenarios 
        int getNbProblemsSAA() const { return nbproblemsSAA; }
        int getNbScenariosSAA() const { return nbscenariosSAA; }
        int getNbValidateScenarios() const { return nbvalidatescenarios; }

        // Auxiliar
        std::string doubleToString(double);

        // Solution
        void writeHead(std::string);     /* Write head of solution file. */
        void write(std::string);         /* Write information of solution file. */

        void writeHeadComp(std::string);   /* Write head of comparison file. */
        void writeComp(std::string);       /* Write information of comparison file. */

        // Display
        void display();                  /* Display information. */
};

std::ostream& operator<<(std::ostream&, const Input::FollowerBehavior &);
std::ostream& operator<<(std::ostream&, const Input::TypesDepStrongWeak &);
std::ostream& operator<<(std::ostream&, const Input::TypesDepGeneral &);
std::ostream& operator<<(std::ostream&, const Input::SolverApproach &);

#endif