#include "input.hpp"

void Input::defaultParams(){
    // Problem definition
    is_follower_near_optimal = true; 
    eps_near_optimal = 0.01;
    is_near_optimal_mult = false;
    
    // Parameters for Fixed Strong-Weak
    fix_cooperation_level = 0.5;

    // Parameters for Type of Decision-Dependent Stong Weak case
    threshold_param = 2.0;  // This value corresponds to the special case of softmax
    strong_param = 3.0;
    fragile_param = 3.0;
    strong_power_param = 2.0;
    fragile_power_param = 2.0;
    strwk_nb_pwl_intervals = 100;
    rel_error_pwl = 1;
    abs_error_pwl = 1e-4;

    // Parameters for Type of Decision-Dependent General case
    gen_nb_intervals = 4;
    gen_scaling_param = 1.0;
    defineGeneralIntervals();
    int_method_dep_gen = 1;
    method_dep_gen = MethodDepGeneral(int_method_dep_gen);
    gen_step_explicit = true;
    gen_step_explicit_grbminmax = false;
    
    // Problem solution
    int_solver_approach = 0;
    solver_approach = SolverApproach(int_solver_approach);

    // General parameters
    time_limit = 3600;
    nb_threads = 1;
    verbose = 0;
    eps_bigm = 0.001;
                  
    // Parameters for SAA or evaluation analysis
    nbproblemsSAA = 5;
    nbscenariosSAA = 100;
    nbthinningSAA = 2;

    nbvalidateproblems = 5;
    nbvalidatescenarios = 10000;
    nbvalidatethinning = 50;

    coordinate_har = true;
    metropolis_hastings = false;
}

Input::Input(int argc, char* argv[]){
    defaultParams();
    int mandatory = 0;
    bool mandatory_fix_strwk = false;
    bool mandatory_dep_strwk = false;
    bool mandatory_dep_gen = false;
    for (int i = 1; i < argc; i += 2){
        if(std::string(argv[i]) == "-instancefile"){ // Mandatory parameters
            instance_file = argv[i+1];
            mandatory += 1;
        }else if(std::string(argv[i]) == "-follower"){
            int_follower_behavior = std::stoi(argv[i+1]);
            follower_behavior = FollowerBehavior(int_follower_behavior);
            mandatory += 1;
        }else if(std::string(argv[i]) == "-fix_strwk"){
            fix_cooperation_level = std::stod(argv[i+1]);  
            mandatory_fix_strwk = true;
        }else if(std::string(argv[i]) == "-dep_strwk"){
            int_type_dep_strongweak = std::stoi(argv[i+1]);
            type_dep_strongweak = TypesDepStrongWeak(int_type_dep_strongweak);
            mandatory_dep_strwk = true;
        }else if(std::string(argv[i]) == "-dep_gen"){
            int_type_dep_general = std::stoi(argv[i+1]);
            type_dep_general = TypesDepGeneral(int_type_dep_general);
            mandatory_dep_gen = true;
        }
        else if(std::string(argv[i]) == "-near_opt")  // Optional :: Instance and Behavior
            is_follower_near_optimal = std::stoi(argv[i+1]);  
        else if(std::string(argv[i]) == "-eps_near_opt")  
            eps_near_optimal = std::stod(argv[i+1]); 
        else if(std::string(argv[i]) == "-mult_near_opt")  
            is_near_optimal_mult = std::stoi(argv[i+1]); 
        else if(std::string(argv[i]) == "-thr_param")  // Optional :: Strong Weak Decision-Dependent parameters
            threshold_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-str_param")
            strong_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-frag_param")
            fragile_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-str_pow_param")
            strong_power_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-frag_pow_param")
            fragile_power_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-nb_intv_pwl")
            strwk_nb_pwl_intervals = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nb_intv"){ // Optional :: General Decision-Dependent parameters
            gen_nb_intervals = std::stoi(argv[i+1]);
            defineGeneralIntervals();
        }else if(std::string(argv[i]) == "-scal_param")
            gen_scaling_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-method_dep_gen"){
            int_method_dep_gen = std::stoi(argv[i+1]);
            method_dep_gen = Input::MethodDepGeneral(int_method_dep_gen);
        }else if(std::string(argv[i]) == "-expl_step")
            gen_step_explicit = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-expl_step_minmax")
            gen_step_explicit_grbminmax = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-sol"){ // Optional :: Solution 
            int_solver_approach = std::stoi(argv[i+1]);
            solver_approach = SolverApproach(int_solver_approach);
        }else if(std::string(argv[i]) == "-timelimit") // Optional :: General parameters
            time_limit = std::stod(argv[i+1]);   
        else if(std::string(argv[i]) == "-nbthreads")
            nb_threads = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-verbose")
            verbose = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbproblemsSAA") // Optional :: SAA problem and validation
            nbproblemsSAA = std::stoi(argv[i+1]); 
        else if(std::string(argv[i]) == "-nbscenariosSAA")
            nbscenariosSAA = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbthinningSAA")
            nbthinningSAA = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbvalidateproblems")
            nbvalidateproblems = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbvalidatescenarios")
            nbvalidatescenarios = std::stoi(argv[i+1]); 
        else if(std::string(argv[i]) == "-nbvalidatethinning")
            nbvalidatethinning = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-coord_har")
            coordinate_har = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-metrop_has")
            metropolis_hastings = std::stoi(argv[i+1]);
        else{
            std::cerr << "ERROR: Argument '" << argv[i] << "' not defined." << std::endl;
            throw std::runtime_error(std::string("Incorrect line of command"));
        }
    }
    testParameters();
    if(mandatory != 2){
        std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
        std::cerr << "You need to define: -instancefile -follower and -fix_coop_strweak (if -follower = 0), -dep_coop_strweak (if -follower = 1), or -dep_coop_gen (if -follower = 2)." << std::endl;
        throw std::runtime_error(std::string("Incorrect line of command"));
    }
    if(follower_behavior == Input::FollowerBehavior::FixStrongWeak && mandatory_fix_strwk == false){
        std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
        std::cerr << "You need to define: -fix_coop_strweak when -follower = 0." << std::endl;
        throw std::runtime_error(std::string("Incorrect line of command"));
    }
    if(follower_behavior == Input::FollowerBehavior::DepStrongWeak && mandatory_dep_strwk == false){
        std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
        std::cerr << "You need to define: -dep_coop_strweak when -follower = 1." << std::endl;
        throw std::runtime_error(std::string("Incorrect line of command"));     
    }
    if(follower_behavior == Input::FollowerBehavior::DepGeneral && mandatory_dep_gen == false){
        std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
        std::cerr << "You need to define: -dep_coop_gen when -follower = 2." << std::endl;
        throw std::runtime_error(std::string("Incorrect line of command"));  
    }
    defineSolutionFile();
}

// Auxiliar input for computing bounds.
Input::Input(const Input & input, double fix_strwk){
    // Problem definition.
    is_follower_near_optimal = input.isFollowerNearOpt(); 
    eps_near_optimal = input.getEpsNearOpt();
    is_near_optimal_mult = input.isNearOptMult();

    // Define it as fix strong weak.
    int_follower_behavior = 0;
    follower_behavior = FollowerBehavior(int_follower_behavior);

    // Parameters for Fixed Strong-Weak
    fix_cooperation_level = fix_strwk;

    // Define instance.
    instance_file = input.getInstanceFile();

    // Problem solution
    int_solver_approach = 1;
    solver_approach = SolverApproach(int_solver_approach);

    // General parameters
    time_limit = 600;
    nb_threads = input.getNbThreads();
    verbose = 1;
    eps_bigm = 0.001;
}

void Input::defineGeneralIntervals(){
    if(gen_nb_intervals > 0){
        gen_intervals.resize(gen_nb_intervals+1);
        gen_coeff_intervals.resize(gen_nb_intervals);
        double interval = 1.0 / static_cast<double>(gen_nb_intervals);
        for(int i = 0; i <= gen_nb_intervals; ++i) 
            gen_intervals[i] = i*interval;
        gen_max_coeff_intervals = 0.0;
        for(int i = 0; i < gen_nb_intervals; ++i){
            // gen_coeff_intervals[i] = 0.5 - (gen_intervals[i]+gen_intervals[i+1])/2.0;
            gen_coeff_intervals[i] = 1.0 - (gen_intervals[i]+gen_intervals[i+1]);
            gen_max_coeff_intervals = std::max(gen_max_coeff_intervals,gen_coeff_intervals[i]);
        }
        gen_max_coeff_intervals += 1e-8;
    }
}

void Input::defineSolutionFile(){
    std::string instfile = instance_file;

    while(instfile.find("/") != std::string::npos)
        instfile = instfile.substr(instfile.find("/")+1);
    
    solution_file = "../results/";
    
    // Configuration corresponding to follower behavior.
    if(follower_behavior == Input::FollowerBehavior::FixStrongWeak){
        solution_file += "fixstrongweak/param" + doubleToString(fix_cooperation_level) + "_";
    }
    if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
        solution_file += "depstrongweak/";
        if(solver_approach == Input::SolverApproach::TR_SAA)
            solution_file += "saa/";
        if(solver_approach == Input::SolverApproach::DEP)
            solution_file += "dep/";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Proportional)
            solution_file += "proportional/";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Threshold)
            solution_file += "threshold/param" + doubleToString(threshold_param) + "_";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Strong)
            solution_file += "strong/param" + doubleToString(strong_param) + "_";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Fragile)
            solution_file += "fragile/param" + doubleToString(fragile_param) + "_";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::StrongPower)
            solution_file += "strong_power/param" + doubleToString(strong_power_param) + "_";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::FragilePower)
            solution_file += "fragile_power/param" + doubleToString(fragile_power_param) + "_";
    }
    if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        solution_file += "depgeneral/";
        if(type_dep_general == Input::TypesDepGeneral::Neutral)
            solution_file += "neutral/";
        if(type_dep_general == Input::TypesDepGeneral::GenProportional)
            solution_file += "proportional/nbintv" + std::to_string(gen_nb_intervals) + "_scal" + doubleToString(gen_scaling_param) + "_"; 
        if(type_dep_general == Input::TypesDepGeneral::GenFragile)
            solution_file += "fragile/nbintv" + std::to_string(gen_nb_intervals) + "_scal" + doubleToString(gen_scaling_param) + "_";
        if(type_dep_general == Input::TypesDepGeneral::GenFragilePower)
            solution_file += "fragile_power/nbintv" + std::to_string(gen_nb_intervals) + "_scal" + doubleToString(gen_scaling_param) + "_";
        if(type_dep_general == Input::TypesDepGeneral::GenStrongPower)
            solution_file += "strong_power/nbintv" + std::to_string(gen_nb_intervals) + "_scal" + doubleToString(gen_scaling_param) + "_";

        if(method_dep_gen == Input::MethodDepGeneral::SingleLevel)
            solution_file += "sglevel_explstep" + std::to_string(gen_step_explicit) + "_grbminmax" + std::to_string(gen_step_explicit_grbminmax) + "_";
        if(method_dep_gen == Input::MethodDepGeneral::Relax)
            solution_file+= "relax_";
        if(method_dep_gen == Input::MethodDepGeneral::Penalty)
            solution_file+= "penalty_";
        if(method_dep_gen == Input::MethodDepGeneral::LocalSearch)
            solution_file+= "ls_";
        if(method_dep_gen == Input::MethodDepGeneral::Callback)
            solution_file+= "callback_";
    }

    // Configuration corresponding to follower optimiality.
    solution_file += "nearopt" + std::to_string(is_follower_near_optimal);
    if(is_follower_near_optimal == true) {
        solution_file += "_epsnearopt" + doubleToString(eps_near_optimal);
        solution_file += "_multnearopt" + std::to_string(is_near_optimal_mult);
    }

    // Configuration corresponding to instance.
    solution_file += "_" + instfile.substr(0, instfile.size()-4);

    // Configuration corresponding to SAA method: SAA problem and evaluation problem.
    if(solver_approach != Input::SolverApproach::DEP) {   
        solution_file += "_nbsaapr" + std::to_string(nbproblemsSAA) + 
                         "_nbsaasc" + std::to_string(nbscenariosSAA) + 
                         "_nbsaath" + std::to_string(nbthinningSAA);

        if(follower_behavior == Input::FollowerBehavior::DepGeneral) {
            // Decision dependent general case is evaluated with markov
            // chain monte carlo hit and run. Strong-weak cases are evaluated
            // by enumerating the two endogenous scenarios, no need for MCMC.
            solution_file += "_nbvalpr" + std::to_string(nbvalidateproblems) + 
                             "_nbvalsc" + std::to_string(nbvalidatescenarios) + 
                             "_nbvalth" + std::to_string(nbvalidatethinning);
            solution_file += "_coordhar" + std::to_string(coordinate_har);
            solution_file += "_metrophas" + std::to_string(metropolis_hastings);
        }
    }else{
        if(follower_behavior == Input::FollowerBehavior::DepStrongWeak)
            solution_file += "_errpwl" + std::to_string(abs_error_pwl);
            // solution_file += "_errpwl" + std::to_string(rel_error_pwl);
            // solution_file += "_intpwl" + std::to_string(strwk_nb_pwl_intervals);
    }

    // Initializing comparison file with information on model solved to obtain leader decision. 
    comparison_file = solution_file;

    if(follower_behavior != Input::FollowerBehavior::DepGeneral){
        // Configuration to evaluate other configurations with the decision obtained
        // by the considered problem.
        comparison_file += "_nbvalpr" + std::to_string(nbvalidateproblems) + 
                            "_nbvalsc" + std::to_string(nbvalidatescenarios) + 
                            "_nbvalth" + std::to_string(nbvalidatethinning);
        comparison_file += "_coordhar" + std::to_string(coordinate_har);
        comparison_file += "_metrophas" + std::to_string(metropolis_hastings);
    }

    // Configuration corresponding to general parameters.
    int aux_time = time_limit;
    solution_file +=  "_time" + std::to_string(aux_time);
    comparison_file +=  "_time" + std::to_string(aux_time);

    comparison_file += "_comparison.csv";
    solution_file += "_solution.csv"; 
    
    // Open files
    solFile.open(solution_file, std::ios_base::out);
    compFile.open(comparison_file, std::ios_base::out);
}

std::string Input::doubleToString(double param){
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(2);
    oss << param;
    return oss.str();
}

void Input::testParameters(){
    if(int_follower_behavior > 2){
        throw std::runtime_error(std::string("Wrong follower behavior value."));
        exit(0);
    }
    if(follower_behavior == Input::FollowerBehavior::FixStrongWeak){
        if(fix_cooperation_level <= (-1e-8) || fix_cooperation_level >= (1.0 + 1e-8)){
            throw std::runtime_error(std::string("Wrong value for fixed cooperation parameter: it should be defined in the interval [0.0,1.0]"));
            exit(0);
        }
    }
    if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
        if(int_type_dep_strongweak > 5){
            throw std::runtime_error(std::string("Wrong type of strong weak decision-dependent cooperation value."));
            exit(0);
        }
    }
    if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        if(int_type_dep_general > 4){
            throw std::runtime_error(std::string("Wrong type of general decision-dependent cooperation value."));
            exit(0);
        }
        if(int_method_dep_gen > 4){
            throw std::runtime_error(std::string("Wrong type of general decision-dependent method value."));
            exit(0);
        }
        if(int_method_dep_gen == 2){
            throw std::runtime_error(std::string("Penalty method for general decision-dependent method is obsolete. Change method."));
            exit(0);
        }
    }

    if(follower_behavior == Input::FollowerBehavior::FixStrongWeak){ // Use deterministic equivalent
        int_solver_approach = 1;
        solver_approach = SolverApproach(int_solver_approach);
        std::cout << " Fixed Strong Weak type always solved with DEP approach." << std::endl;
    }else if(follower_behavior == Input::FollowerBehavior::DepGeneral){ // Use transformed SAA
        int_solver_approach = 0;
        solver_approach = SolverApproach(int_solver_approach);
        std::cout << " Decision-dependent general type always solved with TR-SAA approach." << std::endl;
    }else if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
        if(solver_approach == Input::SolverApproach::TR_SAA){
            // Thinning always 1 in this case. We use monte carlo sampling in this case, not 
            // markov chain monte carlo.
            nbthinningSAA = 1; 
        }
    }

    if(int_solver_approach > 1){
        throw std::runtime_error(std::string("Wrong solver approach value."));
        exit(0);
    }

    // If SAA is not used, then there is no SAA problem, only evaluation.
    if(solver_approach == Input::SolverApproach::DEP){ 
        nbproblemsSAA = 0;
        nbscenariosSAA = 0;
        nbthinningSAA = 0;
    }

    if(is_follower_near_optimal == true){
        if(eps_near_optimal <= (-1e-8) || eps_near_optimal >= (1.0 + 1e-8)){
            throw std::runtime_error(std::string("Wrong value for epsilon near-optimality: it should be defined in the interval [0.0,1.0]"));
            exit(0);
        }
    }

    // Verify parameters.
    if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Threshold){
            if(threshold_param <= -1e-8){
                std::cout << "Threshold parameter must be > 0. Setting parameter to default value 2.0." << std::endl;
                threshold_param = 2.0;
            }
        }
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Strong){
            if(strong_param <= -1e-8){
                std::cout << "Strong parameter must be > 0. Setting parameter to default value 3.0." << std::endl;
                strong_param = 3.0;
            }
        }
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Fragile){
            if(fragile_param <= -1e-8){
                std::cout << "Fragile parameter must be > 0. Setting parameter to default value 3.0." << std::endl;
                fragile_param = 3.0;
            }
        }
        if(type_dep_strongweak == Input::TypesDepStrongWeak::StrongPower){
            if(strong_power_param <= (1.0 - 1e-8)){
                std::cout << "Strong power parameter must be > 1.0. Setting parameter to default value 2.0." << std::endl;
                strong_power_param = 2.0;
            }
        }
        if(type_dep_strongweak == Input::TypesDepStrongWeak::FragilePower){
            if(fragile_param <= (1.0 - 1e-8)){
                std::cout << "Fragile power parameter must be > 1.0. Setting parameter to default value 2.0." << std::endl;
                fragile_power_param = 2.0;
            }
        }
    }

    if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        if(type_dep_general != Input::TypesDepGeneral::Neutral){
            if(gen_nb_intervals <= (-1 - 1e-8)){
                std::cout << "Nb. intervals parameter must be >= -1. Setting parameter to default value 4." << std::endl;
                gen_nb_intervals = 4;
            }

            if((method_dep_gen != Input::MethodDepGeneral::Callback && method_dep_gen != Input::MethodDepGeneral::LocalSearch) && gen_nb_intervals <= 1e-8){
                std::cout << "Nb. intervals parameter must be > 0 for methods different from callback and local search. Setting parameter to default value 4." << std::endl;
                gen_nb_intervals = 4;
            }

            if(type_dep_general == Input::TypesDepGeneral::GenProportional){
                if(gen_scaling_param <= (1.0 - 1e-8) || gen_scaling_param >= (1.0 + 1e-8) ){
                    std::cout << "Scaling parameter must be equal to 1.0 for general proportional case. Setting parameter to value 1.0." << std::endl;
                    gen_scaling_param = 1.0;
                }
            }

            if(type_dep_general == Input::TypesDepGeneral::GenFragile || type_dep_general == Input::TypesDepGeneral::GenFragilePower){
                if(gen_scaling_param <= -1e-8){
                    std::cout << "Scaling parameter must be > 0 for general fragile or general fragile power case. Setting parameter to default value 1.0." << std::endl;
                    gen_scaling_param = 1.0;
                }
            }

            if(type_dep_general == Input::TypesDepGeneral::GenStrongPower){
                if(gen_scaling_param <= (1.0 - 1e-8)){
                    std::cout << "Scaling parameter must be >= 1.0 for general strong power case. Setting parameter to default value 1.0." << std::endl;
                    gen_scaling_param = 1.0;
                }
            }
        }

        if(type_dep_general != Input::TypesDepGeneral::Neutral){
            if(method_dep_gen == Input::MethodDepGeneral::Penalty || method_dep_gen == Input::MethodDepGeneral::Relax || method_dep_gen == Input::MethodDepGeneral::LocalSearch || method_dep_gen == Input::MethodDepGeneral::Callback){
                if(metropolis_hastings == true){
                    std::cout << "Relax, Penalty, Local Search and Callback methods works with slice sampling formulation. Setting metropolis hastings to false." << std::endl;
                    metropolis_hastings = false;  
                }
            }
        }

        if(method_dep_gen == Input::MethodDepGeneral::SingleLevel){
            if(gen_step_explicit == false){
                gen_step_explicit_grbminmax = false;
            }
        }
    }
}

void Input::display(){
    std::cout << "----------- INPUT PARAMS -------------" << std::endl;
	std::cout << "INSTANCE  FILE   : '" << instance_file << "'" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "SOLUTION  FILE   : '" << solution_file << "'" << std::endl;
    std::cout << "COMPARISON  FILE : '" << comparison_file << "'" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "FOLLOWER BEHAVIOR:  " << follower_behavior << std::endl;
    if(follower_behavior == Input::FollowerBehavior::FixStrongWeak)
        std::cout << "COOPERATION TYPE :  " << fix_cooperation_level << std::endl;
    if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
        std::cout << "COOPERATION TYPE :  " << type_dep_strongweak << std::endl;
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Threshold)
            std::cout << "THRESHOLD PARAM  :  " << threshold_param << std::endl;
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Strong)
            std::cout << "STRONG PARAM     :  " << strong_param << std::endl;
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Fragile)
            std::cout << "FRAGILE PARAM    :  " << fragile_param << std::endl;
        if(type_dep_strongweak == Input::TypesDepStrongWeak::StrongPower)
            std::cout << "STRONG POW PARAM :  " << strong_power_param << std::endl;
        if(type_dep_strongweak == Input::TypesDepStrongWeak::FragilePower)
            std::cout << "FRAGILE POW PARAM:  " << fragile_power_param << std::endl;
    }if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        std::cout << "COOPERATION TYPE :  " << type_dep_general << std::endl;
        if(type_dep_general != Input::TypesDepGeneral::Neutral){
            std::cout << "NB INTERVALS     :  " << gen_nb_intervals << std::endl;
            std::cout << "SCALING PARAM    :  " << gen_scaling_param << std::endl;
        }
        std::cout << "METHOD           :  " << method_dep_gen << std::endl;
        if(method_dep_gen == Input::MethodDepGeneral::SingleLevel){
            std::cout << "EXPL STEP        :  " << gen_step_explicit << std::endl;
            std::cout << "GRB MAX/MIN      :  " << gen_step_explicit_grbminmax << std::endl;
        }
    }std::cout << "--------------------------------------" << std::endl;
    std::cout << "NEAR OPTIMAL     :  " << is_follower_near_optimal << std::endl;
    if(is_follower_near_optimal == true) {
        std::cout << "EPS NEAR OPTIMAL :  " << eps_near_optimal << std::endl;
        std::cout << "MULT NEAR OPTIMAL:  " << is_near_optimal_mult << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "SOLUTION APPROACH:  " << solver_approach << std::endl;  
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "TIME LIMIT       :  " << time_limit << std::endl;
    std::cout << "NUMBER THREADS   :  " << nb_threads << std::endl;
    std::cout << "VERBOSE          :  " << verbose << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    if(solver_approach == Input::SolverApproach::TR_SAA) {
        std::cout << "NB SAA PROBLEMS  :  " << nbproblemsSAA << std::endl;
        std::cout << "NB SAA SCENARIOS :  " << nbscenariosSAA << std::endl;
        std::cout << "NB SAA THINNING  :  " << nbthinningSAA << std::endl;
        std::cout << "--------------------------------------" << std::endl;
    }
    std::cout << "NB VAL PROBLEMS  :  " << nbvalidateproblems << std::endl;
    std::cout << "NB VAL SCENARIOS :  " << nbvalidatescenarios << std::endl;
    std::cout << "NB VAL THINNING  :  " << nbvalidatethinning << std::endl;
    std::cout << "COORDINATE HAR.  :  " << coordinate_har << std::endl;
    std::cout << "METROPOLIS HAS.  :  " << metropolis_hastings << std::endl;
    std::cout << "--------------------------------------" << std::endl;
}

void Input::write(std::string output){
    if(!solFile.fail()){
        solFile << follower_behavior << ";";
        if(follower_behavior == Input::FollowerBehavior::FixStrongWeak)
            solFile << fix_cooperation_level << ";";
        if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
            solFile << type_dep_strongweak << ";";
            if(type_dep_strongweak == Input::TypesDepStrongWeak::Threshold) solFile << threshold_param << ";";
            if(type_dep_strongweak == Input::TypesDepStrongWeak::Strong) solFile << strong_param << ";";
            if(type_dep_strongweak == Input::TypesDepStrongWeak::Fragile) solFile << fragile_param << ";";
            if(type_dep_strongweak == Input::TypesDepStrongWeak::StrongPower) solFile << strong_power_param << ";";
            if(type_dep_strongweak == Input::TypesDepStrongWeak::FragilePower) solFile << fragile_power_param << ";";
        }if(follower_behavior == Input::FollowerBehavior::DepGeneral){
            solFile << type_dep_general << ";";
            if(type_dep_general != Input::TypesDepGeneral::Neutral){
                solFile << gen_nb_intervals << ";";
                solFile << gen_scaling_param << ";";
            }
            solFile << method_dep_gen << ";";
            if(method_dep_gen == Input::MethodDepGeneral::SingleLevel){
                solFile << gen_step_explicit << ";";
                solFile << gen_step_explicit_grbminmax << ";";

            }
        }
        solFile << solver_approach << ";";
        solFile << is_follower_near_optimal << ";";
        solFile << eps_near_optimal << ";";
        solFile << is_near_optimal_mult << ";";
        solFile << instance_file << ";";
        if(solver_approach == Input::SolverApproach::TR_SAA){
            solFile << nbproblemsSAA << ";";
            solFile << nbscenariosSAA << ";";
            solFile << nbthinningSAA << ";";
            if(follower_behavior == Input::FollowerBehavior::DepGeneral){
                solFile << nbvalidateproblems << ";";
                solFile << nbvalidatescenarios << ";";
                solFile << nbvalidatethinning << ";";
                solFile << coordinate_har << ";";
                solFile << metropolis_hastings << ";";
            }
        }
        solFile << instance_align << ";";
        solFile << output;
        solFile << std::endl;
    }
}

void Input::writeHead(std::string output){
    if(!solFile.fail()){
        solFile << "FOLLOWER_BEHAVIOR;COOPERATION_TYPE;";
        if(follower_behavior == Input::FollowerBehavior::DepStrongWeak &&
           type_dep_strongweak != Input::TypesDepStrongWeak::Proportional) 
            solFile << "PARAM;";
        if(follower_behavior == Input::FollowerBehavior::DepGeneral){
            if(type_dep_general != Input::TypesDepGeneral::Neutral)
                solFile << "NB_INTERVALS;SCALING;";
            solFile << "METHOD;";
            if(method_dep_gen == Input::MethodDepGeneral::SingleLevel)
                solFile << "EXPL_STEP;GRB_MINMAX;";
        }
        solFile << "SOLUTION;NEAR_OPT;EPS_NEAR_OPT;MULT_NEAR_OPT;INSTANCE_FILE;";
        if(solver_approach == Input::SolverApproach::TR_SAA){
            solFile << "NB_SAA_PROBLEMS;NB_SAA_SCENARIOS;NB_SAA_THINNING;";
            if(follower_behavior == Input::FollowerBehavior::DepGeneral){
                solFile << "NB_VAL_PROBLEMS;NB_VAL_SCENARIOS;NB_VAL_THINNING;COORDINATE_HAR;METROPOLIS_HAS;";
            }
        }
        solFile << "ALIGNMENT;";
        solFile << output;
        solFile << std::endl;
    }
}

void Input::writeComp(std::string output){
    if(!compFile.fail()){
        compFile << output;
        compFile << std::endl;
    }
}

void Input::writeHeadComp(std::string output){
    if(!compFile.fail()){
        compFile << output;
        compFile << std::endl;
    }
}

std::ostream& operator<<(std::ostream& lhs, const Input::FollowerBehavior & follower_behavior) {
    switch(follower_behavior) {
        case Input::FollowerBehavior::FixStrongWeak: {
            lhs << "FixedStrongWeak";
            break;
        }
        case Input::FollowerBehavior::DepStrongWeak: {
            lhs << "DependentStrongWeak";
            break;
        }
        case Input::FollowerBehavior::DepGeneral: {
            lhs << "DependentGeneral";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::TypesDepStrongWeak & type_dep_strongweak) {
    switch(type_dep_strongweak) {
        case Input::TypesDepStrongWeak::Proportional: {
            lhs << "Proportional";
            break;
        }
        case Input::TypesDepStrongWeak::Threshold: {
            lhs << "Threshold";
            break;
        }
        case Input::TypesDepStrongWeak::Strong: {
            lhs << "Strong";
            break;
        }
        case Input::TypesDepStrongWeak::Fragile: {
            lhs << "Fragile";
            break;
        }
        case Input::TypesDepStrongWeak::StrongPower: {
            lhs << "StrongPower";
            break;
        }
        case Input::TypesDepStrongWeak::FragilePower: {
            lhs << "FragilePower";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::TypesDepGeneral & type_dep_general) {
    switch(type_dep_general) {
        case Input::TypesDepGeneral::Neutral: {
            lhs << "Neutral";
            break;
        }
        case Input::TypesDepGeneral::GenProportional: {
            lhs << "GenProportional";
            break;
        }
        case Input::TypesDepGeneral::GenFragile: {
            lhs << "GenFragile";
            break;
        }
        case Input::TypesDepGeneral::GenFragilePower: {
            lhs << "GenFragilePower";
            break;
        }
        case Input::TypesDepGeneral::GenStrongPower: {
            lhs << "GenStrongPower";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::SolverApproach & solver_approach) {
    switch(solver_approach) {
        case Input::SolverApproach::TR_SAA: {
            lhs << "Transformed_SAA";
            break;
        }
        case Input::SolverApproach::DEP: {
            lhs << "Deterministic_Equivalent";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::MethodDepGeneral & method) {
    switch(method) {
        case Input::MethodDepGeneral::SingleLevel: {
            lhs << "SingleLevel";
            break;
        }
        case Input::MethodDepGeneral::Relax: {
            lhs << "Relax";
            break;
        }
        case Input::MethodDepGeneral::Penalty: {
            lhs << "Penalty";
            break;
        }
        case Input::MethodDepGeneral::LocalSearch: {
            lhs << "LocalSearch";
            break;
        }
        case Input::MethodDepGeneral::Callback: {
            lhs << "Callback";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}