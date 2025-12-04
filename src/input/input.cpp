#include "input.hpp"

void Input::defaultParams(){
    // Problem definition
    is_follower_near_optimal = false; // Follower is optimal
    eps_near_optimal = 0.05;
    is_near_optimal_mult = true;
    fix_cooperation_level = 0.5;

    // Parameters for Type of Decision-Dependent Stong Weak case
    threshold_param = 2.0;  // This value corresponds to the special case of softmax
    strong_param = 3.0;
    fragile_param = 3.0;
    strwk_nb_pwl_intervals = 200;

    // Parameters for Type of Decision-Dependent General case
    gen_nb_intervals = 4;
    gen_scaling_param = 2.0;
    gen_step_explicit = true;
    gen_intervals.resize(gen_nb_intervals+1);
    gen_coeff_intervals.resize(gen_nb_intervals);
    double interval = 1.0 / static_cast<double>(gen_nb_intervals);
    for(int i = 0; i <= gen_nb_intervals; ++i){ 
        gen_intervals[i] = i*interval;
        if(i > 0) gen_coeff_intervals[i-1] = 0.5 - (gen_intervals[i-1]+gen_intervals[i])/2.0;
    }

    // Problem solution
    int_solver_approach = 0;
    solver_approach = SolverApproach(int_solver_approach);

    // General parameters
    time_limit = 3600;
    nb_threads = 1;
    seed = 0;
    verbose = 0;
    eps_bigm = 0.001;
                  
    // Parameters for SAA or evaluation analysis
    nbproblemsSAA = 5;
    nbscenariosSAA = 100;
    nbthinningSAA = 5;

    nbvalidateproblems = 1;
    nbvalidatescenarios = 10000;
    nbvalidatethinning = 50;
}

Input::Input(int argc, char* argv[]){
    defaultParams();
    int mandatory = 0;
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
            if(follower_behavior == FollowerBehavior::FixStrongWeak) mandatory += 1;
        }else if(std::string(argv[i]) == "-dep_strwk"){
            int_type_dep_strongweak = std::stoi(argv[i+1]);
            type_dep_strongweak = TypesDepStrongWeak(int_type_dep_strongweak);
            if(follower_behavior == FollowerBehavior::DepStrongWeak) mandatory += 1;
        }else if(std::string(argv[i]) == "-dep_gen"){
            int_type_dep_general = std::stoi(argv[i+1]);
            type_dep_general = TypesDepGeneral(int_type_dep_general);
            if(follower_behavior == FollowerBehavior::DepGeneral) mandatory += 1;
        }else if(std::string(argv[i]) == "-near_opt")  // Optional :: Instance and Behavior
            is_follower_near_optimal = std::stoi(argv[i+1]);  
        else if(std::string(argv[i]) == "-eps_near_opt")  
            eps_near_optimal = std::stod(argv[i+1]); 
        else if(std::string(argv[i]) == "-thr_param")  // Optional :: Strong Weak Decision-Dependent parameters
            threshold_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-str_param")
            strong_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-frag_param")
            fragile_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-nb_intv_pwl")
            strwk_nb_pwl_intervals = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nb_intv"){ // Optional :: General Decision-Dependent parameters
            gen_nb_intervals = std::stoi(argv[i+1]);
            gen_intervals.clear();
            gen_coeff_intervals.clear();
            gen_intervals.resize(gen_nb_intervals+1);
            gen_coeff_intervals.resize(gen_nb_intervals);
            double interval = 1.0 / static_cast<double>(gen_nb_intervals);
            for(int i = 0; i <= gen_nb_intervals; ++i){ 
                gen_intervals[i] = i*interval;
                if(i > 0) gen_coeff_intervals[i-1] = 0.5 - (gen_intervals[i-1]+gen_intervals[i])/2.0;
            }
        }else if(std::string(argv[i]) == "-scal_param")
            gen_scaling_param = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-expl_step")
            gen_step_explicit = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-sol"){ // Optional :: Solution 
            int_solver_approach = std::stoi(argv[i+1]);
            solver_approach = SolverApproach(int_solver_approach);
        }else if(std::string(argv[i]) == "-timelimit") // Optional :: General parameters
            time_limit = std::stod(argv[i+1]);   
        else if(std::string(argv[i]) == "-nbthreads") 
            nb_threads = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-seed") 
            seed = std::stod(argv[i+1]);
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
        else{
            std::cerr << "ERROR: Argument '" << argv[i] << "' not defined." << std::endl;
            throw std::runtime_error(std::string("Incorrect line of command"));
        }
    }
    testParameters();
    if(mandatory < 3){
        std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
        std::cerr << "You need to define: -instancefile -follower and -fix_coop_strweak (if -follower = 0), -dep_coop_strweak (if -follower = 1), or -dep_coop_gen (if -follower = 2)." << std::endl;
        throw std::runtime_error(std::string("Incorrect line of command"));
    }
    defineSolutionFile();
}

void Input::defineSolutionFile(){
    std::string instfile = instance_file;

    while(instfile.find("/") != std::string::npos)
        instfile = instfile.substr(instfile.find("/")+1);
    
    solution_file = "../results/";
    
    if(follower_behavior == Input::FollowerBehavior::FixStrongWeak){
        solution_file += std::string("fixstrongweak/") + std::string("param") + doubleToString(fix_cooperation_level) + "_";
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
            solution_file += std::string("threshold/") + std::string("param") + doubleToString(threshold_param) + "_";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Strong)
            solution_file += std::string("strong/") + std::string("param") + doubleToString(strong_param) + "_";
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Fragile)
            solution_file += std::string("fragile/") + std::string("param") + doubleToString(fragile_param) + "_";
    }
    if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        solution_file += "depgeneral/";
        if(type_dep_general == Input::TypesDepGeneral::Neutral)
            solution_file += std::string("neutral/") + std::string("explstep") + std::to_string(gen_step_explicit) + "_";
        if(type_dep_general == Input::TypesDepGeneral::GenProportional)
            solution_file += std::string("proportional/") + std::string("explstep") + std::to_string(gen_step_explicit) + std::string("_nbintv") + std::to_string(gen_nb_intervals) + "_scal" + doubleToString(gen_scaling_param) + "_";
        if(type_dep_general == Input::TypesDepGeneral::StrongFragile)
            solution_file += std::string("strongfragile/") + std::string("explstep") + std::to_string(gen_step_explicit) + "_nbintv" + std::to_string(gen_nb_intervals) + "_scal" + doubleToString(gen_scaling_param) + "_";
    }

    solution_file += "nearopt" + std::to_string(is_follower_near_optimal);
    if(is_follower_near_optimal == true) solution_file += "_epsnearopt" + doubleToString(eps_near_optimal);

    // Configuration corresponding to instance
    int seed_aux = seed;
    solution_file += "_" + instfile.substr(0, instfile.size()-4) + "_seed" + std::to_string(seed_aux);
    
    // Configuration corresponding to SAA 
    if(solver_approach == Input::SolverApproach::DEP) solution_file += "_nbval" + std::to_string(nbvalidatescenarios);
    else solution_file += "_nbpr" + std::to_string(nbproblemsSAA) + "_nbsc" + std::to_string(nbscenariosSAA) + "_nbval" + std::to_string(nbvalidatescenarios); 
    
    // Configuration corresponding to general parameters
    int aux_time = time_limit;
    solution_file +=  "_time" + std::to_string(aux_time);
    
    comparison_file = solution_file + "_comparison.csv";
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
    if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
        if(int_type_dep_strongweak > 3){
            throw std::runtime_error(std::string("Wrong type of strong weak decision-dependent cooperation value."));
            exit(0);
        }
    }
    if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        if(int_type_dep_general > 2){
            throw std::runtime_error(std::string("Wrong type of general decision-dependent cooperation value."));
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
    }
    if(int_solver_approach > 1){
        throw std::runtime_error(std::string("Wrong solver approach value."));
        exit(0);
    }

    // If SAA is not used, then there is no SAA problem, only evaluation.
    if(solver_approach == Input::SolverApproach::DEP){ 
        nbproblemsSAA = 0;
        nbscenariosSAA = 0;
    }

    // Verify parameters.
    if(follower_behavior == Input::FollowerBehavior::DepStrongWeak){
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Threshold){
            if(threshold_param <= 0.0){
                std::cout << "Threshold parameter must be > 0. Setting parameter to default value 2.0." << std::endl;
                threshold_param = 2.0;
            }
        }
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Strong){
            if(strong_param <= 0.0){
                std::cout << "Strong parameter must be > 0. Setting parameter to default value 3.0." << std::endl;
                strong_param = 3.0;
            }
        }
        if(type_dep_strongweak == Input::TypesDepStrongWeak::Fragile){
            if(fragile_param <= 0.0){
                std::cout << "Fragile parameter must be > 0. Setting parameter to default value 3.0." << std::endl;
                fragile_param = 3.0;
            }
        }
    }
    if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        if(type_dep_general != Input::TypesDepGeneral::Neutral){
            if(gen_scaling_param <= 0.0){
                std::cout << "Scaling parameter must be > 0. Setting parameter to default value 2.0." << std::endl;
                gen_scaling_param = 2.0;
            }
            if(gen_nb_intervals <= 0){
                std::cout << "Nb. intervals parameter must be > 0. Setting parameter to default value 4." << std::endl;
                gen_nb_intervals = 4;

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
    }if(follower_behavior == Input::FollowerBehavior::DepGeneral){
        std::cout << "COOPERATION TYPE :  " << type_dep_general << std::endl;
        std::cout << "EXPL STEP        :  " << gen_step_explicit << std::endl;
        if(type_dep_general == Input::TypesDepGeneral::GenProportional || 
            type_dep_general == Input::TypesDepGeneral::StrongFragile){
            std::cout << "NB INTERVALS     :  " << gen_nb_intervals << std::endl;
            std::cout << "SCALING PARAM    :  " << gen_scaling_param << std::endl;
        }
    }std::cout << "--------------------------------------" << std::endl;
    std::cout << "NEAR OPTIMAL     :  " << is_follower_near_optimal << std::endl;
    if(is_follower_near_optimal == true)
        std::cout << "EPS NEAR OPTIMAL :  " << eps_near_optimal << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "SOLUTION APPROACH:  " << solver_approach << std::endl;  
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "SEED             :  " << seed << std::endl;
    std::cout << "TIME LIMIT       :  " << time_limit << std::endl;
    std::cout << "NUMBER THREADS   :  " << nb_threads << std::endl;
    std::cout << "VERBOSE          :  " << verbose << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "NB PROBLEMS      :  " << nbproblemsSAA << std::endl;
    std::cout << "NB SCENARIOS     :  " << nbscenariosSAA << std::endl;
    std::cout << "NB VAL SCENARIOS :  " << nbvalidatescenarios << std::endl;
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
        }if(follower_behavior == Input::FollowerBehavior::DepGeneral){
            solFile << type_dep_general << ";";
            solFile << gen_step_explicit << ";";
            if(type_dep_general == Input::TypesDepGeneral::GenProportional || 
                type_dep_general == Input::TypesDepGeneral::StrongFragile){
                solFile << gen_nb_intervals << ";";
                solFile << gen_scaling_param << ";";
            }
        }
        solFile << solver_approach << ";";
        solFile << is_follower_near_optimal << ";";
        solFile << eps_near_optimal << ";";
        solFile << instance_file << ";";
        solFile << seed << ";";
        solFile << nbproblemsSAA << ";";
        solFile << nbscenariosSAA << ";";
        solFile << nbvalidatescenarios << ";";
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
            solFile << "EXPL_STEP;";
            if(type_dep_general != Input::TypesDepGeneral::Neutral)
                solFile << "NB_INTERVALS;SCALING;";
        }
        solFile << "SOLUTION;";
        solFile << "NEAR_OPT;EPS_NEAR_OPT;";
        solFile << "INSTANCE_FILE;SEED;";
        solFile << "NB_PROBLEMS;NB_SCENARIOS;NB_TEST_SCENARIOS;";
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
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::TypesDepGeneral & type_dep_general) {
    switch(type_dep_general) {
        case Input::TypesDepGeneral::GenProportional: {
            lhs << "Proportional";
            break;
        }
        case Input::TypesDepGeneral::StrongFragile: {
            lhs << "StrongFragile";
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