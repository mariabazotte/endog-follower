#include "instance.hpp"

Instance::Instance(const Input & input): input(input), instance_file(input.getInstanceFile()), 
                                    nb_saa_problems(input.getNbProblemsSAA()), nb_saa_scenarios(input.getNbScenariosSAA()), nb_saa_thinning(input.getNbThinningSAA()),
                                    nb_val_problems(input.getNbValidateProblems()), nb_val_scenarios(input.getNbValidateScenarios()), nb_val_thinning(input.getNbValidateThinning()),
                                    uniform_eps(1e-12,1.0-1e-12), uniform(0.0,1.0), normal(0.0,1.0),
                                    coordinate_mcmc(input.coordinateHitAndRun()), constr_normal_transform(false), polyround_transform(true), 
                                    polyround_transform_specific(true), transform_specific_executed(false), metropolis_acceptance(input.metropolisHastings()){
    
    defineName();
    defineCriticalValues(); 

    for(int i = 0; i < std::max(nb_saa_problems,nb_val_problems); ++i)
        rd_engine.push_back(std::mt19937(1+7919*i));
    
    // Read instance.
    bilevelmodel = new BilevelModel(instance_file);
    uniform_int = std::uniform_int_distribution<int>(0,bilevelmodel->nb_follower_vars-1);
    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral) defineDepGenParams();
    defineFollowerProblemEigenMatrix();
    std::cout << "Read instance finished..." << std::endl;
    
    double time_ = time(NULL);
    generateScenarios();
    std::cout << "Time to generate scenarios: " << time(NULL) - time_ << std::endl;
    std::cout << "Generate scenarios finished..." << std::endl;
}

void Instance::defineName(){
    name = instance_file;
    while(name.find("/") != std::string::npos)
        name = name.substr(name.find("/")+1);
    name.substr(0,name.size()-4);
} 

void Instance::defineCriticalValues(){
    critical_normal = 1.64;
    if(nb_saa_problems == 1) critical_tstudent = 6.31;
    if(nb_saa_problems == 5) critical_tstudent = 2.02;
    if(nb_saa_problems == 10) critical_tstudent = 1.81;
    else if(nb_saa_problems == 11) critical_tstudent = 1.80;
    else if(nb_saa_problems == 12) critical_tstudent = 1.78;
    else if(nb_saa_problems == 13) critical_tstudent = 1.77;
    else if(nb_saa_problems == 14) critical_tstudent = 1.76;
    else if(nb_saa_problems == 15 || nb_saa_problems == 16) critical_tstudent = 1.75;
    else if(nb_saa_problems == 17) critical_tstudent = 1.74;
    else if(nb_saa_problems == 18 || nb_saa_problems == 19) critical_tstudent = 1.73;
    else if(nb_saa_problems >= 20 && nb_saa_problems <= 22) critical_tstudent = 1.72;
    else if(nb_saa_problems >= 23 && nb_saa_problems <= 26) critical_tstudent = 1.71;
    else if(nb_saa_problems >= 27 && nb_saa_problems <= 39) critical_tstudent = 1.70;
    else if(nb_saa_problems >= 40 && nb_saa_problems <= 59) critical_tstudent = 1.68;
    else if(nb_saa_problems >= 60 && nb_saa_problems <= 79) critical_tstudent = 1.67;
    else if(nb_saa_problems >= 80 && nb_saa_problems <= 199) critical_tstudent = 1.66;
}

// TODO REMOVE PI_C
void Instance::defineDepGenParams(){
    double Fc = (bilevelmodel->leader_lb + bilevelmodel->leader_ub)/2.0;
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional){
        double delta = (bilevelmodel->leader_ub - bilevelmodel->leader_lb)/2.0;
        gen_ref_pt_probab = input.getMaxIntCoeffGeneral()*delta;
        gen_min_probab = gen_ref_pt_probab - input.getMaxIntCoeffGeneral()*(bilevelmodel->leader_ub-Fc);
        gen_max_probab = gen_ref_pt_probab + input.getMaxIntCoeffGeneral()*(bilevelmodel->leader_ub-Fc);
    }if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragile){
        gen_ref_pt_probab = 1.0;
        gen_min_probab = std::exp(-input.getMaxIntCoeffGeneral()*(bilevelmodel->leader_ub-Fc)/(bilevelmodel->leader_ub-bilevelmodel->leader_lb));
        gen_max_probab = std::exp(input.getMaxIntCoeffGeneral()*(bilevelmodel->leader_ub-Fc)/(bilevelmodel->leader_ub-bilevelmodel->leader_lb));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Define A_follower (Eigen::MatrixXd), Aeq_follower (Eigen::MatrixXd), and b_follower (Eigen::VectorXd): 
// follower problem matrix considering only follower variables. A_follower consider only inequality constraints, 
// which are defined in standard form as Ay <= b. Aeq_follower consider only equality constraints.
void Instance::defineFollowerProblemEigenMatrix(){
    A_follower.resize(bilevelmodel->nb_follower_constrs-bilevelmodel->nb_follower_eq_constrs,bilevelmodel->nb_follower_vars);
    b_follower.resize(bilevelmodel->nb_follower_constrs-bilevelmodel->nb_follower_eq_constrs);

    // Define all inequality constraints.
    int j = 0;
    for(int c = 0; c < bilevelmodel->nb_follower_constrs; ++c){
        if(bilevelmodel->follower_constrs[c].sense == '=') continue;

        if(bilevelmodel->follower_constrs[c].sense == '<') b_follower(j) = bilevelmodel->follower_constrs[c].rhs;
        else if(bilevelmodel->follower_constrs[c].sense == '>') b_follower(j) = -bilevelmodel->follower_constrs[c].rhs;
        
        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
            BilevelVariable var = bilevelmodel->follower_vars[i];
            if(bilevelmodel->follower_constrs[c].sense == '<') A_follower(j,i) = bilevelmodel->follower_constrs[c].coeffs[var.id];
            else if(bilevelmodel->follower_constrs[c].sense == '>') A_follower(j,i) = - bilevelmodel->follower_constrs[c].coeffs[var.id];
        }
        ++j;   
    }

    // Include (near)-optimality constraint.
    if(input.isFollowerNearOpt() == true){
        int nb_rows = A_follower.rows();

        A_follower.conservativeResize(nb_rows+1, A_follower.cols()); 
        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
            A_follower(nb_rows,i) = bilevelmodel->follower_vars[i].obj_follower;
        }

        b_follower.conservativeResize(nb_rows+1);
        b_follower(nb_rows) = bilevelmodel->follower_ub;
    }

    // Include lower and upper bounds on variables. 
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if(!(var.lb <= -bilevelmodel->inf)){
            Eigen::VectorXd lb = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);
            lb[i] = -1.0;

            int nb_rows = A_follower.rows();
            A_follower.conservativeResize(nb_rows+1, A_follower.cols());
            A_follower.row(nb_rows) = lb;

            b_follower.conservativeResize(nb_rows+1);
            b_follower(nb_rows) = -var.lb;
        }

        if(!(var.ub >= bilevelmodel->inf)){
            Eigen::VectorXd ub = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);
            ub[i] = 1.0;

            int nb_rows = A_follower.rows();
            A_follower.conservativeResize(nb_rows+1, A_follower.cols());
            A_follower.row(nb_rows) = ub;

            b_follower.conservativeResize(nb_rows+1);
            b_follower(nb_rows) = var.ub;
        }
    }

    // Number of equality constraints follower problem + optimality constraint when follower is optimal.
    int nb_total_eq = bilevelmodel->nb_follower_eq_constrs;
    if(input.isFollowerNearOpt() == false) nb_total_eq += 1;
    
    if(nb_total_eq > 0){
        // Define set of equality constraints (A_eq).
        Aeq_follower.resize(nb_total_eq, bilevelmodel->nb_follower_vars);
        
        int j = 0;
        for(int c = 0; c < bilevelmodel->nb_follower_constrs; ++c){
            if(bilevelmodel->follower_constrs[c].sense == '='){
                for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                    BilevelVariable var = bilevelmodel->follower_vars[i];
                    Aeq_follower(j,i) = bilevelmodel->follower_constrs[c].coeffs[var.id];
                }
                ++j;
            }
        }
        if(input.isFollowerNearOpt() == false){
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                Aeq_follower(j,i) = bilevelmodel->follower_vars[i].obj_follower;
        }
    }
}

// Compute transformation matrix using constraint normal idea for the sampling in the hit-and-run algorithm
// used for General Decision-dependent cases.
void Instance::computeDepGenConstrNormTransform(){
    // Compute auxiliar matrix M = sum_i (a_i a_i^T / ||a_i||^2).
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(bilevelmodel->nb_follower_vars, bilevelmodel->nb_follower_vars);
    for (int i = 0; i < A_follower.rows(); ++i) {
        Eigen::VectorXd ai = A_follower.row(i).transpose();
        double norm2 = ai.squaredNorm();
        if (norm2 > 1e-12) {
            M += (ai * ai.transpose()) / norm2;
        }
    }

    // Compute eigen-decomposition of matrix M.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("Eigen decomposition failed");
    }
 
    Eigen::VectorXd D = es.eigenvalues();   // eigenvalues
    Eigen::MatrixXd U = es.eigenvectors();  // eigenvectors
 
    // Construct  M^(-1/2).
    Eigen::VectorXd D_sqrt = D.array().sqrt();
    Eigen::VectorXd D_inv_sqrt = D_sqrt.array().inverse();
    Eigen::MatrixXd Dinv = D_inv_sqrt.asDiagonal();
    constrnorm_tr = U * Dinv * U.transpose();     // M^(-1/2)
    
    // Construct M^(1/2).
    // Eigen::MatrixXd Dsqrt = D_sqrt.asDiagonal();
    // inv_transform = U * Dsqrt * U.transpose();     // M^(1/2)
}

// Compute transformation matrix using polyround library for the sampling in the hit-and-run algorithm
// used for General Decision-dependent cases.
void Instance::computeDepGenPolyRoundTransform(const Eigen::VectorXd & b){
    py::module pr  = py::module::import("PolyRound.api");
    py::module pc  = py::module::import("PolyRound.mutable_classes.polytope");
    py::module np  = py::module::import("numpy");

    // Make row-major copy of A_follower and contiguous copy of b_follower.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A_row = A_follower;
    Eigen::VectorXd b_copy = b;

    // Convert to NumPy arrays.
    py::array_t<double> A_np({A_row.rows(), A_row.cols()}, A_row.data());
    py::array_t<double> b_np(b_copy.size(), b_copy.data());

    // Create polytope.
    py::object Polytope = pc.attr("Polytope");
    py::object poly = Polytope(A_np, b_np);

    // Simplify and Transform polytope.
    // if(do_simplify) poly = pr.attr("PolyRoundApi").attr("simplify_polytope")(poly);
    // if(do_transform) poly = pr.attr("PolyRoundApi").attr("transform_polytope")(poly);

    // Round polytope.
    try {
        poly = pr.attr("PolyRoundApi").attr("round_polytope")(poly);
    }
    catch(py::error_already_set &e) {
        std::string msg = e.what();
        std::cout << msg << " Proceeding anyway." << std::endl;
    }

    // Extract A_new, b_new, transformation S and shift t (numpy is row-major, eigen is column-major)
    py::array A_new_np = poly.attr("A").attr("values").cast<py::array>();
    py::array b_new_np = poly.attr("b").attr("values").cast<py::array>();
    py::array S_np     = poly.attr("transformation").attr("values").cast<py::array>();
    py::array t_np     = poly.attr("shift").cast<py::array>();

    // Copy matrices on row-major as numpy.
    auto A_new_map = Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>((double*)A_new_np.data(),A_new_np.shape(0),A_new_np.shape(1));
    auto S_map = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>((double*)S_np.data(), S_np.shape(0), S_np.shape(1));
    
    // Converts matrices to column-major safely
    Eigen::MatrixXd A_new = A_new_map;  
    Eigen::MatrixXd S = S_map; 
    
    // Copy vectors.
    Eigen::VectorXd b_new = Eigen::Map<Eigen::VectorXd>((double*)b_new_np.data(),b_new_np.shape(0));
    Eigen::VectorXd t = Eigen::Map<Eigen::VectorXd>((double*)t_np.data(),t_np.shape(0));

    // Define polyround transformation.
    polyround_tr = S;

    // Tests for verifying we obtain correct information.
    Eigen::MatrixXd testA = A_row*S;
    if (!(testA.isApprox(A_new, 1e-12))) 
        throw std::runtime_error("Rounding is incorrect. Test for matrix A.");

    Eigen::VectorXd testb = b_copy - A_row*t;
    if (!(testb.isApprox(b_new, 1e-12))) 
        throw std::runtime_error("Rounding is incorrect. Test for vector b.");

    ////////////////////////////////////////////////////////////////////////////////////

    // Create problem and round it.
    // py::module hp  = py::module::import("hopsy");

    // py::object problem = hp.attr("Problem")(A_np, b_np);
    // problem = hp.attr("round")(problem, "simplify"_a=false);

    // py::array test_A_new_np = problem.attr("A").cast<py::array>();
    // py::array test_b_new_np = problem.attr("b").cast<py::array>();
    // py::array test_S_np     = problem.attr("transformation").cast<py::array>();
    // py::array test_t_np     = problem.attr("shift").cast<py::array>();

    // auto test_A_new_map = Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>((double*)test_A_new_np.data(),test_A_new_np.shape(0),test_A_new_np.shape(1));
    // auto test_S_map = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>((double*)test_S_np.data(), test_S_np.shape(0), test_S_np.shape(1));
    
    // Eigen::MatrixXd test_A_new = test_A_new_map; 
    // Eigen::MatrixXd test_S = test_S_map; 
    
    // Eigen::VectorXd test_b_new = Eigen::Map<Eigen::VectorXd>((double*)test_b_new_np.data(),test_b_new_np.shape(0));
    // Eigen::VectorXd test_t = Eigen::Map<Eigen::VectorXd>((double*)test_t_np.data(),test_t_np.shape(0));

    // if (test_S.isApprox(S, 1e-12)) {  // tolerance = 1e-12
    //     std::cout << "Transformation matrices are approximately equal" << std::endl;
    // } else {
    //     std::cout << "Transformation matrices are different" << std::endl;
    // }

    // if (test_t.isApprox(t, 1e-12)) {  // tolerance = 1e-12
    //     std::cout << "Transformation shift are approximately equal" << std::endl;
    // } else {
    //     std::cout << "Transformation shift are different" << std::endl;
    // }

    // Eigen::MatrixXd other_A = A_row*test_S;
    // if (other_A.isApprox(test_A_new, 1e-12)) {  // tolerance = 1e-12
    //     std::cout << "Matrices are approximately equal" << std::endl;
    // } else {
    //     std::cout << "Matrices are different" << std::endl;
    // }
}

// Compute transformation matrix for follower problem considering specific leader decision x_ using polyround 
// library for the sampling in the hit-and-run algorithm used for General Decision-dependent cases.
void Instance::computeDepGenPolyRoundTransform(const double * x_, double fs_){
    // Compute specific transformation.
    Eigen::VectorXd b = b_follower;
    for(int c = 0; c < bilevelmodel->nb_follower_constrs; ++c){
        for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
            if(bilevelmodel->follower_constrs[c].sense == '<')
                b(c) -= bilevelmodel->follower_constrs[c].coeffs[bilevelmodel->leader_vars[i].id]*x_[i];
            if(bilevelmodel->follower_constrs[c].sense == '>')
                b(c) += bilevelmodel->follower_constrs[c].coeffs[bilevelmodel->leader_vars[i].id]*x_[i];     
        }
    }
    if(input.isFollowerNearOpt() == true){
        if(input.isNearOptMult() == true){
            b(bilevelmodel->nb_follower_constrs) = (1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub;
        }else{
            b(bilevelmodel->nb_follower_constrs) = fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub-bilevelmodel->follower_lb);
        }
    }
    computeDepGenPolyRoundTransform(b);
}

// Compute null space of equality constraints of the follower problem.
void Instance::nullSpaceFollowerEqConstrs(bool evaluation = false){
    // Auxiliar matrix with transformed Aeq.
    Eigen::MatrixXd Aeq_follower_transf = Aeq_follower;

    // Transform the matrix of Aeq_follower. Compute A(eq)T.
    if(evaluation == false || polyround_transform_specific == false){
        if(constr_normal_transform == true) 
            Aeq_follower_transf = Aeq_follower * constrnorm_tr;
        
        if(polyround_transform == true) 
            Aeq_follower_transf = Aeq_follower * polyround_tr;
    }else{
        Aeq_follower_transf = Aeq_follower * polyround_tr;
    }

    // Compute null space of EqConstrs
    Eigen::FullPivLU<Eigen::MatrixXd> lu(Aeq_follower_transf);
    nullspace = lu.kernel();
}

// Update information on step and direction when using specific transformation: transformation for a given leader deicison x_.
void Instance::updInfoEvaluateDepGeneral(const double * x_, double fs_){
    if(polyround_transform_specific == false)
        throw std::runtime_error("Should not modify sampling for evaluation: specific polyround transform is false.");
    if(transform_specific_executed == true)
        throw std::runtime_error("Sampling modfication with specific polyround transform already executed.");

    // Transform
    computeDepGenPolyRoundTransform(x_,fs_);

    // Null space equality constraints
    if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false) 
        nullSpaceFollowerEqConstrs(true);

    // Update direction.
    directionDepGeneralScenarioEvalModify();

    // Update step info
    for(int pr = 0; pr < nb_val_problems; ++pr){
        eval_coeff_obj_follower_alpha[pr].clear();
        eval_coeff_obj_leader_alpha[pr].clear();
        eval_coeff_alpha[pr].clear();
        eval_constrs_alpha[pr].clear();
        eval_bdlbs_alpha[pr].clear();
        eval_bdubs_alpha[pr].clear();

        defineGeneralStepInfo(true,pr,eval_coeff_obj_follower_alpha[pr],eval_coeff_obj_leader_alpha[pr],
                    eval_coeff_alpha[pr],eval_constrs_alpha[pr],eval_bdlbs_alpha[pr],eval_bdubs_alpha[pr]);
    }
    
    // Specific polyround transformation.
    transform_specific_executed = true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::generateScenarios(){
    // Define SAA scenarios for Decision-Dependent Strong-Weak case using SAA method.
    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepStrongWeak &&
        input.getSolverApproach() == Input::SolverApproach::TR_SAA){
        scenarios_dep_strwk.resize(nb_saa_problems);
        for(int pr = 0; pr < nb_saa_problems; ++pr){
            for(int s = 0; s < nb_saa_scenarios; ++s){
                scenarios_dep_strwk[pr].push_back(inverseDepStrongWeakScenario(pr));
            }
            std::sort(scenarios_dep_strwk[pr].begin(),scenarios_dep_strwk[pr].end(),std::greater<double>());
        }
    }

    // Compute transformations for sampling Decision-Dependent General case.
    if(constr_normal_transform == true) 
        computeDepGenConstrNormTransform();

    if(polyround_transform == true) 
        computeDepGenPolyRoundTransform(b_follower);

    // Compute null space for equality constraints for Decision-Dependent General case.
    if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false) 
        nullSpaceFollowerEqConstrs();

    // Define SAA scenarios for Decision-Dependent General 
    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral){
        dir_scenarios_dep_gen.resize(nb_saa_problems);
        step_scenarios_dep_gen.resize(nb_saa_problems);
        accep_scenarios_dep_gen.resize(nb_saa_problems);
        
        coeff_obj_follower_alpha.resize(nb_saa_problems);
        coeff_obj_leader_alpha.resize(nb_saa_problems);
        coeff_alpha.resize(nb_saa_problems);
        constrs_alpha.resize(nb_saa_problems);
        bdlbs_alpha.resize(nb_saa_problems);
        bdubs_alpha.resize(nb_saa_problems);

        for(int pr = 0; pr < nb_saa_problems; ++pr){
            dir_scenarios_dep_gen[pr].resize(nb_saa_scenarios*nb_saa_thinning);
            for(int s = 0; s < nb_saa_scenarios*nb_saa_thinning; ++s){
                directionDepGeneralScenario(pr,dir_scenarios_dep_gen[pr][s],false);
                step_scenarios_dep_gen[pr].push_back(uniform(rd_engine[pr]));
                accep_scenarios_dep_gen[pr].push_back(uniform(rd_engine[pr]));
            }
            defineGeneralStepInfo(false,pr,coeff_obj_follower_alpha[pr],coeff_obj_leader_alpha[pr],
                                coeff_alpha[pr],constrs_alpha[pr],bdlbs_alpha[pr],bdubs_alpha[pr]);
        }
    }

    // Define evaluation scenarios for the Decision-Dependent General case.
    // Decision-Dpendent Strong-Weak case has only 2 endogenous scenarios, that can be enumerated for evaluation.
    dir_eval_dep_gen.resize(nb_val_problems);
    step_eval_dep_gen.resize(nb_val_problems);
    accep_eval_dep_gen.resize(nb_val_problems);
    for(int pr = 0; pr < nb_val_problems; ++pr){
        dir_eval_dep_gen[pr].resize(nb_val_scenarios*nb_val_thinning);
        for(int s = 0; s < nb_val_scenarios*nb_val_thinning; ++s){
            directionDepGeneralScenario(pr,dir_eval_dep_gen[pr][s],true);
            step_eval_dep_gen[pr].push_back(uniform(rd_engine[pr]));
            accep_eval_dep_gen[pr].push_back(uniform(rd_engine[pr]));
        }
    }

    eval_coeff_obj_follower_alpha.resize(nb_val_problems);
    eval_coeff_obj_leader_alpha.resize(nb_val_problems);
    eval_coeff_alpha.resize(nb_val_problems);
    eval_constrs_alpha.resize(nb_val_problems);
    eval_bdlbs_alpha.resize(nb_val_problems);
    eval_bdubs_alpha.resize(nb_val_problems);

    if(polyround_transform_specific == false){
        for(int pr = 0; pr < nb_val_problems; ++pr){
            defineGeneralStepInfo(true,pr,eval_coeff_obj_follower_alpha[pr],eval_coeff_obj_leader_alpha[pr],
                    eval_coeff_alpha[pr],eval_constrs_alpha[pr],eval_bdlbs_alpha[pr],eval_bdubs_alpha[pr]);
        }
    }
}

double Instance::inverseDepStrongWeakScenario(int pr){
    if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional)
        return bilevelmodel->follower_ub - uniform_eps(rd_engine[pr])*(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Threshold)
        return bilevelmodel->follower_lb + (bilevelmodel->follower_ub - bilevelmodel->follower_lb)*(0.5 + (1.0/input.getParamThreshold())*std::log((1.0/uniform_eps(rd_engine[pr]) - 1.0)));
        // return (bilevelmodel->follower_ub + bilevelmodel->follower_lb)/2.0 + (1.0/input.getParamThreshold())*std::log((1.0/uniform_eps(rd_engine[pr]) - 1.0));
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Strong)
        return bilevelmodel->follower_lb + (bilevelmodel->follower_ub - bilevelmodel->follower_lb)*(1.0 + (1.0/input.getParamStrong())*std::log((1.0 - uniform_eps(rd_engine[pr])*(1.0 - std::exp(-input.getParamStrong())))));
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Fragile)
        return bilevelmodel->follower_lb - ((bilevelmodel->follower_ub - bilevelmodel->follower_lb)/input.getParamFragile())*std::log((std::exp(-input.getParamFragile()) + uniform_eps(rd_engine[pr])*(1.0 - std::exp(-input.getParamFragile()))));
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::StrongPower)
        return bilevelmodel->follower_lb + (bilevelmodel->follower_ub - bilevelmodel->follower_lb)*(std::pow((1.0 - uniform_eps(rd_engine[pr])), 1.0/input.getParamStrongPower()));
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::FragilePower)
        return bilevelmodel->follower_lb + (bilevelmodel->follower_ub - bilevelmodel->follower_lb)*(1.0 - std::pow(uniform_eps(rd_engine[pr]), 1.0/input.getParamFragilePower()));
    else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior."));
}

void Instance::directionDepGeneralScenario(int pr, std::vector<double> & dir, bool evaluation){
    // Auxiliar vector used to compute direction.
    Eigen::VectorXd aux = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);

    // Initial sampling.
    if(coordinate_mcmc == true){ 
        // Coordinate hit-and-run.
        int chosen = uniform_int(rd_engine[pr]);
        double direct = uniform(rd_engine[pr]);

        aux = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);
        if(direct <= 0.5) aux(chosen) = 1.0;
        else aux(chosen) = -1.0; 
    }else{ 
        // Full hit-and-run.
        aux = Eigen::VectorXd::NullaryExpr(bilevelmodel->nb_follower_vars, [&](){ return normal(rd_engine[pr]); }); 
    }  

    if(evaluation == false || polyround_transform_specific == false){
        // Project to null space of equality constraints.
        if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false){
            if(nullspace.cols() == 0) {
                aux = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);
                return;
            }else{
                // Sample a random vector in the null space.
                aux = nullspace * aux; 
            }
        }

        // Normalize to unit length
        aux.normalize();

        // Transform the direction: new_d = T*old_d (T is nxn and d is nx1)
        if(constr_normal_transform == true)
            aux = constrnorm_tr * aux;

        if(polyround_transform == true)
            aux = polyround_tr * aux;
    }

    // Pass information to vector 
    dir.assign(aux.data(), aux.data() + aux.size());
}

void Instance::directionDepGeneralScenarioEvalModify(){
    if(polyround_transform_specific == false)
        throw std::runtime_error("Should not modify sampling for evaluation: specific polyround transform is false.");
    if(transform_specific_executed == true)
        throw std::runtime_error("Sampling modfication with specific polyround transform already executed.");


    // Modify sampling direction according to specific polyround transformation.
    for(int pr = 0; pr < nb_val_problems; ++pr){
        for(int s = 0; s < nb_val_scenarios*nb_val_thinning; ++s){
            Eigen::VectorXd tmp = Eigen::VectorXd::Map(dir_eval_dep_gen[pr][s].data(), dir_eval_dep_gen[pr][s].size());

            // Project to null space of equality constraints.
            if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false){
                if(nullspace.cols() == 0) {
                    tmp = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);
                    std::copy(tmp.data(), tmp.data() + tmp.size(), dir_eval_dep_gen[pr][s].begin());
                    continue;
                }else{
                    // Sample a random vector in the null space.
                    tmp = nullspace * tmp; 
                }
            }

            // Normalize to unit length.
            tmp.normalize();

            // Transform the direction.
            tmp = polyround_tr * tmp;

            // Copy to vector
            std::copy(tmp.data(), tmp.data() + tmp.size(), dir_eval_dep_gen[pr][s].begin());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::defineGeneralStepInfo(bool eval, int pr, 
                                    std::vector<double> & coeff_obj_follower,
                                    std::vector<double> & coeff_obj_leader,
                                    std::vector<std::vector<double>> & coeff,
                                    std::vector<std::vector<std::vector<int>>> & constrs, 
                                    std::vector<std::vector<std::vector<int>>> & bdlbs, 
                                    std::vector<std::vector<std::vector<int>>> & bdubs) {
    int num_scenarios = nb_saa_scenarios*nb_saa_thinning;
    if(eval == true) num_scenarios = nb_val_scenarios*nb_val_thinning;
    
    // Coefficients step variables (alpha_min and alpha_max).
    coeff.resize(num_scenarios);
    coeff_obj_follower.resize(num_scenarios);
    coeff_obj_leader.resize(num_scenarios);
    for(int s = 0; s < num_scenarios; ++s)
        coeff[s].resize(bilevelmodel->nb_follower_constrs,0.0);
    
    // Constraints to be used to define alpha_min and alpha_max.
    constrs.resize(num_scenarios);
    bdlbs.resize(num_scenarios);
    bdubs.resize(num_scenarios);
    for(int s = 0; s < num_scenarios; ++s){
        constrs[s].resize(2);
        bdlbs[s].resize(2);
        bdubs[s].resize(2);
    }

    // Define constraints for step problem (to define alpha variables).
    for(int c = 0; c < bilevelmodel->nb_follower_constrs; ++c){
        BilevelConstraint constr = bilevelmodel->follower_constrs[c];

        // Skip equalities: they are respected by the sampling method of the direction.
        if(constr.sense == '=') continue; 

        for(int s = 0; s < num_scenarios; ++s){
            coeff[s][c] = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                BilevelVariable var = bilevelmodel->follower_vars[i];

                // Compute coefficient alpha variable for constraint and scenario.
                if(eval == true) coeff[s][c] += constr.coeffs[var.id]*dir_eval_dep_gen[pr][s][i];
                else coeff[s][c] += constr.coeffs[var.id]*getGeneralDirection(pr,s,i);
            }
            if(constr.sense == '<'){
                if(coeff[s][c] >= 1e-10) constrs[s][1].push_back(c);
                else if(coeff[s][c] <= -1e-10) constrs[s][0].push_back(c);
            }
            if(constr.sense == '>'){
                if(coeff[s][c] >= 1e-10) constrs[s][0].push_back(c);
                else if(coeff[s][c] <= -1e-10) constrs[s][1].push_back(c);
            }
        }
    }

    // Define coefficients considering follower objective.
    for(int s = 0; s < num_scenarios; ++s){
        coeff_obj_follower[s] = 0.0;
        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
            BilevelVariable var = bilevelmodel->follower_vars[i];

            if(eval == true) coeff_obj_follower[s] += var.obj_follower*dir_eval_dep_gen[pr][s][i];
            else coeff_obj_follower[s] += var.obj_follower*getGeneralDirection(pr,s,i);
        }
        // Also consider near optimality constraint in this case
        if(input.isFollowerNearOpt() == true){ 
            // This constraints is <=.
            if(coeff_obj_follower[s] >= 1e-10) constrs[s][1].push_back(bilevelmodel->nb_follower_constrs);
            else if(coeff_obj_follower[s] <= -1e-10) constrs[s][0].push_back(bilevelmodel->nb_follower_constrs);
        }
    }

    // Define coefficients considering leader objective.
    for(int s = 0; s < num_scenarios; ++s){
        coeff_obj_leader[s] = 0.0;
        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
            BilevelVariable var = bilevelmodel->follower_vars[i];

            if(eval == true) coeff_obj_leader[s] += var.obj_leader*dir_eval_dep_gen[pr][s][i];
            else coeff_obj_leader[s] += var.obj_leader*getGeneralDirection(pr,s,i);
        }
    }

    // Define bound constraints for step problem (to define alpha variables).
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        for(int s = 0; s < num_scenarios; ++s){
            double dir = 0.0;
            if(eval == true) dir = dir_eval_dep_gen[pr][s][i];
            else dir = getGeneralDirection(pr,s,i);
            
            if(!(var.lb <= -bilevelmodel->inf)){
                if(dir >= 1e-10) bdlbs[s][0].push_back(i);
                else if(dir <= -1e-10) bdlbs[s][1].push_back(i);
            }
            if(!(var.ub >=  bilevelmodel->inf)){
                if(dir >= 1e-10) bdubs[s][1].push_back(i);
                else if(dir <= -1e-10) bdubs[s][0].push_back(i);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::getStrongWeakProbab(double fs_) const{
    if(input.getFollowerBehavior() == Input::FollowerBehavior::FixStrongWeak) 
        return input.getFixCoopLevel();
    else if(input.getFollowerBehavior() == Input::FollowerBehavior::DepStrongWeak){
        if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional)
            return (bilevelmodel->follower_ub - fs_)/(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Threshold){
            double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
            return 1.0 / (1.0 + std::exp(input.getParamThreshold() * (t - 0.5)));
            // double phi_c = (bilevelmodel->follower_lb + bilevelmodel->follower_ub) / 2.0;
            // return 1.0 / (1.0 + std::exp(input.getParamThreshold() * (fs_ - phi_c)));
        }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Strong){
            double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
            return (1.0 - std::exp(-input.getParamStrong() * (1.0 - t))) / (1.0 - std::exp(-input.getParamStrong()));
        }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Fragile){
            double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
            return (1.0 - (1.0 - std::exp(-input.getParamFragile() * t)) / (1.0 - std::exp(-input.getParamFragile())));
        }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::StrongPower){
            double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
            return (1.0 - std::pow(t,input.getParamStrongPower()));
        }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::FragilePower){
            double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
            return std::pow((1.0 - t),input.getParamFragilePower());
        }else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior."));
    }
    else throw std::runtime_error(std::string("Incorrect choice of Behavior."));
}

double Instance::getEvalDepStrongWeakProbab(double fs_, Input::TypesDepStrongWeak eval_behavior, double param) const {
    if(eval_behavior == Input::TypesDepStrongWeak::Proportional)
        return (bilevelmodel->follower_ub - fs_)/(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    else if(eval_behavior== Input::TypesDepStrongWeak::Threshold){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return 1.0 / (1.0 + std::exp(param * (t - 0.5)));
        // double phi_c = (bilevelmodel->follower_lb + bilevelmodel->follower_ub) / 2.0;
        // return 1.0 / (1.0 + std::exp(param * (fs_ - phi_c)));
    }else if(eval_behavior == Input::TypesDepStrongWeak::Strong){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return (1.0 - std::exp(-param * (1.0 - t))) / (1.0 - std::exp(-param));
    }else if(eval_behavior == Input::TypesDepStrongWeak::Fragile){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return 1.0 - (1.0 - std::exp(-param * t)) / (1.0 - std::exp(-param));
    }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::StrongPower){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return (1.0 - std::pow(t,param));
    }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::FragilePower){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return std::pow((1.0 - t),param);
    }else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior."));
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::defineEvalMinStep(int pr, int s, const double * x_, const double * y_, double fs_) const {
    // Lower bound on step (alpha_min).
    double alpha_min = -getStepSizeInterval();

    // Constraints defining alpha_min.
    for(int c : eval_constrs_alpha[pr][s][0]){
        // Usual follower constraints.
        if(c != bilevelmodel->nb_follower_constrs){
            BilevelConstraint constr = bilevelmodel->follower_constrs[c];

            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
                BilevelVariable var = bilevelmodel->leader_vars[i];
                value += constr.coeffs[var.id]*x_[i]; 
            }
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                BilevelVariable var = bilevelmodel->follower_vars[i];
                value += constr.coeffs[var.id]*y_[i];
            }

            if(constr.sense == '>' && (constr.rhs - value) >= 1e-8){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }
            if(constr.sense == '<' && (constr.rhs - value) <= -1e-8){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }

            double diff = (constr.rhs - value);
            // if(diff >= -1e-8 && diff <= 1e-8) diff = 0.0;

            // Update alpha value.
            alpha_min = std::max(alpha_min, diff/eval_coeff_alpha[pr][s][c]);

            // std::cout << "min \t" << constr.rhs << " " << value << " " << eval_coeff_alpha[pr][s][c] << " " << diff/eval_coeff_alpha[pr][s][c] << std::endl;
        }
        // Near optimality constraint.
        else{
            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                BilevelVariable var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            if(input.isNearOptMult() == true){
                if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-8){
                    throw std::runtime_error("Constraint on near optimality was not respected by current y.");
                }
            }else{
                if((fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value) <= -1e-8){
                    throw std::runtime_error("Constraint on near optimality was not respected by current y.");
                }
            }

            double diff = 0.0;
            if(input.isNearOptMult() == true){
                diff = ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value);
            }else{
                diff = (fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value);
            }
            // if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha min.
            alpha_min = std::max(alpha_min, diff/eval_coeff_obj_follower_alpha[pr][s]);

            // std::cout << "min obj\t" << ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub) << " " << value << " " << eval_coeff_obj_alpha[pr][s] << " " << diff/eval_coeff_obj_alpha[pr][s] << std::endl;
        }
    }

    // Lower bound constraints.
    for(int i: eval_bdlbs_alpha[pr][s][0]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.lb - y_[i]) >= 1e-8){
            throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y.");
        }
        
        // Update alpha min.
        alpha_min = std::max(alpha_min, (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i]);
        
        // std::cout << "min lb\t" << var.lb << " " << y_[i] << " " << dir_eval_dep_gen[pr][s][i] << " " << (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    // Upper bound constraints.
    for(int i: eval_bdubs_alpha[pr][s][0]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.ub - y_[i]) <= -1e-8){
            throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Update alpha min.
        alpha_min = std::max(alpha_min, (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i]);
        
        // std::cout << "min ub\t" << var.ub << " " << y_[i] << " " << dir_eval_dep_gen[pr][s][i] << " " << (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    if(alpha_min >= 1e-8) {
        throw std::runtime_error(std::string("Incorrect min step value: Infeasible current point."));
    }

    return alpha_min;
}

double Instance::defineEvalMaxStep(int pr, int s, const double * x_, const double * y_, double fs_) const {
    // Upper bound on step (alpha_max).
    double alpha_max = getStepSizeInterval();

    // Constraints defining alpha_max.
    for(int c : eval_constrs_alpha[pr][s][1]){
        // Usual follower constraints.
        if(c != bilevelmodel->nb_follower_constrs){
            BilevelConstraint constr = bilevelmodel->follower_constrs[c];

            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
                BilevelVariable var = bilevelmodel->leader_vars[i];
                value += constr.coeffs[var.id]*x_[i]; 
            }
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                BilevelVariable var = bilevelmodel->follower_vars[i];
                value += constr.coeffs[var.id]*y_[i];
            }

            if(constr.sense == '>' && (constr.rhs - value) >= 1e-8){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }
            if(constr.sense == '<' && (constr.rhs - value) <= -1e-8){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }

            double diff = (constr.rhs - value);
            // if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha max.
            alpha_max = std::min(alpha_max, diff/eval_coeff_alpha[pr][s][c]);

            // std::cout << "max \t" << constr.rhs << " " << value << " " << eval_coeff_alpha[pr][s][c] << " " << diff/eval_coeff_alpha[pr][s][c] << std::endl;
        }
        // Near optimality constraint.
        else{
            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                BilevelVariable var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            if(input.isNearOptMult() == true){
                if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-8){
                    throw std::runtime_error("Constraint on near optimality was not respected by current y.");
                }
            }else{
                if((fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value) <= -1e-8){
                    throw std::runtime_error("Constraint on near optimality was not respected by current y.");
                }
            }

            double diff = 0.0;
            if(input.isNearOptMult() == true){
                diff = ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value);
            }else{
                diff = (fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value);
            }
            // if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha_max value.
            alpha_max = std::min(alpha_max, diff/eval_coeff_obj_follower_alpha[pr][s]);

            // std::cout << "max obj\t" << (1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub << " " << value << " " << eval_coeff_obj_alpha[pr][s] << " " << diff/eval_coeff_obj_alpha[pr][s] << std::endl;
        }
    }

    // Lower bound constraints.
    for(int i: eval_bdlbs_alpha[pr][s][1]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.lb - y_[i]) >= 1e-8){
            throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Update alpha max.
        alpha_max = std::min(alpha_max, (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i]);

        // std::cout << "max lb\t" << var.lb << " " << y_[i] << " " << dir_eval_dep_gen[pr][s][i] << " " << (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    // Upper bound constraints.
    for(int i: eval_bdubs_alpha[pr][s][1]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.ub - y_[i]) <= -1e-8){
            throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Define alpha max.
        alpha_max = std::min(alpha_max, (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i]);

        // std::cout << "max ub\t" << var.ub << " " << y_[i] << " " << dir_eval_dep_gen[pr][s][i] << " " << (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    if(alpha_max <= -1e-8) throw std::runtime_error(std::string("Incorrect max step value: Infeasible current point."));

    return alpha_max;
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::sampleEvalAlpha(int pr, int s, double alpha_min, double alpha_max, double y_obj_leader, double F_c, double F_delta, double pi_c, double m, Input::TypesDepGeneral eval_behavior) const {
    // Metropolis-Acceptance samples step (alpha) according to uniform distribution.
    if(metropolis_acceptance == false || eval_behavior == Input::TypesDepGeneral::Neutral){
        return sampleEvalUniformAlpha(pr, s, alpha_min, alpha_max);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenProportional){
        double A = pi_c + m*(y_obj_leader - F_c);
        double B = m*computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data());
        return sampleEvalLinearAlpha(pr, s, alpha_min, alpha_max, A, B);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
        double k = (m*computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data()))/(F_delta);
        return sampleEvalExponentialAlpha(pr, s, alpha_min, alpha_max, k);
    }
    else{
        throw std::runtime_error("Wrong option on Type for General Decision-dependent case.");
        exit(0);
    }
}

double Instance::sampleEvalUniformAlpha(int pr, int s, double alpha_min, double alpha_max) const {
    return alpha_min + step_eval_dep_gen[pr][s] * (alpha_max - alpha_min);
}

double Instance::sampleEvalLinearAlpha(int pr, int s, double alpha_min, double alpha_max, double A, double B) const {   
    if (B == 0) return alpha_min + step_eval_dep_gen[pr][s] * (alpha_max - alpha_min);

    double denom = A * (alpha_max - alpha_min) + 0.5 * B * (alpha_max*alpha_max - alpha_min*alpha_min);
    double discrim = A*A + 2*B*(0.5 * B * alpha_min * alpha_min + A * alpha_min + step_eval_dep_gen[pr][s] * denom);

    if (discrim < 0) throw std::runtime_error("Negative discriminant in linear sampling");

    double sqrt_disc = std::sqrt(discrim);
    double alpha1 = (-A + sqrt_disc)/B;
    double alpha2 = (-A - sqrt_disc)/B;

    if (alpha1 >= alpha_min && alpha1 <= alpha_max) return alpha1;
    if (alpha2 >= alpha_min && alpha2 <= alpha_max) return alpha2;

    throw std::runtime_error("No valid root in [alpha_min, alpha_max]");
}

double Instance::sampleEvalExponentialAlpha(int pr, int s, double alpha_min, double alpha_max, double k) const {
    if (k == 0) return alpha_min + step_eval_dep_gen[pr][s] * (alpha_max - alpha_min);

    double exp_min = std::exp(k * alpha_min);
    double exp_max = std::exp(k * alpha_max);
    return std::log(step_eval_dep_gen[pr][s] * (exp_max - exp_min) + exp_min) / k;
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::getEvalDepGeneralProb(double y_obj_leader, double F_c, double F_delta, double pi_c, double m, Input::TypesDepGeneral eval_behavior) const {
    // if(eval_behavior == Input::TypesDepGeneral::GenProportional)
    //     std::cout << "Proportional: " << (pi_c + m*(y_obj_leader - F_c)) << std::endl;
    // if(eval_behavior == Input::TypesDepGeneral::GenFragile)
    //     std::cout << "Strong-Fragile: " << (pi_c*std::exp(m*(y_obj_leader - F_c)/(F_delta))) << std::endl;

    if(eval_behavior == Input::TypesDepGeneral::GenProportional) return (pi_c + m*(y_obj_leader - F_c));
    else if(eval_behavior == Input::TypesDepGeneral::GenFragile) return (pi_c*std::exp(m*(y_obj_leader - F_c)/(F_delta)));
    else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent General Behavior."));
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::computeLeaderObjLeaderVars(const double * x_) const {
    double x_obj_leader = 0.0;
    for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
        x_obj_leader += bilevelmodel->leader_vars[i].obj_leader*x_[i];
    }
    return x_obj_leader;
}

double Instance::computeLeaderObjFollowerVars(const double * y_) const {
    double y_obj_leader = 0.0;
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        y_obj_leader += bilevelmodel->follower_vars[i].obj_leader*y_[i];
    }
    return y_obj_leader;
}

double Instance::computeFollowerObj(const double * y_) const {
    double obj_follower = 0.0;
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        obj_follower += bilevelmodel->follower_vars[i].obj_follower*y_[i];
    }
    return obj_follower;
}

void Instance::computeParamsDepGeneral(double Fs_, double Fw_, double fs_, Input::TypesDepGeneral eval_behavior, int nb_int, double scaling, double & F_c, double & F_delta, double & pi_c, double & m) const {
    // Centering point for leader objective and probability.
    F_c = (Fs_ + Fw_)/2.0;
    F_delta = Fw_ - Fs_;

    // Coefficient.
    double interval = 1.0 / nb_int;
    double fs_prop = (bilevelmodel->follower_ub - fs_) /
                     (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    int i = std::min(static_cast<int>(fs_prop / interval), nb_int - 1);
    m = scaling * (1.0 - 2.0*((i + 0.5) * interval));

    std::cout << "fsp:" << fs_prop << " int:" << interval << std::endl;
    std::cout << "i:" << i << std::endl;

    // Probability centering point.
    pi_c = 1.0;
    if(eval_behavior == Input::TypesDepGeneral::GenProportional) {
        double max_m = 0.0;
        for(int ii = 0; ii < nb_int; ++ ii) 
            max_m = std::max(max_m, (scaling * (1.0 - 2.0*((ii + 0.5) * interval))) );
        std::cout << "max m:" << max_m << std::endl;
        pi_c = (max_m + 1e-8)*((Fw_ - Fs_)/2.0);  // Include small epsilon 1e-4 to avoid probability zero.
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::evaluateFixStrongWeak(double & eval, double & f_eval, 
                                     double Fs_, double Fw_, double fs_eps_, double fw_eps_,
                                     double cooperation) const {
    eval = cooperation*Fs_ + (1.0-cooperation)*Fw_;
    f_eval = cooperation*fs_eps_ + (1.0-cooperation)*fw_eps_;
}

void Instance::evaluateDepStrongWeak(double & eval, double & f_eval, 
                                     double Fs_, double Fw_, double fs_, double fs_eps_, double fw_eps_,
                                     Input::TypesDepStrongWeak eval_behavior, double param) const {
    double beta = getEvalDepStrongWeakProbab(fs_,eval_behavior,param);
    eval = beta*Fs_ + (1.0-beta)*Fw_;
    f_eval = beta*fs_eps_ + (1.0-beta)*fw_eps_;
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::evaluateDepGeneral(double & eval, double & var_eval, 
                                  double & f_eval, double & f_var_eval,
                                  double & ess_, double & rhat_, 
                                  const double * x_, const double * y_, 
                                  double Fs_, double Fw_, double fs_, 
                                  Input::TypesDepGeneral eval_behavior, int nb_int, double scaling) {

    // Apply specific polyround transformation for leader decision x_ if necessary.
    if(polyround_transform_specific == true && transform_specific_executed == false) 
        updInfoEvaluateDepGeneral(x_,fs_);

    // Compute parameters for current evaluation.
    double F_c, F_delta, pi_c, m;
    computeParamsDepGeneral(Fs_,Fw_,fs_,eval_behavior,nb_int,scaling,F_c,F_delta,pi_c,m);
    std::cout  << "Fc: " << F_c << " Fdelta: " << F_delta << " pi_c: " << pi_c << " m: " << m << std::endl;

    // Initialize samples.
    std::vector<Eigen::MatrixXd> samples(nb_val_problems);

    // Acceptance rate for metropolis option.
    std::vector<double> acceptance(nb_val_problems, 0.0);
    double mean_acceptance = 0.0;

    // Evaluation values.
    eval = 0.0; var_eval = 0.0;
    f_eval = 0.0; f_var_eval = 0.0;
    std::vector<Eigen::VectorXd> y_obj_leader(nb_val_problems);
    std::vector<Eigen::VectorXd> obj_follower(nb_val_problems);

    for(int pr = 0; pr < nb_val_problems; ++pr){
        // Initialize samples.
        samples[pr].resize(nb_val_scenarios, bilevelmodel->nb_follower_vars);  
        y_obj_leader[pr].resize(nb_val_scenarios);
        obj_follower[pr].resize(nb_val_scenarios);
        
        samples[pr].setZero();
        y_obj_leader[pr].setZero();
        obj_follower[pr].setZero();

        // Initialize first current sample.
        Eigen::VectorXd current_sample = Eigen::Map<const Eigen::VectorXd>(y_, bilevelmodel->nb_follower_vars);

        for(int s = 0; s < nb_val_scenarios; ++s){
            for(int t = 0; t < nb_val_thinning; ++t){
        
                // Computing alpha_min and alpha_max.
                double alpha_min = defineEvalMinStep(pr,(s*nb_val_thinning+t),x_,current_sample.data(),fs_);
                double alpha_max = defineEvalMaxStep(pr,(s*nb_val_thinning+t),x_,current_sample.data(),fs_);

                // Sampling alpha according to appropriate distribution. (only consider leader objective with respect to follower variables).
                double alpha = sampleEvalAlpha(pr,(s*nb_val_thinning+t),alpha_min,alpha_max,computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,eval_behavior);

                // std::cout << s << " " << t << ": alpha bound: " << alpha_min << " " << alpha_max << std::endl;
                // std::cout << s << " " << t << ": alpha value: " << alpha << std::endl;

                if(metropolis_acceptance == false || eval_behavior == Input::TypesDepGeneral::Neutral){
                    // New point is always accepted.
                    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
                        current_sample(i) = current_sample(i) + alpha*dir_eval_dep_gen[pr][s*nb_val_thinning+t][i];
                }else{
                    // Compute new point to be tested if accepted or not.
                    Eigen::VectorXd test(bilevelmodel->nb_follower_vars);
                    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                        test(i) = current_sample(i) + alpha*dir_eval_dep_gen[pr][s*nb_val_thinning+t][i];

                    // Compute probability to accept new point. (only consider leader objective with respect to follower variables).
                    double h = getEvalDepGeneralProb(computeLeaderObjFollowerVars(test.data()),F_c,F_delta,pi_c,m,eval_behavior)/
                               getEvalDepGeneralProb(computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,eval_behavior);

                    // std::cout << s << ": Leader obj current: " << computeLeaderObjFollowerVars(current_sample.data()) << " " << getEvalDepGeneralProb(computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,eval_behavior) << std::endl;
                    // std::cout << s << ": Leader obj test: " << computeLeaderObjFollowerVars(test.data()) << " " << getEvalDepGeneralProb(computeLeaderObjFollowerVars(test.data()),F_c,F_delta,pi_c,m,eval_behavior) << std::endl;
                    // std::cout << s << ": Possible acceptance: " << (computeLeaderObjFollowerVars(test.data())-Fs_)/(computeLeaderObjFollowerVars(current_sample.data())-Fs_) << std::endl;
                    // std::cout << s << ": Acceptance: " << accep_eval_dep_gen[pr][s*nb_val_thinning+t] << " " << h << " " <<  (accep_eval_dep_gen[pr][s*nb_val_thinning+t] <= h) << std::endl;
                    
                    // Accept or not new point according to this probability.
                    if(accep_eval_dep_gen[pr][s*nb_val_thinning+t] <= h){
                        acceptance[pr] += 1.0;
                        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                            current_sample(i) = test(i);
                    }
                }
            }

            // Include sample.
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                samples[pr](s,i) = current_sample(i);
            
            // Compute leader and follower objective from this scenario.
            y_obj_leader[pr](s) = computeLeaderObjFollowerVars(current_sample.data());
            obj_follower[pr](s) = computeFollowerObj(current_sample.data());

            // std::cout << s << ": Leader obj: " << y_obj_leader[pr](s) << std::endl;
            // std::cout << s << ": Follower obj: " << obj_follower[pr](s) << std::endl;
        }

        // Metropolis acceptance.
        acceptance[pr] /= (nb_val_scenarios*nb_val_thinning);

        // Compute mean for leader and follower objectives.
        eval += y_obj_leader[pr].mean();
        f_eval += obj_follower[pr].mean();
    }

    mean_acceptance = std::accumulate(acceptance.begin(),acceptance.end(),0.0)/nb_val_problems;
     
    // Compute statistics for leader and follower objectives.
    eval /= nb_val_problems;
    f_eval /= nb_val_problems;

    var_eval = mcse(y_obj_leader, true);                 // mcse(y_obj_leader, true, "bulk", -1)
    f_var_eval = mcse(obj_follower, true);               // mcse(obj_follower, true, "bulk", -1)

    Eigen::VectorXd samples_mcse = mcse(samples,true);   // mcse(samples, true, "bulk", -1)

    double ess_eval = ess(y_obj_leader)/(nb_val_problems*nb_val_scenarios);     // ess(y_obj_leader, "bulk", -1)
    double f_ess_eval = ess(obj_follower)/(nb_val_problems*nb_val_scenarios);   // ess(obj_follower, "bulk", -1)

    Eigen::VectorXd samples_ess = ess(samples,"bulk");                               // ess(samples, "bulk", -1)
    double min_ess_dim = samples_ess.minCoeff()/(nb_val_problems*nb_val_scenarios);  
    double max_ess_dim = samples_ess.maxCoeff()/(nb_val_problems*nb_val_scenarios);

    ess_ = min_ess_dim;

    Eigen::VectorXd samples_rhat = rhat(samples);

    rhat_ = 0.0;
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
        rhat_ += std::pow(samples_rhat(i),2);
    rhat_ = std::sqrt(rhat_);

    // Print.
    std::cout << "Metropolis acceptance: " << mean_acceptance << std::endl;

    std::cout << "Effective Sample Size (per dimension):\n" << samples_ess.transpose() << "\n\n";
    std::cout << "Effective Sample Size (follower):\n" << ess(obj_follower) << "\n\n";
    std::cout << "Effective Sample Size (leader):\n" << ess(y_obj_leader) << "\n\n";

    std::cout << "Markov Chain Variance (per dimension):\n" << samples_mcse.transpose() << "\n\n";
    std::cout << "Markov Chain Variance (follower):\n" << f_var_eval << "\n\n";
    std::cout << "Markov Chain Variance (leader):\n" << var_eval << "\n\n";

    std::cout << "R-hat (per dimension):\n" << samples_rhat.transpose() << "\n\n"; 

    // Testing with baseline computation.
    py::module np = py::module::import("numpy");
    py::module hp  = py::module::import("hopsy");

    // Convert to NumPy arrays.
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> copy_samples = samples[0];
    // py::array_t<double> samples_np({copy_samples.rows(), copy_samples.cols()}, copy_samples.data());
    // samples_np = samples_np.reshape({1, nb_val_scenarios, bilevelmodel->nb_follower_vars});

    py::list samples_list;
    for(int pr = 0; pr < nb_val_problems; ++pr){
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> copy_samples = samples[pr];
        py::array_t<double> sample_np({copy_samples.rows(), copy_samples.cols()}, copy_samples.data());
        samples_list.append(sample_np);
    }
    py::array samples_np_3d = np.attr("stack")(samples_list, 0);

    // ESS
    py::object ess = hp.attr("ess")(samples_np_3d,"method"_a="bulk");
    std::cout << "ESS: " << py::str(ess).cast<std::string>() << std::endl;

    std::cout << "--------------------------------------\n";

    // evaluateDepGeneralLibrary(x_,fs_,F_c,F_delta,pi_c,m,eval_behavior);
}

void Instance::evaluateDepGeneralLibrary(const double * x_, double fs_, double F_c, double F_delta, double pi_c, double m, Input::TypesDepGeneral eval_behavior){
    // Make row-major copy of A_follower and contiguous copy of b_follower.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A = A_follower;
    Eigen::VectorXd b = b_follower;
    for(int c = 0; c < bilevelmodel->nb_follower_constrs; ++c){
        for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
            if(bilevelmodel->follower_constrs[c].sense == '<')
                b(c) -= bilevelmodel->follower_constrs[c].coeffs[bilevelmodel->leader_vars[i].id]*x_[i];
            if(bilevelmodel->follower_constrs[c].sense == '>')
                b(c) += bilevelmodel->follower_constrs[c].coeffs[bilevelmodel->leader_vars[i].id]*x_[i];     
        }
    }
    if(input.isFollowerNearOpt() == true){
        if(input.isNearOptMult() == true){
            b(bilevelmodel->nb_follower_constrs) = (1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub;
        }else{
            b(bilevelmodel->nb_follower_constrs) = fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub-bilevelmodel->follower_lb);
        }
    }

    // Make Eigen::VectorXd of leader objective corresponding to follower variables.
    Eigen::VectorXd c(bilevelmodel->nb_follower_vars);
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        c(i) = m*bilevelmodel->follower_vars[i].obj_leader;
        if(eval_behavior == Input::TypesDepGeneral::GenFragile){
            c(i) = c(i)/(F_delta);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    py::module hp  = py::module::import("hopsy");
    py::module np  = py::module::import("numpy");
    py::module models = py::module::import("models");

    // Convert to NumPy arrays.
    py::array_t<double> A_np({A.rows(), A.cols()}, A.data());
    py::array_t<double> b_np(b.size(), b.data());
    py::array_t<double> c_np(c.size(), c.data());

    // Create problem.
    py::object problem; 
    if(eval_behavior == Input::TypesDepGeneral::Neutral){
        problem = hp.attr("Problem")(A_np, b_np);
    }else if(eval_behavior == Input::TypesDepGeneral::GenProportional){
        double a = pi_c + m*(-F_c);
        py::object py_model = models.attr("LinearModel")(a, c_np);
        problem = hp.attr("Problem")(A_np, b_np,"model"_a=py_model);
    }else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
        double a = m*(-F_c)/(F_delta);
        py::object py_model = models.attr("ExponentialModel")(a, c_np);
        problem = hp.attr("Problem")(A_np, b_np,"model"_a=py_model);
    }

    // Compute chebyshev center original.
    py::object original_point = hp.attr("compute_chebyshev_center")(problem);
    // py::print(original_point);

    // Round problem.
    problem = hp.attr("round")(problem,"simplify"_a=false);
    // problem = hp.attr("round")(problem);

    // Create proposal.
    py::object proposal = hp.attr("UniformHitAndRunProposal");

    // Computing starting point.
    py::object starting_point = hp.attr("compute_chebyshev_center")(problem);

    py::object transformed_point = hp.attr("compute_chebyshev_center")(problem,true);
    // py::print(transformed_point);

    // Build MarkovChain list.
    py::list markov_chains;
    for(int i = 0; i < 4; ++i){
        markov_chains.append(hp.attr("MarkovChain")(problem, proposal, "starting_point"_a = starting_point));
    }

    // RNG list.
    py::list rng;
    for(int i = 0; i < 4; ++i){
        rng.append(hp.attr("RandomNumberGenerator")("seed"_a = i * 42));
    }

    // Sampling.
    py::tuple result = hp.attr("sample")(markov_chains,rng,"n_samples"_a = 10000,"thinning"_a = 50,"n_procs"_a = 1);

    py::object acceptance_rates = result[0];
    py::array samples = result[1].cast<py::array>();

    // Print information.
    std::cout << "Acceptance rates: " << py::str(acceptance_rates).cast<std::string>() << std::endl;

    // Rhat
    py::object rhat = hp.attr("rhat")(samples);
    std::cout << "Rhat: " << py::str(rhat).cast<std::string>() << std::endl;

    // ESS
    py::object ess = hp.attr("ess")(samples);
    std::cout << "ESS: " << py::str(ess).cast<std::string>() << std::endl;

    // Convergence: rhat <= 1.01
    bool verify = false;
    if(verify == true){
        py::object rhat_le = rhat.attr("__le__")(py::float_(1.01));
        if (!py::bool_(rhat_le.attr("all")())) {
            throw std::runtime_error("Convergence test failed (rhat > 1.01)");
        }
        py::object ess_ge = ess.attr("__ge__")(py::int_(400));
        if (!py::bool_(ess_ge.attr("all")())) {
            throw std::runtime_error("ESS test failed (ess < 400)");
        }
    }

    // Verify if returns original or transformed samples.
    py::array chain0_np = samples.attr("__getitem__")(0).cast<py::array>();
    auto chain0_map = Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>((double*)chain0_np.data(),chain0_np.shape(0),chain0_np.shape(1));
    Eigen::MatrixXd chain0 = chain0_map;

    for(int i = 0; i < chain0.rows(); ++i){
        Eigen::VectorXd y = chain0.row(i);
        Eigen::VectorXd slack = b - A_follower*y;

        if (slack.size() > 0 && (slack.array() < 1e-12).any()) {
            std::cout << "(array) found negative\n";
            throw std::runtime_error("Sample that does not respect constraints.");
            exit(0);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::display() const{
    std::cout << "----------- INSTANCE PARAMS -------------" << std::endl;
    std::cout << "NB LEADER VARS      : " << bilevelmodel->nb_leader_vars << std::endl;
    std::cout << "NB FOLLOWER VARS    : " << bilevelmodel->nb_follower_vars << std::endl;
    std::cout << "NB LEADER CONSTRS   : " << bilevelmodel->nb_leader_constrs << std::endl;
    std::cout << "NB FOLLOWER CONSTRS : " << bilevelmodel->nb_follower_constrs << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////