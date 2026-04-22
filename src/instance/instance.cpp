#include "instance.hpp"

Instance::Instance(const Input & input): input(input), instance_file(input.getInstanceFile()), 
                                    nb_saa_problems(input.getNbProblemsSAA()), nb_saa_scenarios(input.getNbScenariosSAA()), nb_saa_thinning(input.getNbThinningSAA()),
                                    nb_val_problems(input.getNbValidateProblems()), nb_val_scenarios(input.getNbValidateScenarios()), nb_val_thinning(input.getNbValidateThinning()),
                                    nb_warm_start_scenarios(1000), uniform_eps(1e-12,1.0-1e-12), uniform(0.0,1.0), normal(0.0,1.0), 
                                    coordinate_mcmc(input.coordinateHitAndRun()), metropolis_acceptance(input.metropolisHastings()), eval_slice_sampling(false),
                                    constr_normal_transform(false), covariance_transform(false), polyround_transform(false), polyround_transform_specific(true){
    defineName();
    defineCriticalValues(); 

    for(int i = 0; i < std::max(nb_saa_problems,nb_val_problems); ++i)
        rd_engine.push_back(std::mt19937(1+7919*i));
    
    // Read instance.
    bilevelmodel = new BilevelModel(instance_file);
    uniform_int = std::uniform_int_distribution<int>(0,bilevelmodel->nb_follower_vars-1);
    defineFollowerProblemEigenMatrix();
    std::cout << "Read instance finished..." << std::endl;
    
    double time_ = time(NULL);
    generateScenarios();
    std::cout << "Time to generate scenarios: " << time(NULL) - time_ << std::endl;
    std::cout << "Generate scenarios finished..." << std::endl;
    
    // Leader decision for polyround specific
    if(polyround_transform_specific == true){
        prev_x_ = new double[bilevelmodel->nb_leader_vars];
        eval_prev_x_ = new double[bilevelmodel->nb_leader_vars];

        for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
            prev_x_[i] = bilevelmodel->leader_vars[i].lb - 1;
            eval_prev_x_[i] = bilevelmodel->leader_vars[i].lb - 1;
        }

        prev_yi_ = new double[bilevelmodel->nb_follower_vars];
    }

    int sum_transform = constr_normal_transform + covariance_transform + polyround_transform + polyround_transform_specific;
    if(sum_transform > 1) throw std::runtime_error("Only one transformation for hit and run is allowed.");
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

////////////////////////////////////////////////////////////////////////////////////////////////////

// Define A_follower (Eigen::MatrixXd), Aeq_follower (Eigen::MatrixXd), and b_follower (Eigen::VectorXd): 
// follower problem matrix considering only follower variables. A_follower and b_follower consider only inequality constraints, 
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
            const BilevelVariable & var = bilevelmodel->follower_vars[i];
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
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

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
                    const BilevelVariable & var = bilevelmodel->follower_vars[i];
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
}

// Compute transformation matrix using polyround library for the sampling in the hit-and-run algorithm
// used for General Decision-dependent cases.
void Instance::computeDepGenPolyRoundTransform(const Eigen::VectorXd & b){
    // Inlude libraries.
    py::module np  = py::module::import("numpy");
    py::module pr  = py::module::import("PolyRound.api");
    py::module pc  = py::module::import("PolyRound.mutable_classes.polytope");
    py::module st  = py::module::import("PolyRound.settings");
    py::module ut  = py::module::import("PolyRound.static_classes.lp_utils");

    // Print output polyround.
    // if(input.getVerbose() > 1) py::scoped_ostream_redirect stream(std::cout,py::module_::import("sys").attr("stdout"));
    
    // Ignore warnings.
    py::module warnings = py::module::import("warnings");
    warnings.attr("simplefilter")("ignore");
    // warnings.attr("filterwarnings")("error");

    // Make row-major copy of A_follower and contiguous copy of b_follower.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A_row = A_follower;
    Eigen::VectorXd b_copy = b;

    // Convert to NumPy arrays.
    py::array_t<double> A_np({A_row.rows(), A_row.cols()}, A_row.data());
    py::array_t<double> b_np(b_copy.size(), b_copy.data());

    delete[] prev_yi_;
    prev_yi_ = NULL;

    try {
        double time_ = time(NULL);
        if(input.getVerbose() > 1) std::cout << "Running polyround..." << std::endl;

        // Create polytope.
        py::object Polytope = pc.attr("Polytope");
        py::object poly = Polytope(A_np, b_np);
        py::object settings = st.attr("PolyRoundSettings")(); //"verbose"_a=true,"sgp"_a=true,"presolve"_a=true);

        // poly = pr.attr("PolyRoundApi").attr("simplify_transform_and_round")(poly,settings);//,"normalize"_a=false);
        
        // Simplify polytope.
        bool do_simplify = true;
        if(do_simplify) poly = pr.attr("PolyRoundApi").attr("simplify_polytope")(poly,settings); //"normalize"_a=false);

        // Extract A and b after simplify.
        bool verify = false;
        py::array A_simp_np;
        py::array b_simp_np;
        if((input.getVerbose() > 1) && verify == true) {
            A_simp_np = poly.attr("A").attr("values").attr("copy")().cast<py::array>(); 
            b_simp_np = poly.attr("b").attr("values").attr("copy")().cast<py::array>();
        }
        
        // Transform polytope.
        bool do_transform = (!poly.attr("S").is_none());
        if(do_transform == true && do_simplify == true) {
            if(input.getVerbose() > 1) std::cout << "Running transform..." << std::endl;
            poly = pr.attr("PolyRoundApi").attr("transform_polytope")(poly,settings);
        }

        // Round polytope.
        poly = pr.attr("PolyRoundApi").attr("round_polytope")(poly,settings);

        // Extract transformation S (numpy is row-major, eigen is column-major)
        py::array tr_np = poly.attr("transformation").attr("values").cast<py::array>();

        // Copy matrices on row-major as numpy.
        auto tr_map = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>((double*)tr_np.data(), tr_np.shape(0), tr_np.shape(1));
        
        // Converts matrices to column-major safely.
        Eigen::MatrixXd tr = tr_map; 
    
        // Define polyround transformation.
        polyround_tr = tr;

        // Center to start sampling.
        py::object cv_fd = ut.attr("ChebyshevFinder");
        py::tuple center = cv_fd.attr("chebyshev_center")(poly, settings);

        py::array y_np      = center[0].cast<py::array>();
        py::array s_np      = poly.attr("shift").cast<py::array>();

        Eigen::VectorXd s   = Eigen::Map<Eigen::VectorXd>((double*)s_np.data(),s_np.shape(0));
        Eigen::VectorXd yt  = Eigen::Map<Eigen::VectorXd>((double*)y_np.data(),y_np.shape(0));

        Eigen::VectorXd yo = tr * yt + s;

        Eigen::VectorXd residual = A_row * yo - b_copy;
        double max_violation = residual.maxCoeff();
        if(max_violation <= 1e-8){
            prev_yi_ = new double[bilevelmodel->nb_follower_vars];
            std::copy(yo.data(), yo.data() + yo.size(), prev_yi_);
        }

        // Time for polyround.
        if(input.getVerbose() > 1) std::cout << "Finishing polyround... Time:" << (time(NULL) - time_) << std::endl;

        // Tests for verifying we obtain correct information.
        if(input.getVerbose() > 1 && verify == true) {
            // Verify Chebyshev Center: Compute residual: A*x - b
            Eigen::VectorXd residual = A_row * yo - b_copy;

            // Maximum violation
            double max_violation = residual.maxCoeff();

            // Number of violated constraints
            int num_violated = (residual.array() >= 1e-8).count();

            // Print info
            std::cout << "Max violation of A*x <= b: " << max_violation << std::endl;
            std::cout << "Number of violated constraints: " << num_violated << std::endl;
            
            // Extract A_new, b_new, and shift t (numpy is row-major, eigen is column-major)
            py::array A_new_np = poly.attr("A").attr("values").cast<py::array>();
            py::array b_new_np = poly.attr("b").attr("values").cast<py::array>();
            py::array t_np     = poly.attr("shift").cast<py::array>();

            // Copy matrices on row-major as numpy.
            auto A_simp_map = Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>((double*)A_simp_np.data(),A_simp_np.shape(0),A_simp_np.shape(1));
            auto A_new_map = Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>((double*)A_new_np.data(),A_new_np.shape(0),A_new_np.shape(1));
            
            Eigen::MatrixXd A_simp = A_simp_map;
            Eigen::MatrixXd A_new = A_new_map;  

            // Copy vectors.
            Eigen::VectorXd b_simp = Eigen::Map<Eigen::VectorXd>((double*)b_simp_np.data(),b_simp_np.shape(0));
            Eigen::VectorXd b_new = Eigen::Map<Eigen::VectorXd>((double*)b_new_np.data(),b_new_np.shape(0));
            Eigen::VectorXd t = Eigen::Map<Eigen::VectorXd>((double*)t_np.data(),t_np.shape(0));

            Eigen::MatrixXd test_A = A_simp * tr;
            if (!(test_A.isApprox(A_new, 1e-12))) 
                throw std::runtime_error("Rounding is incorrect. Test for matrix A.");

            Eigen::VectorXd test_b = (b_simp - A_simp * t);
            if (!(test_b.isApprox(b_new, 1e-12))) 
                throw std::runtime_error("Rounding is incorrect. Test for vector b.");

            std::cout << "Transformation/Rounding tests passed..." << std::endl;
        }
    }
    catch(py::error_already_set &e) {
        std::string msg = e.what();
        if(input.getVerbose() > 1) std::cout << "Error during polyround... " << msg  << std::endl;
        polyround_tr = Eigen::MatrixXd::Identity(bilevelmodel->nb_follower_vars, bilevelmodel->nb_follower_vars);
    }
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
void Instance::nullSpaceFollowerEqConstrs(bool specific = false){
    // Auxiliar matrix with transformed Aeq.
    Eigen::MatrixXd Aeq_follower_transf = Aeq_follower;

    // Transform the matrix of Aeq_follower. Compute A(eq)T.
    if(specific == false){
        if(constr_normal_transform == true) 
            Aeq_follower_transf = Aeq_follower * constrnorm_tr;
        
        if(polyround_transform == true) 
            Aeq_follower_transf = Aeq_follower * polyround_tr;
        
    }else{
        if(polyround_transform_specific == true)
            Aeq_follower_transf = Aeq_follower * polyround_tr;
        
        if(covariance_transform == true)
            Aeq_follower_transf = Aeq_follower * covariance_tr;
    }

    // Compute null space of EqConstrs
    Eigen::FullPivLU<Eigen::MatrixXd> lu(Aeq_follower_transf);
    nullspace = lu.kernel();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Update information on step and direction when using specific transformation: transformation for a given leader deicison x_.
void Instance::updDirDepGeneral(bool evaluation){

    // Null space equality constraints
    if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false) 
        nullSpaceFollowerEqConstrs(true);

    // Update direction.
    transformDirDepGeneralScenario(evaluation);

    // Update step info.
    if(evaluation == true){
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
    }else{
        for(int pr = 0; pr < nb_saa_problems; ++pr){
            coeff_obj_follower_alpha[pr].clear();
            coeff_obj_leader_alpha[pr].clear();
            coeff_alpha[pr].clear();
            constrs_alpha[pr].clear();
            bdlbs_alpha[pr].clear();
            bdubs_alpha[pr].clear();

            defineGeneralStepInfo(false,pr,coeff_obj_follower_alpha[pr],coeff_obj_leader_alpha[pr],
                        coeff_alpha[pr],constrs_alpha[pr],bdlbs_alpha[pr],bdubs_alpha[pr]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::generateScenarios(){
    // Define SAA scenarios for Decision-Dependent Strong-Weak case using SAA method.
    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepStrongWeak &&
        input.getSolverApproach() == Input::SolverApproach::TR_SAA){
        scenarios_dep_strwk.resize(nb_saa_problems);
        weight_dep_strwk.resize(nb_saa_problems);
        for(int pr = 0; pr < nb_saa_problems; ++pr){
            for(int s = 0; s < nb_saa_scenarios; ++s){
                scenarios_dep_strwk[pr].push_back(inverseDepStrongWeakScenario(pr));   
            }
            for(int s = 0; s < nb_saa_scenarios; ++s){
                if(scenarios_dep_strwk[pr][s] <= bilevelmodel->follower_lb - 1e-12) scenarios_dep_strwk[pr][s] = bilevelmodel->follower_lb;
                if(scenarios_dep_strwk[pr][s] >= bilevelmodel->follower_ub + 1e-12) scenarios_dep_strwk[pr][s] = bilevelmodel->follower_ub;
            }
            std::sort(scenarios_dep_strwk[pr].begin(),scenarios_dep_strwk[pr].end(),std::greater<double>());

            std::vector<double> values;
            std::vector<double> weights;

            int s = 0; double w = 0.0;
            while(s < nb_saa_scenarios) {
                double base_val = scenarios_dep_strwk[pr][s];
                values.push_back(scenarios_dep_strwk[pr][s]);
                weights.push_back(w+1.0);

                int j = s + 1;
                if(j < nb_saa_scenarios){
                    while(std::abs(base_val - scenarios_dep_strwk[pr][j]) <= 0.01) {
                        weights[weights.size()-1] += 1.0;
                        j++;
                        if(j >= nb_saa_scenarios) break;
                    }
                }
                w = weights[weights.size()-1]; 
                s = j; 
            }

            nb_scenarios_dep_strwk.push_back(values.size());
            scenarios_dep_strwk[pr].clear();
            scenarios_dep_strwk[pr].resize(values.size());
            std::copy(values.begin(), values.end(), scenarios_dep_strwk[pr].begin());
            weight_dep_strwk[pr].resize(values.size());
            std::copy(weights.begin(), weights.end(), weight_dep_strwk[pr].begin());
        }
    }

    // Compute transformations for sampling Decision-Dependent General case.
    if(constr_normal_transform == true){
        computeDepGenConstrNormTransform();
    }

    if(polyround_transform == true) {
        computeDepGenPolyRoundTransform(b_follower);
    }

    // Compute null space for equality constraints for Decision-Dependent General case.
    if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false) 
        nullSpaceFollowerEqConstrs(false);

    // Define SAA scenarios for Decision-Dependent General 
    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral){
        dir_scenarios_dep_gen.resize(nb_saa_problems);
        step_scenarios_dep_gen.resize(nb_saa_problems);
        accep_scenarios_dep_gen.resize(nb_saa_problems);

        if((input.getMethodDepGeneral() == Input::MethodDepGeneral::LocalSearch || 
            input.getMethodDepGeneral() == Input::MethodDepGeneral::Callback)
            && (covariance_transform == true || polyround_transform_specific == true)){
            dir_scenarios_dep_gen_base.resize(nb_saa_scenarios*nb_saa_thinning);
        }
        
        coeff_obj_follower_alpha.resize(nb_saa_problems);
        coeff_obj_leader_alpha.resize(nb_saa_problems);
        coeff_alpha.resize(nb_saa_problems);
        constrs_alpha.resize(nb_saa_problems);
        bdlbs_alpha.resize(nb_saa_problems);
        bdubs_alpha.resize(nb_saa_problems);

        for(int pr = 0; pr < nb_saa_problems; ++pr){

            dir_scenarios_dep_gen[pr].resize(nb_saa_scenarios*nb_saa_thinning);

            if((input.getMethodDepGeneral() == Input::MethodDepGeneral::LocalSearch || 
                input.getMethodDepGeneral() == Input::MethodDepGeneral::Callback)
                && (covariance_transform == true || polyround_transform_specific == true)){
                dir_scenarios_dep_gen_base[pr].resize(nb_saa_scenarios*nb_saa_thinning);
            }

            for(int s = 0; s < nb_saa_scenarios*nb_saa_thinning; ++s){

                if((input.getMethodDepGeneral() == Input::MethodDepGeneral::LocalSearch || 
                    input.getMethodDepGeneral() == Input::MethodDepGeneral::Callback)
                    && (covariance_transform == true || polyround_transform_specific == true)){
                    dirDepGeneralScenario(pr,dir_scenarios_dep_gen_base[pr][s]);
                    dir_scenarios_dep_gen[pr][s].resize(bilevelmodel->nb_follower_vars);
                }else normDirDepGeneralScenario(pr,dir_scenarios_dep_gen[pr][s]);

                step_scenarios_dep_gen[pr].push_back(std::round(uniform(rd_engine[pr])*1000.0)/1000.0);
                accep_scenarios_dep_gen[pr].push_back(std::round(uniform(rd_engine[pr])*1000.0)/1000.0);
            }

            if(!((input.getMethodDepGeneral() == Input::MethodDepGeneral::LocalSearch || 
                input.getMethodDepGeneral() == Input::MethodDepGeneral::Callback)
                && (covariance_transform == true || polyround_transform_specific == true))){
                defineGeneralStepInfo(false,pr,coeff_obj_follower_alpha[pr],coeff_obj_leader_alpha[pr],
                        coeff_alpha[pr],constrs_alpha[pr],bdlbs_alpha[pr],bdubs_alpha[pr]);
            }
        }
    }

    // Define evaluation scenarios for the Decision-Dependent General case.
    // Decision-Dpendent Strong-Weak case has only 2 endogenous scenarios, 
    // that can be enumerated for evaluation.

    // Define samples.
    dir_eval_dep_gen.resize(nb_val_problems);
    step_eval_dep_gen.resize(nb_val_problems);
    accep_eval_dep_gen.resize(nb_val_problems);

    // Save base sample direction to be transformed and normalized according to a fixed leader decision.
    if((covariance_transform == true) || (polyround_transform_specific == true)) 
        dir_eval_dep_gen_base.resize(nb_val_problems); 

    for(int pr = 0; pr < nb_val_problems; ++pr){
        dir_eval_dep_gen[pr].resize(nb_val_scenarios*nb_val_thinning);
        
        if((covariance_transform == true) || (polyround_transform_specific == true)) 
            dir_eval_dep_gen_base[pr].resize(nb_val_scenarios*nb_val_thinning);

        for(int s = 0; s < nb_val_scenarios*nb_val_thinning; ++s){
            if((covariance_transform == true) || (polyround_transform_specific == true)){ 
                dirDepGeneralScenario(pr,dir_eval_dep_gen_base[pr][s]);
                dir_eval_dep_gen[pr][s].resize(bilevelmodel->nb_follower_vars);
            }else normDirDepGeneralScenario(pr,dir_eval_dep_gen[pr][s]);

            step_eval_dep_gen[pr].push_back(uniform(rd_engine[pr]));
            accep_eval_dep_gen[pr].push_back(uniform(rd_engine[pr]));
        }
    }

    // Define auxiliar information for alpha computation.
    eval_coeff_obj_follower_alpha.resize(nb_val_problems);
    eval_coeff_obj_leader_alpha.resize(nb_val_problems);
    eval_coeff_alpha.resize(nb_val_problems);
    eval_constrs_alpha.resize(nb_val_problems);
    eval_bdlbs_alpha.resize(nb_val_problems);
    eval_bdubs_alpha.resize(nb_val_problems);

    // Only prcomputed if samples are not transformed according to specific leader decision.
    if((covariance_transform == false) && (polyround_transform_specific == false)){
        for(int pr = 0; pr < nb_val_problems; ++pr){
            defineGeneralStepInfo(true,pr,eval_coeff_obj_follower_alpha[pr],eval_coeff_obj_leader_alpha[pr],
                    eval_coeff_alpha[pr],eval_constrs_alpha[pr],eval_bdlbs_alpha[pr],eval_bdubs_alpha[pr]);
        }
    }

    // Create new chain for warm start, which is used for computing covariance and specific transformation
    // for a given leader decision.
    if(covariance_transform == true) {
        dir_eval_dep_gen.resize(nb_val_problems+1);
        step_eval_dep_gen.resize(nb_val_problems+1);
        accep_eval_dep_gen.resize(nb_val_problems+1);

        dir_eval_dep_gen[nb_val_problems].resize(nb_warm_start_scenarios);
        for(int s = 0; s < nb_warm_start_scenarios; ++s){
            normDirDepGeneralScenario(0,dir_eval_dep_gen[nb_val_problems][s]);
            step_eval_dep_gen[nb_val_problems].push_back(uniform(rd_engine[0]));
            accep_eval_dep_gen[nb_val_problems].push_back(uniform(rd_engine[0]));
        }
        
        eval_coeff_obj_follower_alpha.resize(nb_val_problems+1);
        eval_coeff_obj_leader_alpha.resize(nb_val_problems+1);
        eval_coeff_alpha.resize(nb_val_problems+1);
        eval_constrs_alpha.resize(nb_val_problems+1);
        eval_bdlbs_alpha.resize(nb_val_problems+1);
        eval_bdubs_alpha.resize(nb_val_problems+1);

        defineGeneralStepInfo(true,nb_val_problems,eval_coeff_obj_follower_alpha[nb_val_problems],eval_coeff_obj_leader_alpha[nb_val_problems],
                eval_coeff_alpha[nb_val_problems],eval_constrs_alpha[nb_val_problems],eval_bdlbs_alpha[nb_val_problems],eval_bdubs_alpha[nb_val_problems]);
    }
}

double Instance::inverseDepStrongWeakScenario(int pr){
    if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional)
        return bilevelmodel->follower_ub - uniform_eps(rd_engine[pr])*(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Threshold)
        return bilevelmodel->follower_lb + (bilevelmodel->follower_ub - bilevelmodel->follower_lb)*(0.5 + (1.0/input.getParamThreshold())*std::log((1.0/uniform_eps(rd_engine[pr]) - 1.0)));
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

void Instance::dirEigenDepGeneralScenario(int pr, Eigen::VectorXd & aux){
    // Auxiliar vector used to compute direction.
    aux = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);

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

        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
            aux(i) = std::round(aux(i) * 1000.0)/1000.0;
        }
    } 
}

void Instance::dirDepGeneralScenario(int pr, std::vector<double> & dir){
    // Sample using Eigen::VectorXd
    Eigen::VectorXd aux;
    dirEigenDepGeneralScenario(pr,aux);
    
    // Pass information to vector 
    dir.assign(aux.data(), aux.data() + aux.size());
}

void Instance::normDirDepGeneralScenario(int pr, std::vector<double> & dir){
    // Sample using Eigen::VectorXd
    Eigen::VectorXd aux;
    dirEigenDepGeneralScenario(pr,aux);

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

    if(polyround_transform == true){
        if(polyround_tr.cols() < bilevelmodel->nb_follower_vars){
            aux.conservativeResize(polyround_tr.cols());
        }
    }

    // Normalize to unit length
    aux.normalize();

    // Transform the direction: new_d = T*old_d (T is nxn and d is nx1)
    if(constr_normal_transform == true)
        aux = constrnorm_tr * aux;

    if(polyround_transform == true){
        aux = polyround_tr * aux;
    }

    // Pass information to vector 
    dir.assign(aux.data(), aux.data() + aux.size());
}
    
////////////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::transformDirDepGeneralScenario(bool evaluation){
    if(polyround_transform_specific == false && covariance_transform == false)
        throw std::runtime_error("Should not modify sampling for evaluation: specific polyround transform and covariance transform are false.");

    int num_problems = nb_val_problems;
    if(evaluation == false) num_problems = nb_saa_problems;

    int num_scenarios = nb_val_scenarios*nb_val_thinning;
    if(evaluation == false) num_scenarios = nb_saa_scenarios*nb_saa_thinning;

    // Modify sampling direction according to specific polyround transformation.
    for(int pr = 0; pr < num_problems; ++pr){
        for(int s = 0; s < num_scenarios; ++s){
            
            Eigen::VectorXd tmp;
            if(evaluation == true) tmp = Eigen::VectorXd::Map(dir_eval_dep_gen_base[pr][s].data(), dir_eval_dep_gen_base[pr][s].size());
            else tmp = Eigen::VectorXd::Map(dir_scenarios_dep_gen_base[pr][s].data(), dir_scenarios_dep_gen_base[pr][s].size());
            
            // Project to null space of equality constraints.
            if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false){
                if(nullspace.cols() == 0) {
                    tmp = Eigen::VectorXd::Zero(bilevelmodel->nb_follower_vars);
                    if(evaluation == true) std::copy(tmp.data(), tmp.data() + tmp.size(), dir_eval_dep_gen[pr][s].begin());
                    else std::copy(tmp.data(), tmp.data() + tmp.size(), dir_scenarios_dep_gen[pr][s].begin());
                    continue;
                }else{
                    // Sample a random vector in the null space.
                    tmp = nullspace * tmp; 
                }
            }

            if(polyround_transform_specific == true){
                if(polyround_tr.cols() < bilevelmodel->nb_follower_vars){ // tranform for removing redundent equalities was used
                    if(coordinate_mcmc == true){
                        std::uniform_int_distribution<int> uniform_int_tr(0,polyround_tr.cols()-1);
                        int chosen = uniform_int_tr(rd_engine[pr]);
                        double direct = uniform(rd_engine[pr]);
                        tmp = Eigen::VectorXd::Zero(polyround_tr.cols());
                        if(direct <= 0.5) tmp(chosen) = 1.0;
                        else tmp(chosen) = -1.0;
                    }else tmp.conservativeResize(polyround_tr.cols());  
                }
            }

            // Normalize to unit length.
            tmp.normalize();

            // Transform the direction.
            if(polyround_transform_specific == true) tmp = polyround_tr * tmp;
            if(covariance_transform == true) tmp = covariance_tr * tmp;

            // Copy to vector
            if(evaluation == true) std::copy(tmp.data(), tmp.data() + tmp.size(), dir_eval_dep_gen[pr][s].begin());
            else std::copy(tmp.data(), tmp.data() + tmp.size(), dir_scenarios_dep_gen[pr][s].begin());
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
    if(eval == true) {
        if(pr == nb_val_problems) num_scenarios = nb_warm_start_scenarios;
        else num_scenarios = nb_val_scenarios*nb_val_thinning;
    }
    
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
        const BilevelConstraint & constr = bilevelmodel->follower_constrs[c];

        // Skip equalities: they are respected by the sampling method of the direction.
        if(constr.sense == '=') continue; 

        for(int s = 0; s < num_scenarios; ++s){
            coeff[s][c] = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];

                // Compute coefficient alpha variable for constraint and scenario.
                if(eval == true) coeff[s][c] += constr.coeffs.at(var.id)*dir_eval_dep_gen[pr][s][i];
                else coeff[s][c] += constr.coeffs.at(var.id)*getGeneralDirection(pr,s,i);
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
            const BilevelVariable & var = bilevelmodel->follower_vars[i];

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
            const BilevelVariable & var = bilevelmodel->follower_vars[i];

            if(eval == true) coeff_obj_leader[s] += var.obj_leader*dir_eval_dep_gen[pr][s][i];
            else coeff_obj_leader[s] += var.obj_leader*getGeneralDirection(pr,s,i);
        }
    }

    // Define bound constraints for step problem (to define alpha variables).
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

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
    }else if(eval_behavior == Input::TypesDepStrongWeak::Strong){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return (1.0 - std::exp(-param * (1.0 - t))) / (1.0 - std::exp(-param));
    }else if(eval_behavior == Input::TypesDepStrongWeak::Fragile){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return 1.0 - (1.0 - std::exp(-param * t)) / (1.0 - std::exp(-param));
    }else if(eval_behavior == Input::TypesDepStrongWeak::StrongPower){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return (1.0 - std::pow(t,param));
    }else if(eval_behavior == Input::TypesDepStrongWeak::FragilePower){
        double t = (fs_ - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return std::pow((1.0 - t),param);
    }else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior."));
}

std::string Instance::getStrongWeakStringFunction() const {
    if(input.getTypeDepStrongWeak()== Input::TypesDepStrongWeak::Threshold)
        return ":(1 / (1 + exp(" + std::to_string(input.getParamThreshold()) + "*(x - 0.5))))";
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Strong)
        return ":((1 - exp(-" + std::to_string(input.getParamStrong()) + "*(1 - x)))/(1 - exp(-" + std::to_string(input.getParamStrong()) + ")))";
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::StrongPower)
        return ":(1 - x^" + std::to_string(input.getParamStrongPower()) + ")";
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Fragile)
        return ":(1 - (1 - exp(-" + std::to_string(input.getParamFragile()) + "*x))/(1 - exp(-" + std::to_string(input.getParamFragile()) + ")))";
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::FragilePower)
        return ":((1 - x)^" + std::to_string(input.getParamFragilePower()) + ")";
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional)
        throw std::invalid_argument("Proportional case does not need piecewise linear approximation.");
    else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior.")); 
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::defineEvalMinStep(int pr, int s, const double * x_, const double * y_, double fs_) const {
    // Lower bound on step (alpha_min).
    double alpha_min = -getStepSizeInterval();

    // Constraints defining alpha_min.
    for(int c : eval_constrs_alpha[pr][s][0]){
        // Usual follower constraints.
        if(c != bilevelmodel->nb_follower_constrs){
            const BilevelConstraint & constr = bilevelmodel->follower_constrs[c];

            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
                const BilevelVariable & var = bilevelmodel->leader_vars[i];
                value += constr.coeffs.at(var.id)*x_[i]; 
            }
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += constr.coeffs.at(var.id)*y_[i];
            }

            // if(constr.sense == '>' && (constr.rhs - value) >= 1e-8){
            //     throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            // }
            // if(constr.sense == '<' && (constr.rhs - value) <= -1e-8){
            //     throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            // }

            double diff = (constr.rhs - value);
            // if(diff >= -1e-8 && diff <= 1e-8) diff = 0.0;

            // Update alpha value.
            alpha_min = std::max(alpha_min, diff/eval_coeff_alpha[pr][s][c]);

            // std::cout << "min \t" << constr.rhs << " " << value << " " << 
            //                          eval_coeff_alpha[pr][s][c] << " " << 
            //                          diff/eval_coeff_alpha[pr][s][c] << std::endl;
        }
        // Near optimality constraint.
        else{
            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            // if(input.isNearOptMult() == true){
            //     if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-8){
            //         throw std::runtime_error("Constraint on near optimality was not respected by current y.");
            //     }
            // }else{
            //     if((fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value) <= -1e-8){
            //         throw std::runtime_error("Constraint on near optimality was not respected by current y.");
            //     }
            // }

            double diff = 0.0;
            if(input.isNearOptMult() == true){
                diff = ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value);
            }else{
                diff = (fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value);
            }
            // if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha min.
            alpha_min = std::max(alpha_min, diff/eval_coeff_obj_follower_alpha[pr][s]);

            // std::cout << "min obj\t" << ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub) << " " << value << " " << 
            //                             eval_coeff_obj_alpha[pr][s] << " " << diff/eval_coeff_obj_alpha[pr][s] << std::endl;
        }
    }

    // Lower bound constraints.
    for(int i: eval_bdlbs_alpha[pr][s][0]){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        // if((var.lb - y_[i]) >= 1e-8){
        //     throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y.");
        // }
        
        // Update alpha min.
        alpha_min = std::max(alpha_min, (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i]);
        
        // std::cout << "min lb\t" << var.lb << " " << y_[i] << " " << 
        //               dir_eval_dep_gen[pr][s][i] << " " << (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    // Upper bound constraints.
    for(int i: eval_bdubs_alpha[pr][s][0]){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        // if((var.ub - y_[i]) <= -1e-8){
        //     throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        // }

        // Update alpha min.
        alpha_min = std::max(alpha_min, (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i]);
        
        // std::cout << "min ub\t" << var.ub << " " << y_[i] << " " << dir_eval_dep_gen[pr][s][i] << " " << (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    if(alpha_min >= 1e-2) {
        throw std::runtime_error(std::string("Incorrect min step value: Infeasible current point."));
    }
    alpha_min = std::min(alpha_min, 0.0);

    return alpha_min;
}

double Instance::defineEvalMaxStep(int pr, int s, const double * x_, const double * y_, double fs_) const {
    // Upper bound on step (alpha_max).
    double alpha_max = getStepSizeInterval();

    // Constraints defining alpha_max.
    for(int c : eval_constrs_alpha[pr][s][1]){
        // Usual follower constraints.
        if(c != bilevelmodel->nb_follower_constrs){
            const BilevelConstraint & constr = bilevelmodel->follower_constrs[c];

            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
                const BilevelVariable & var = bilevelmodel->leader_vars[i];
                value += constr.coeffs.at(var.id)*x_[i]; 
            }
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += constr.coeffs.at(var.id)*y_[i];
            }

            // if(constr.sense == '>' && (constr.rhs - value) >= 1e-8){
            //     throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            // }
            // if(constr.sense == '<' && (constr.rhs - value) <= -1e-8){
            //     throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            // }

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
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            // if(input.isNearOptMult() == true){
            //     if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-8){
            //         throw std::runtime_error("Constraint on near optimality was not respected by current y.");
            //     }
            // }else{
            //     if((fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value) <= -1e-8){
            //         throw std::runtime_error("Constraint on near optimality was not respected by current y.");
            //     }
            // }

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
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        // if((var.lb - y_[i]) >= 1e-8){
        //     throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y.");
        // }

        // Update alpha max.
        alpha_max = std::min(alpha_max, (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i]);

        // std::cout << "max lb\t" << var.lb << " " << y_[i] << " " << dir_eval_dep_gen[pr][s][i] << " " << (var.lb - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    // Upper bound constraints.
    for(int i: eval_bdubs_alpha[pr][s][1]){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        // if((var.ub - y_[i]) <= -1e-8){
        //     throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        // }

        // Define alpha max.
        alpha_max = std::min(alpha_max, (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i]);

        // std::cout << "max ub\t" << var.ub << " " << y_[i] << " " << dir_eval_dep_gen[pr][s][i] << " " << (var.ub - y_[i])/dir_eval_dep_gen[pr][s][i] << std::endl;
    }

    if(alpha_max <= -1e-2) {
        throw std::runtime_error(std::string("Incorrect max step value: Infeasible current point. " + std::to_string(alpha_max) + "  pr:" + std::to_string(pr) + "  s:" + std::to_string(s)));
    }
    alpha_max = std::max(alpha_max, 0.0);
    
    return alpha_max;
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::defineMinStep(int & j_min, int pr, int s, const double * x_, const double * y_, double fs_, double Fs_, double Fw_, 
                                Input::TypesDepGeneral eval_behavior, int nb_int, double scaling, bool with_slice_sampling) const {
    // Lower bound on step (alpha_min).
    double alpha_min = -getStepSizeInterval();
    j_min = getGeneralNbConstrsAlpha(pr,s,0) + getGeneralNbBdlbsAlpha(pr,s,0) + getGeneralNbBdubsAlpha(pr,s,0);

    // Constraints defining alpha_min.
    int k = 0;
    for(int c : constrs_alpha[pr][s][0]){
        // Usual follower constraints.
        if(c != bilevelmodel->nb_follower_constrs){
            const BilevelConstraint & constr = bilevelmodel->follower_constrs[c];

            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
                const BilevelVariable & var = bilevelmodel->leader_vars[i];
                value += constr.coeffs.at(var.id)*x_[i]; 
            }
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += constr.coeffs.at(var.id)*y_[i];
            }

            if(constr.sense == '>' && (constr.rhs - value) >= 1e-2){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y. " + std::to_string((constr.rhs - value)));
            }
            if(constr.sense == '<' && (constr.rhs - value) <= -1e-2){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y. " + std::to_string((constr.rhs - value)));
            }

            double diff = (constr.rhs - value);
            // if(diff >= -1e-8 && diff <= 1e-8) diff = 0.0;

            // Update alpha value.
            if((diff/coeff_alpha[pr][s][c]) >= alpha_min + 1e-10) j_min = k;
            alpha_min = std::max(alpha_min, diff/coeff_alpha[pr][s][c]);

            // std::cout << "min \t" << constr.rhs << " " << value << " " << 
            //                          coeff_alpha[pr][s][c] << " " << 
            //                          diff/coeff_alpha[pr][s][c] << std::endl;
        }
        // Near optimality constraint.
        else{
            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            if(input.isNearOptMult() == true){
                if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-2){
                    throw std::runtime_error("Constraint on near optimality was not respected by current y.");
                }
            }else{
                if((fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value) <= -1e-2){
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
            if(diff/coeff_obj_follower_alpha[pr][s] >= alpha_min + 1e-10) j_min = k;
            alpha_min = std::max(alpha_min, diff/coeff_obj_follower_alpha[pr][s]);

            // std::cout << "min obj\t" << ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub) << " " << value << " " << 
            //                             coeff_obj_alpha[pr][s] << " " << diff/coeff_obj_alpha[pr][s] << std::endl;
        }
        ++k;
    }

    // Lower bound constraints.
    for(int i: bdlbs_alpha[pr][s][0]){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        if((var.lb - y_[i]) >= 1e-2){
            throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y.");
        }
        
        // Update alpha min.
        if((var.lb - y_[i])/dir_scenarios_dep_gen[pr][s][i] >= alpha_min + 1e-10) j_min = k;
        alpha_min = std::max(alpha_min, (var.lb - y_[i])/dir_scenarios_dep_gen[pr][s][i]);
        
        // std::cout << "min lb\t" << var.lb << " " << y_[i] << " " << 
        //              dir_scenarios_dep_gen[pr][s][i] << " " << (var.lb - y_[i])/dir_scenarios_dep_gen[pr][s][i] << std::endl;
        ++k;
    }

    // Upper bound constraints.
    for(int i: bdubs_alpha[pr][s][0]){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        if((var.ub - y_[i]) <= -1e-2){
            throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Update alpha min.
        if((var.ub - y_[i])/dir_scenarios_dep_gen[pr][s][i] >= alpha_min + 1e-10) j_min = k;
        alpha_min = std::max(alpha_min, (var.ub - y_[i])/dir_scenarios_dep_gen[pr][s][i]);

        // std::cout << "min ub\t" << var.ub << " " << y_[i] << " " << 
        //              dir_scenarios_dep_gen[pr][s][i] << " " << (var.ub - y_[i])/dir_scenarios_dep_gen[pr][s][i] << std::endl;
        ++k;
    }

    // Slice sampling.
    if((with_slice_sampling == true) && (eval_behavior != Input::TypesDepGeneral::Neutral)){
        int index = getGeneralNbConstrsAlpha(pr,s,0) + getGeneralNbBdlbsAlpha(pr,s,0) + getGeneralNbBdubsAlpha(pr,s,0) + 1;
        // Compute parameters for current evaluation.
        double F_c, F_delta, pi_c, m, p;
        computeParamsDepGeneral(Fs_,Fw_,fs_,eval_behavior,nb_int,scaling,F_c,F_delta,pi_c,m,p);

        double y_obj_leader = computeLeaderObjFollowerVars(y_);

        if(eval_behavior == Input::TypesDepGeneral::GenProportional){
            double A = pi_c + (m*(y_obj_leader - F_c))/(F_delta);
            double B = (m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            if(B >= 1e-12){
                if(((accep_scenarios_dep_gen[pr][s] - 1.0)*A)/(B) >= alpha_min + 1e-10) j_min = index;
                alpha_min = std::max(alpha_min, ((accep_scenarios_dep_gen[pr][s] - 1.0)*A)/(B));
            }
        }else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
            double B = (m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            if(B >= 1e-12) {
                if(std::log(accep_scenarios_dep_gen[pr][s])/(B) >= alpha_min + 1e-10) j_min = index;
                alpha_min = std::max(alpha_min, std::log(accep_scenarios_dep_gen[pr][s])/(B));
            }
        }else if(eval_behavior == Input::TypesDepGeneral::GenFragilePower || eval_behavior == Input::TypesDepGeneral::GenStrongPower){
            double A = 0.0;
            double B = 0.0;
            if(m <= -1e-12){
                A = 0.5 - (y_obj_leader - F_c)/(F_delta);
                B = (-computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            }else{
                A = 0.5 + (y_obj_leader - F_c)/(F_delta);
                B = (+computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            }
            if(B >= 1e-12) {
                if(((std::pow(accep_scenarios_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B) >= alpha_min + 1e-10) j_min = index;
                alpha_min = std::max(alpha_min, ((std::pow(accep_scenarios_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B));
            }
        }
    }

    if(alpha_min >= 1e-2) throw std::runtime_error(std::string("Incorrect min step value: Infeasible current point."));
    alpha_min = std::min(alpha_min, 0.0);

    return alpha_min;
}

double Instance::defineMaxStep(int & j_max, int pr, int s, const double * x_, const double * y_, double fs_, double Fs_, double Fw_, 
                            Input::TypesDepGeneral eval_behavior, int nb_int, double scaling, bool with_slice_sampling) const {
    // Upper bound on step (alpha_max).
    double alpha_max = getStepSizeInterval();
    j_max = getGeneralNbConstrsAlpha(pr,s,1) + getGeneralNbBdlbsAlpha(pr,s,1) + getGeneralNbBdubsAlpha(pr,s,1);

    int k = 0;
    // Constraints defining alpha_max.
    for(int c : constrs_alpha[pr][s][1]){
        // Usual follower constraints.
        if(c != bilevelmodel->nb_follower_constrs){
            const BilevelConstraint & constr = bilevelmodel->follower_constrs[c];

            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
                const BilevelVariable & var = bilevelmodel->leader_vars[i];
                value += constr.coeffs.at(var.id)*x_[i]; 
            }
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += constr.coeffs.at(var.id)*y_[i];
            }

            if(constr.sense == '>' && (constr.rhs - value) >= 1e-2){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }
            if(constr.sense == '<' && (constr.rhs - value) <= -1e-2){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }

            double diff = (constr.rhs - value);
            // if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha max.
            if((diff/coeff_alpha[pr][s][c]) <= alpha_max - 1e-10) j_max = k;
            alpha_max = std::min(alpha_max, diff/coeff_alpha[pr][s][c]);

            // std::cout << "max \t" << constr.rhs << " " << value << " " << 
            //              coeff_alpha[pr][s][c] << " " << diff/coeff_alpha[pr][s][c] << std::endl;
        }
        // Near optimality constraint.
        else{
            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                const BilevelVariable & var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            if(input.isNearOptMult() == true){
                if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-2){
                    throw std::runtime_error("Constraint on near optimality was not respected by current y.");
                }
            }else{
                if((fs_ + input.getEpsNearOpt()*(bilevelmodel->follower_ub - bilevelmodel->follower_lb) - value) <= -1e-2){
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
            if((diff/coeff_obj_follower_alpha[pr][s]) <= alpha_max - 1e-10) j_max = k;
            alpha_max = std::min(alpha_max, diff/coeff_obj_follower_alpha[pr][s]);

            // std::cout << "max obj\t" << (1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub << " " << value << " " << 
            //              coeff_obj_alpha[pr][s] << " " << diff/coeff_obj_alpha[pr][s] << std::endl;
        }
        ++k;
    }

    // Lower bound constraints.
    for(int i: bdlbs_alpha[pr][s][1]){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        if((var.lb - y_[i]) >= 1e-2){
            throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y. " + std::to_string((var.lb - y_[i])));
        }

        // Update alpha max.
        if(((var.lb - y_[i])/dir_scenarios_dep_gen[pr][s][i]) <= alpha_max - 1e-10) j_max = k;
        alpha_max = std::min(alpha_max, (var.lb - y_[i])/dir_scenarios_dep_gen[pr][s][i]);

        // std::cout << "max lb\t" << var.lb << " " << y_[i] << " " << 
        //              dir_scenarios_dep_gen[pr][s][i] << " " << (var.lb - y_[i])/dir_scenarios_dep_gen[pr][s][i] << std::endl;
        ++k;
    }

    // Upper bound constraints.
    for(int i: bdubs_alpha[pr][s][1]){
        const BilevelVariable & var = bilevelmodel->follower_vars[i];

        if((var.ub - y_[i]) <= -1e-2){
            throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Define alpha max.
        if(((var.ub - y_[i])/dir_scenarios_dep_gen[pr][s][i]) <= alpha_max - 1e-10) j_max = k;
        alpha_max = std::min(alpha_max, (var.ub - y_[i])/dir_scenarios_dep_gen[pr][s][i]);

        // std::cout << "max ub\t" << var.ub << " " << y_[i] << " " << 
        //              dir_scenarios_dep_gen[pr][s][i] << " " << (var.ub - y_[i])/dir_scenarios_dep_gen[pr][s][i] << std::endl;
        ++k;
    }

    // Slice sampling.
    if((with_slice_sampling == true) && (eval_behavior != Input::TypesDepGeneral::Neutral)){
        int index = getGeneralNbConstrsAlpha(pr,s,1) + getGeneralNbBdlbsAlpha(pr,s,1) + getGeneralNbBdubsAlpha(pr,s,1) + 1;
        
        // Compute parameters for current evaluation.
        double F_c, F_delta, pi_c, m, p;
        computeParamsDepGeneral(Fs_,Fw_,fs_,eval_behavior,nb_int,scaling,F_c,F_delta,pi_c,m,p);

        double y_obj_leader = computeLeaderObjFollowerVars(y_);

        if(eval_behavior == Input::TypesDepGeneral::GenProportional){
            double A = pi_c + (m*(y_obj_leader - F_c))/(F_delta);
            double B = (m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            if(B <= -1e-12){
                if(((accep_scenarios_dep_gen[pr][s] - 1.0)*A)/(B) <= alpha_max - 1e-10) j_max = index;
                alpha_max = std::min(alpha_max, ((accep_scenarios_dep_gen[pr][s] - 1.0)*A)/(B));
            } 
        }else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
            double B = (m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            if(B <= -1e-12) {
                if(std::log(accep_scenarios_dep_gen[pr][s])/(B) <= alpha_max - 1e-10) j_max = index;
                alpha_max = std::min(alpha_max, std::log(accep_scenarios_dep_gen[pr][s])/(B));
            }
        }else if(eval_behavior == Input::TypesDepGeneral::GenFragilePower || eval_behavior == Input::TypesDepGeneral::GenStrongPower){
            double A = 0.0;
            double B = 0.0;
            if(m <= -1e-12){
                A = 0.5 - (y_obj_leader - F_c)/(F_delta);
                B = (-computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            }else{
                A = 0.5 + (y_obj_leader - F_c)/(F_delta);
                B = (+computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
            }
            if(B <= -1e-12) {
                if(((std::pow(accep_scenarios_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B) <= alpha_max - 1e-10) j_max = index;
                alpha_max = std::min(alpha_max, ((std::pow(accep_scenarios_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B));
            }
        }
    }

    if(alpha_max <= -1e-2) throw std::runtime_error(std::string("Incorrect max step value: Infeasible current point.") );
    alpha_max = std::max(alpha_max, 0.0);
    
    return alpha_max;
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::sampleEvalAlpha(int pr, int s, double alpha_min, double alpha_max, double y_obj_leader, double F_c, double F_delta, double pi_c, double m, double p, Input::TypesDepGeneral eval_behavior) const {
    if((alpha_max - alpha_min) <= 1e-8) return alpha_min;
    if(metropolis_acceptance == true){ // Use metropolis hastings in evaluation.
        return sampleEvalUniformAlpha(pr, s, alpha_min, alpha_max);
    }
    else if(eval_behavior == Input::TypesDepGeneral::Neutral){
        return sampleEvalUniformAlpha(pr, s, alpha_min, alpha_max);
    }
    else if(eval_slice_sampling == true){ // Use slice sampling in evaluation.
        defineEvalSlice(pr, s, alpha_min, alpha_max, y_obj_leader, F_c, F_delta, pi_c, m, p, eval_behavior);
        return sampleEvalUniformAlpha(pr, s, alpha_min, alpha_max);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenProportional){ // Use inverse sampling in evaluation.
        double A = pi_c + m*(y_obj_leader - F_c)/(F_delta);
        double B = m*computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data())/(F_delta);
        return sampleEvalLinearAlpha(pr, s, alpha_min, alpha_max, A, B);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
        double k = (m*computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data()))/(F_delta);
        return sampleEvalExponentialAlpha(pr, s, alpha_min, alpha_max, k);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragilePower){
        int signal = (m <= -1e-12) ? -1 : 1;
        double c0 = 1.0/2.0 + signal * (y_obj_leader - F_c)/(F_delta);
        double c1 = signal * computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data())/(F_delta);
        return sampleEvalPowerAlpha(pr, s, alpha_min, alpha_max, c0, c1, p);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenStrongPower){
        int signal = (m <= -1e-12) ? -1 : 1;
        double c0 = 1.0/2.0 + signal * (y_obj_leader - F_c)/(F_delta);
        double c1 = signal * computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data())/(F_delta);
        return sampleEvalPowerAlpha(pr, s, alpha_min, alpha_max, c0, c1, p);
    }
    else{
        throw std::runtime_error("Wrong option on Type for General Decision-dependent case.");
        exit(0);
    }
}

void Instance::defineEvalSlice(int pr, int s, double & alpha_min, double & alpha_max, double y_obj_leader, double F_c, double F_delta, double pi_c, double m, double p, Input::TypesDepGeneral eval_behavior) const {
    if(eval_behavior == Input::TypesDepGeneral::GenProportional){
        double A = pi_c + (m*(y_obj_leader - F_c))/(F_delta);
        double B = (m*computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data()))/(F_delta);

        if(B >= 1e-12) alpha_min = std::max(alpha_min, ((accep_eval_dep_gen[pr][s] - 1.0)*A)/(B));
        else if(B <= -1e-12) alpha_max = std::min(alpha_max, ((accep_eval_dep_gen[pr][s] - 1.0)*A)/(B));
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
        double B = (m*computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data()))/(F_delta);

        if(B >= 1e-12) alpha_min = std::max(alpha_min, std::log(accep_eval_dep_gen[pr][s])/(B));
        else if(B <= -1e-12) alpha_max = std::min(alpha_max, std::log(accep_eval_dep_gen[pr][s])/(B));
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragilePower || eval_behavior == Input::TypesDepGeneral::GenStrongPower){
        double A = 0.0;
        double B = 0.0;
        if(m <= -1e-12){
            A = 0.5 - (y_obj_leader - F_c)/(F_delta);
            B = (-computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data()))/(F_delta);
        }else{
            A = 0.5 + (y_obj_leader - F_c)/(F_delta);
            B = (+computeLeaderObjFollowerVars(dir_eval_dep_gen[pr][s].data()))/(F_delta);
        }

        if(B >= 1e-12) alpha_min = std::max(alpha_min, ((std::pow(accep_eval_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B));
        else if(B <= -1e-12) alpha_max = std::min(alpha_max, ((std::pow(accep_eval_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B));
    }
    else throw std::runtime_error("Wrong option on Type for General Decision-dependent case for slice sampling.");
}

double Instance::sampleEvalUniformAlpha(int pr, int s, double alpha_min, double alpha_max) const {
    return alpha_min + step_eval_dep_gen[pr][s] * (alpha_max - alpha_min);
}

double Instance::sampleEvalLinearAlpha(int pr, int s, double alpha_min, double alpha_max, double A, double B) const {   
    if (std::abs(B) <= 1e-2) return alpha_min + step_eval_dep_gen[pr][s] * (alpha_max - alpha_min);

    double denom = A * (alpha_max - alpha_min) + 0.5 * B * (alpha_max*alpha_max - alpha_min*alpha_min);
    double discrim = A*A + 2.0*B*(0.5 * B * alpha_min * alpha_min + A * alpha_min + step_eval_dep_gen[pr][s] * denom);

    if (discrim < 0) throw std::runtime_error("Negative discriminant in linear sampling");

    double sqrt_disc = std::sqrt(discrim);
    double alpha1 = (-A + sqrt_disc)/B;
    double alpha2 = (-A - sqrt_disc)/B;

    if (alpha1 >= alpha_min && alpha1 <= alpha_max) return alpha1;
    if (alpha2 >= alpha_min && alpha2 <= alpha_max) return alpha2;

    throw std::runtime_error("No valid root in [alpha_min: " + std::to_string(alpha_min) + ", alpha_max: " + std::to_string(alpha_max) + "]. Roots are " + std::to_string(alpha1) + " " + std::to_string(alpha2));

    // if (discrim < 0) discrim = 0.0;
    // double sqrt_disc = std::sqrt(discrim);
    // double alpha = (-A + std::copysign(std::sqrt(discrim), B)) / B;
    // return std::clamp(alpha, alpha_min, alpha_max);
}

double Instance::sampleEvalExponentialAlpha(int pr, int s, double alpha_min, double alpha_max, double k) const {
    if ((k >= -1e-12) && (k <= 1e-12)) return alpha_min + step_eval_dep_gen[pr][s] * (alpha_max - alpha_min);

    double exp_min = std::exp(k * alpha_min);
    double exp_max = std::exp(k * alpha_max);
    return std::log(step_eval_dep_gen[pr][s] * (exp_max - exp_min) + exp_min) / k;
}

double Instance::sampleEvalPowerAlpha(int pr, int s, double alpha_min, double alpha_max, double c0, double c1, double p) const {
    double g_min = c0 + c1*alpha_min;
    double g_max = c0 + c1*alpha_max;

    if((g_min >= -1e-2) && (g_min <= 1e-12)) g_min = 0.0;
    if((g_max >= -1e-2) && (g_max <= 1e-12)) g_max = 0.0;

    if(g_min <= -1e-2 || g_max <= -1e-2) 
        throw std::runtime_error("Negative value for objective values defined by alpha_min and alpha_max.");

    if (std::abs(c1) <= 1e-12) return alpha_min + step_eval_dep_gen[pr][s] * (alpha_max - alpha_min);
    
    return (std::pow((std::pow(g_min,p+1.0) + step_eval_dep_gen[pr][s] * (std::pow(g_max,p+1.0) - std::pow(g_min,p+1.0))), 1.0/(p+1.0)) - c0)/c1;
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::sampleAlpha(int pr, int s, double alpha_min, double alpha_max, double y_obj_leader, double F_c, double F_delta, double pi_c, double m, double p, Input::TypesDepGeneral eval_behavior) const {
    bool use_slice = true;
    if((alpha_max - alpha_min) <= 1e-6) return alpha_min;
    if(metropolis_acceptance == true){ // Use metropolis hastings in evaluation.
        return sampleUniformAlpha(pr, s, alpha_min, alpha_max);
    }
    else if(eval_behavior == Input::TypesDepGeneral::Neutral){
        return sampleUniformAlpha(pr, s, alpha_min, alpha_max);
    }
    else if(use_slice == true){ // Use slice sampling in evaluation.
        defineSlice(pr, s, alpha_min, alpha_max, y_obj_leader, F_c, F_delta, pi_c, m, p, eval_behavior);
        return sampleUniformAlpha(pr, s, alpha_min, alpha_max);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenProportional){ // Use inverse sampling in evaluation.
        double A = pi_c + m*(y_obj_leader - F_c)/(F_delta);
        double B = m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data())/(F_delta);
        return sampleLinearAlpha(pr, s, alpha_min, alpha_max, A, B);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
        double k = (m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
        return sampleExponentialAlpha(pr, s, alpha_min, alpha_max, k);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragilePower){
        int signal = (m <= -1e-12) ? -1 : 1;
        double c0 = 1.0/2.0 + signal * (y_obj_leader - F_c)/(F_delta);
        double c1 = signal * computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data())/(F_delta);
        return samplePowerAlpha(pr, s, alpha_min, alpha_max, c0, c1, p);
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenStrongPower){
        int signal = (m <= -1e-12) ? -1 : 1;
        double c0 = 1.0/2.0 + signal * (y_obj_leader - F_c)/(F_delta);
        double c1 = signal * computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data())/(F_delta);
        return samplePowerAlpha(pr, s, alpha_min, alpha_max, c0, c1, p);
    }
    else{
        throw std::runtime_error("Wrong option on Type for General Decision-dependent case.");
        exit(0);
    }
}

double Instance::sampleUniformAlpha(int pr, int s, double alpha_min, double alpha_max) const {
    return alpha_min + step_scenarios_dep_gen[pr][s] * (alpha_max - alpha_min);
}

void Instance::defineSlice(int pr, int s, double & alpha_min, double & alpha_max, double y_obj_leader, double F_c, double F_delta, double pi_c, double m, double p, Input::TypesDepGeneral eval_behavior) const {
    if(eval_behavior == Input::TypesDepGeneral::GenProportional){
        double A = pi_c + (m*(y_obj_leader - F_c))/(F_delta);
        double B = (m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);

        if(B >= 1e-12) alpha_min = std::max(alpha_min, ((accep_scenarios_dep_gen[pr][s] - 1.0)*A)/(B));
        else if(B <= -1e-12) alpha_max = std::min(alpha_max, ((accep_scenarios_dep_gen[pr][s] - 1.0)*A)/(B));
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragile){
        double B = (m*computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);

        if(B >= 1e-12) alpha_min = std::max(alpha_min, std::log(accep_scenarios_dep_gen[pr][s])/(B));
        else if(B <= -1e-12) alpha_max = std::min(alpha_max, std::log(accep_scenarios_dep_gen[pr][s])/(B));
    }
    else if(eval_behavior == Input::TypesDepGeneral::GenFragilePower || eval_behavior == Input::TypesDepGeneral::GenStrongPower){
        double A = 0.0;
        double B = 0.0;
        if(m <= -1e-12){
            A = 0.5 - (y_obj_leader - F_c)/(F_delta);
            B = (-computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
        }else{
            A = 0.5 + (y_obj_leader - F_c)/(F_delta);
            B = (+computeLeaderObjFollowerVars(dir_scenarios_dep_gen[pr][s].data()))/(F_delta);
        }

        if(B >= 1e-12) alpha_min = std::max(alpha_min, ((std::pow(accep_scenarios_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B));
        else if(B <= -1e-12) alpha_max = std::min(alpha_max, ((std::pow(accep_scenarios_dep_gen[pr][s], 1.0/p) - 1.0)*A)/(B));
    }
    else throw std::runtime_error("Wrong option on Type for General Decision-dependent case for slice sampling.");
}

double Instance::sampleLinearAlpha(int pr, int s, double alpha_min, double alpha_max, double A, double B) const {   
    if (std::abs(B) <= 1e-2) return alpha_min + step_scenarios_dep_gen[pr][s] * (alpha_max - alpha_min);

    double denom = A * (alpha_max - alpha_min) + 0.5 * B * (alpha_max*alpha_max - alpha_min*alpha_min);
    double discrim = A*A + 2.0*B*(0.5 * B * alpha_min * alpha_min + A * alpha_min + step_scenarios_dep_gen[pr][s] * denom);

    if (discrim < 0) throw std::runtime_error("Negative discriminant in linear sampling");

    double sqrt_disc = std::sqrt(discrim);
    double alpha1 = (-A + sqrt_disc)/B;
    double alpha2 = (-A - sqrt_disc)/B;

    if (alpha1 >= alpha_min && alpha1 <= alpha_max) return alpha1;
    if (alpha2 >= alpha_min && alpha2 <= alpha_max) return alpha2;

    throw std::runtime_error("No valid root in [alpha_min: " + std::to_string(alpha_min) + ", alpha_max: " + std::to_string(alpha_max) + "]. Roots are " + std::to_string(alpha1) + " " + std::to_string(alpha2));

    // if (discrim < 0) discrim = 0.0;
    // double sqrt_disc = std::sqrt(discrim);
    // double alpha = (-A + std::copysign(std::sqrt(discrim), B)) / B;
    // return std::clamp(alpha, alpha_min, alpha_max);
}

double Instance::sampleExponentialAlpha(int pr, int s, double alpha_min, double alpha_max, double k) const {
    if ((k >= -1e-2) && (k <= 1e-2)) return alpha_min + step_scenarios_dep_gen[pr][s] * (alpha_max - alpha_min);

    double exp_min = std::exp(k * alpha_min);
    double exp_max = std::exp(k * alpha_max);
    return std::log(step_scenarios_dep_gen[pr][s] * (exp_max - exp_min) + exp_min) / k;
}

double Instance::samplePowerAlpha(int pr, int s, double alpha_min, double alpha_max, double c0, double c1, double p) const {
    double g_min = c0 + c1*alpha_min;
    double g_max = c0 + c1*alpha_max;

    if((g_min >= -1e-2) && (g_min <= 1e-10)) g_min = 0.0;
    if((g_max >= -1e-2) && (g_max <= 1e-10)) g_max = 0.0;

    if(g_min <= -1e-2 || g_max <= -1e-2) 
        throw std::runtime_error("Negative value for objective values defined by alpha_min and alpha_max.");
    
    if (std::abs(c1) <= 1e-8) return alpha_min + step_scenarios_dep_gen[pr][s] * (alpha_max - alpha_min);

    return (std::pow((std::pow(g_min,p+1.0) + step_scenarios_dep_gen[pr][s] * ( std::pow(g_max,p+1.0) - std::pow(g_min,p+1.0))), 1.0/(p+1.0)) - c0)/c1;
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::getEvalDepGeneralProb(double y_obj_leader, double F_c, double F_delta, double pi_c, double m, double p, Input::TypesDepGeneral eval_behavior) const {
    if(eval_behavior == Input::TypesDepGeneral::GenProportional) return (pi_c + m*(y_obj_leader - F_c)/(F_delta));
    else if(eval_behavior == Input::TypesDepGeneral::GenFragile) return (pi_c*std::exp(m*(y_obj_leader - F_c)/(F_delta)));
    else if(eval_behavior == Input::TypesDepGeneral::GenFragilePower || eval_behavior == Input::TypesDepGeneral::GenStrongPower)  {
        if(m <= -1e-12) return (pi_c * std::pow((1.0/2.0 - (y_obj_leader - F_c)/(F_delta)), p));
        else return (pi_c * std::pow((1.0/2.0 + (y_obj_leader - F_c)/(F_delta)), p));
    }else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent General Behavior."));
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

void Instance::computeParamsDepGeneral(double Fs_, double Fw_, double fs_, Input::TypesDepGeneral eval_behavior, int nb_int, double scaling, double & F_c, double & F_delta, double & pi_c, double & m, double & p) const {
    if(nb_int <= (-1.0 - 1e-8)) 
        throw std::invalid_argument("Number of intervals should be >= -1.");
    if(eval_behavior == Input::TypesDepGeneral::GenProportional){
        if(scaling <= (1.0 - 1e-8) || scaling >= (1.0 + 1e-8))
            throw std::invalid_argument("Scaling parameter for proportional case should be = 1.");
    }
    if(eval_behavior == Input::TypesDepGeneral::GenFragile || eval_behavior == Input::TypesDepGeneral::GenFragilePower ){
        if(scaling <= (- 1e-8))
            throw std::invalid_argument("Scaling parameter for fragile case should be > 0.");
    }
    if(eval_behavior == Input::TypesDepGeneral::GenStrongPower ){
        if(scaling <= (1.0 - 1e-8))
            throw std::invalid_argument("Scaling parameter for strong case should be >= 1.0.");
    }

    // Centering point for leader objective and probability.
    F_c = (Fs_ + Fw_)/2.0;
    F_delta = Fw_ - Fs_;

    // Coefficient.
    if(nb_int <= -1){
        double fs_prop = (fs_ - bilevelmodel->follower_lb) /
                     (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        m = scaling * (2.0 * fs_prop - 1.0);
    }else{
        double interval = 1.0 / nb_int;
        double fs_prop = (bilevelmodel->follower_ub - fs_) /
                        (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        int i = std::min(static_cast<int>(fs_prop / interval), nb_int - 1);
        m = scaling * (1.0 - 2.0*((i + 0.5) * interval));
    }

    // Probability centering point and power value.
    p = 0;
    if(eval_behavior == Input::TypesDepGeneral::GenFragile){
        pi_c = 1.0/std::exp((std::abs(m)/2.0));
    }

    if(eval_behavior == Input::TypesDepGeneral::GenFragilePower){
        pi_c = 1.0;
        p = 1.0 + std::abs(m);
    }

    if(eval_behavior == Input::TypesDepGeneral::GenStrongPower){
        pi_c = 1.0;
        p = std::max(1e-8, std::abs(m/(scaling*scaling)));
    }
 
    if(eval_behavior == Input::TypesDepGeneral::GenProportional) {
        double max_m = 0.0;
        if(nb_int <= -1){
            max_m = scaling * 1.0;
        }else{
            double interval = 1.0 / nb_int;
            for(int ii = 0; ii < nb_int; ++ ii) 
                max_m = std::max(max_m, (scaling * (1.0 - 2.0*((ii + 0.5) * interval))) );
        }   
        pi_c = (max_m + 1e-8)/2.0;  // Include small epsilon 1e-4 to avoid probability zero.
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

void Instance::computeWhitening(const Eigen::MatrixXd & pilot_samples) {
    // Compute mean vector
    Eigen::VectorXd mean = pilot_samples.colwise().mean();

    // Subtract mean from all rows
    Eigen::MatrixXd centered = pilot_samples.rowwise() - mean.transpose();

    // Covariance = (centered^T * centered) / (n-1)
    Eigen::MatrixXd cov = (centered.transpose() * centered) / (pilot_samples.rows() - 1);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd > es(cov);
    Eigen::VectorXd eigvals = es.eigenvalues();
    Eigen::MatrixXd eigvecs = es.eigenvectors();

    // Avoid tiny eigenvalues
    for (int i = 0; i < eigvals.size(); ++i)
        eigvals(i) = std::max(eigvals(i), 1e-10);

    Eigen::MatrixXd Dsqrt     = eigvals.cwiseSqrt().asDiagonal();
    covariance_tr = eigvecs * Dsqrt     * eigvecs.transpose();  // Σ^{1/2} for Option B
}

bool Instance::arrays_are_different(const double* x_){
    for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
        if(std::abs(x_[i] - prev_x_[i]) >= 1e-9) return true; 
    }
    return false; 
}

bool Instance::eval_arrays_are_different(const double* x_){
    for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
        if(std::abs(x_[i] - eval_prev_x_[i]) >= 1e-9) return true; 
    }
    return false; 
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::evaluateDepGeneral(double & eval, double & var_eval, 
                                  double & f_eval, double & f_var_eval,
                                  double & ess_, double & rhat_, 
                                  const double * x_, const double * y_, 
                                  double Fs_, double Fw_, double fs_, 
                                  Input::TypesDepGeneral eval_behavior, int nb_int, double scaling) {

    // Apply specific polyround transformation for leader decision x_ if necessary.
    if(polyround_transform_specific == true && eval_arrays_are_different(x_)){
        // Compute polyround transform
        computeDepGenPolyRoundTransform(x_,fs_);
        // Update samples.
        updDirDepGeneral(true);
        std::copy(x_, x_ + bilevelmodel->nb_leader_vars, eval_prev_x_); // keep x_ to compare with in new calls.
    }
    
    // Compute parameters for current evaluation.
    double F_c, F_delta, pi_c, m, p;
    if(eval_behavior != Input::TypesDepGeneral::Neutral) 
        computeParamsDepGeneral(Fs_,Fw_,fs_,eval_behavior,nb_int,scaling,F_c,F_delta,pi_c,m,p);

    // Warming start to compute covariance.
    if(covariance_transform == true){
        Eigen::MatrixXd pilot_samples(nb_warm_start_scenarios, bilevelmodel->nb_follower_vars);
        Eigen::VectorXd current_sample = Eigen::Map<const Eigen::VectorXd>(y_, bilevelmodel->nb_follower_vars);

        for(int s = 0; s < nb_warm_start_scenarios; ++s){
            double alpha_min = defineEvalMinStep(nb_val_problems,s,x_,current_sample.data(),fs_);
            double alpha_max = defineEvalMaxStep(nb_val_problems,s,x_,current_sample.data(),fs_);
            double alpha = sampleEvalAlpha(nb_val_problems,s,alpha_min,alpha_max,computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,p,eval_behavior);

            // Update sample along precomputed direction
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
                current_sample(i) = current_sample(i) + alpha*dir_eval_dep_gen[nb_val_problems][s][i];

            // Store pilot sample
            pilot_samples.row(s) = current_sample.transpose();
        }
        // Compute covariance transformation
        computeWhitening(pilot_samples);
        // Update samples
        updDirDepGeneral(true);
    }

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

    // Time computation.
    double time_sample = 0.0;
    double time_alpha = 0.0;
    double time_upd = 0.0; 
    double time_upd2 = 0.0;                         

    for(int pr = 0; pr < nb_val_problems; ++pr){
        // Initialize samples.
        samples[pr].resize(nb_val_scenarios, bilevelmodel->nb_follower_vars);  
        y_obj_leader[pr].resize(nb_val_scenarios);
        obj_follower[pr].resize(nb_val_scenarios);
        
        samples[pr].setZero();
        y_obj_leader[pr].setZero();
        obj_follower[pr].setZero();

        // Initialize first current sample.
        Eigen::VectorXd current_sample;
        if(prev_yi_ == NULL){
            current_sample = Eigen::Map<const Eigen::VectorXd>(y_, bilevelmodel->nb_follower_vars);
        }else{
            current_sample = Eigen::Map<const Eigen::VectorXd>(prev_yi_, bilevelmodel->nb_follower_vars);
        }

        for(int s = 0; s < nb_val_scenarios; ++s){
            double time_ = 0.0;
            for(int t = 0; t < nb_val_thinning; ++t){
                // Computing alpha_min and alpha_max.
                time_ = time(NULL); 
                double alpha_min = defineEvalMinStep(pr,(s*nb_val_thinning+t),x_,current_sample.data(),fs_);
                double alpha_max = defineEvalMaxStep(pr,(s*nb_val_thinning+t),x_,current_sample.data(),fs_);
                time_alpha += time(NULL) - time_;
                
                // Sampling alpha according to appropriate distribution. (only consider leader objective with respect to follower variables).
                time_ = time(NULL);
                double alpha = sampleEvalAlpha(pr,(s*nb_val_thinning+t),alpha_min,alpha_max,computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,p,eval_behavior);
                time_sample += time(NULL) - time_;

                if(alpha >= alpha_max) alpha = alpha_max;
                if(alpha <= alpha_min) alpha = alpha_min;

                if(metropolis_acceptance == false){
                    // New point is always accepted.
                    time_ = time(NULL);
                    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
                        current_sample(i) = current_sample(i) + alpha*dir_eval_dep_gen[pr][s*nb_val_thinning+t][i];
                    time_upd += time(NULL) - time_;
                }else if(eval_behavior == Input::TypesDepGeneral::Neutral){
                    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
                        current_sample(i) = current_sample(i) + alpha*dir_eval_dep_gen[pr][s*nb_val_thinning+t][i];
                }else{
                    // Compute new point to be tested if accepted or not.
                    Eigen::VectorXd test(bilevelmodel->nb_follower_vars);
                    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                        test(i) = current_sample(i) + alpha*dir_eval_dep_gen[pr][s*nb_val_thinning+t][i];

                    // Compute probability to accept new point. (only consider leader objective with respect to follower variables).
                    double h = getEvalDepGeneralProb(computeLeaderObjFollowerVars(test.data()),F_c,F_delta,pi_c,m,p,eval_behavior)/
                               getEvalDepGeneralProb(computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,p,eval_behavior);

                    // Accept or not new point according to this probability.
                    if(accep_eval_dep_gen[pr][s*nb_val_thinning+t] <= h){
                        acceptance[pr] += 1.0;
                        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                            current_sample(i) = test(i);
                    }
                }
            }

            // Include sample. Compute leader and follower objective from this scenario.
            time_ = time(NULL);
            samples[pr].row(s) = current_sample.transpose();
            y_obj_leader[pr](s) = computeLeaderObjFollowerVars(current_sample.data());
            obj_follower[pr](s) = computeFollowerObj(current_sample.data());
            time_upd2 += time(NULL) - time_;
        }

        // Metropolis acceptance.
        acceptance[pr] /= (nb_val_scenarios*nb_val_thinning);

        // Compute mean for leader and follower objectives.
        eval += y_obj_leader[pr].mean();
        f_eval += obj_follower[pr].mean();
    }

    double time_init = time(NULL);
    mean_acceptance = std::accumulate(acceptance.begin(),acceptance.end(),0.0)/nb_val_problems;

    // Compute statistics for leader and follower objectives.
    eval /= nb_val_problems;
    f_eval /= nb_val_problems;

    var_eval = mcse(y_obj_leader, true);                 // mcse(y_obj_leader, true, "bulk", -1)
    f_var_eval = mcse(obj_follower, true);               // mcse(obj_follower, true, "bulk", -1)

    Eigen::VectorXd samples_ess = ess(samples,"bulk");   // ess(samples, "bulk", -1)
    ess_ = samples_ess.minCoeff()/(nb_val_problems*nb_val_scenarios);  

    Eigen::VectorXd samples_rhat = rhat(samples);
    rhat_ = samples_rhat.maxCoeff();
    double time_est = time(NULL) - time_init;

    // Print.
    if(input.getVerbose() >= 2){
        std::cout << "Leader Eval value: " << eval << std::endl;
        std::cout << "Follower Eval value: " << f_eval << std::endl;
        
        std::cout << "Time for alphas min/max: " << time_alpha << std::endl;
        std::cout << "Time for sample: " << time_sample << std::endl;
        std::cout << "Time for curr upd: " << time_upd << std::endl;
        std::cout << "Time for final upd: " << time_upd2 << std::endl;
        std::cout << "Time for estimators: " << time_est << std::endl;

        Eigen::VectorXd samples_mcse = mcse(samples,true);   // mcse(samples, true, "bulk", -1)
        double ess_eval = ess(y_obj_leader)/(nb_val_problems*nb_val_scenarios);     // ess(y_obj_leader, "bulk", -1)
        double f_ess_eval = ess(obj_follower)/(nb_val_problems*nb_val_scenarios);   // ess(obj_follower, "bulk", -1)

        std::cout << "--------------------------------------\n";
        std::cout << "Metropolis acceptance: " << mean_acceptance << std::endl;

        std::cout << "Effective Sample Size (per dimension):\n" << samples_ess.transpose() << "\n\n";
        std::cout << "Effective Sample Size (follower):\n" << ess_eval << "\n\n";
        std::cout << "Effective Sample Size (leader):\n" << f_ess_eval << "\n\n";

        std::cout << "Markov Chain Variance (per dimension):\n" << samples_mcse.transpose() << "\n\n";
        std::cout << "Markov Chain Variance (follower):\n" << f_var_eval << "\n\n";
        std::cout << "Markov Chain Variance (leader):\n" << var_eval << "\n\n";

        std::cout << "R-hat (per dimension):\n" << samples_rhat.transpose() << "\n\n"; 
        std::cout << "--------------------------------------\n";
    }

    // if(eval_behavior == Input::TypesDepGeneral::Neutral)
    //     evaluateDepGeneralLibrary(x_, fs_, F_c, F_delta, pi_c, m, eval_behavior);


    ////////////////////////////////////////////////////////////////////////////////////////

    // Computing ESS with baseline hopsy.
    // py::module np = py::module::import("numpy");
    // py::module hp  = py::module::import("hopsy");

    // py::list samples_list;
    // for(int pr = 0; pr < nb_val_problems; ++pr){
    //     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> copy_samples = samples[pr];
    //     py::array_t<double> sample_np({copy_samples.rows(), copy_samples.cols()}, copy_samples.data());
    //     samples_list.append(sample_np);
    // }
    // py::array samples_np_3d = np.attr("stack")(samples_list, 0);

    // py::object ess = hp.attr("ess")(samples_np_3d,"method"_a="bulk");
    // std::cout << "ESS: " << py::str(ess).cast<std::string>() << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::evaluateSAADepGeneral(double & obj_value, int pr, const double * x_, const double * y_, double Fs_, double Fw_, double fs_, 
                                    Input::TypesDepGeneral eval_behavior, int nb_int, double scaling) {
    // Single point. Simulation not needed.
    if(std::abs(Fw_ - Fs_) <= 1e-3) { 
        obj_value = Fs_; 
        return;
    }

    // Apply specific polyround transformation for leader decision x_ if necessary.
    if(polyround_transform_specific == true && arrays_are_different(x_)){
        // Compute polyround transform
        computeDepGenPolyRoundTransform(x_,fs_);
        // Update samples.
        updDirDepGeneral(false); // saa
        std::copy(x_, x_ + bilevelmodel->nb_leader_vars, prev_x_); // keep x_ to compare with in new calls.
    }
    
    // Compute parameters for current evaluation.
    double F_c, F_delta, pi_c, m, p;
    if(eval_behavior != Input::TypesDepGeneral::Neutral) 
        computeParamsDepGeneral(Fs_,Fw_,fs_,eval_behavior,nb_int,scaling,F_c,F_delta,pi_c,m,p);

    // Warming start to compute covariance.
    if(covariance_transform == true){
        Eigen::MatrixXd pilot_samples(nb_warm_start_scenarios, bilevelmodel->nb_follower_vars);
        Eigen::VectorXd current_sample = Eigen::Map<const Eigen::VectorXd>(y_, bilevelmodel->nb_follower_vars);

        for(int s = 0; s < nb_warm_start_scenarios; ++s){
            double alpha_min = defineEvalMinStep(nb_val_problems,s,x_,current_sample.data(),fs_);
            double alpha_max = defineEvalMaxStep(nb_val_problems,s,x_,current_sample.data(),fs_);
            double alpha = sampleEvalAlpha(nb_val_problems,s,alpha_min,alpha_max,computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,p,eval_behavior);

            // Update sample along precomputed direction
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
                current_sample(i) = current_sample(i) + alpha*dir_eval_dep_gen[pr][s][i];

            // Store pilot sample
            pilot_samples.row(s) = current_sample.transpose();
        }
        // Compute covariance transformation
        computeWhitening(pilot_samples);
        // Update samples
        updDirDepGeneral(false); // saa
    }

    // Initialize first current sample.
    Eigen::VectorXd current_sample;
    if(prev_yi_ == NULL) {
        current_sample = Eigen::Map<const Eigen::VectorXd>(y_, bilevelmodel->nb_follower_vars);
    }else{
        current_sample = Eigen::Map<const Eigen::VectorXd>(prev_yi_, bilevelmodel->nb_follower_vars);
    }

    obj_value = 0.0;
    for(int s = 0; s < nb_saa_scenarios; ++s){
        for(int t = 0; t < nb_saa_thinning; ++t){
    
            // Computing alpha_min and alpha_max.
            int j_min, j_max;
            double alpha_min = defineMinStep(j_min,pr,(s*nb_saa_thinning+t),x_,current_sample.data(),fs_,Fs_,Fw_,eval_behavior,nb_int,scaling,false);
            double alpha_max = defineMaxStep(j_max,pr,(s*nb_saa_thinning+t),x_,current_sample.data(),fs_,Fs_,Fw_,eval_behavior,nb_int,scaling,false);

            // Sampling alpha according to appropriate distribution. (only consider leader objective with respect to follower variables).
            double alpha = sampleAlpha(pr,(s*nb_saa_thinning+t),alpha_min,alpha_max,computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,p,eval_behavior);

            if(metropolis_acceptance == false){
                // New point is always accepted.
                for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
                    current_sample(i) = current_sample(i) + alpha*dir_scenarios_dep_gen[pr][(s*nb_saa_thinning+t)][i];
            }else if(eval_behavior == Input::TypesDepGeneral::Neutral){
                for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i) 
                    current_sample(i) = current_sample(i) + alpha*dir_scenarios_dep_gen[pr][(s*nb_saa_thinning+t)][i];
            }else{
                // Compute new point to be tested if accepted or not.
                Eigen::VectorXd test(bilevelmodel->nb_follower_vars);
                for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                    test(i) = current_sample(i) + alpha*dir_scenarios_dep_gen[pr][(s*nb_saa_thinning+t)][i];

                // Compute probability to accept new point. (only consider leader objective with respect to follower variables).
                double h = getEvalDepGeneralProb(computeLeaderObjFollowerVars(test.data()),F_c,F_delta,pi_c,m,p,eval_behavior)/
                           getEvalDepGeneralProb(computeLeaderObjFollowerVars(current_sample.data()),F_c,F_delta,pi_c,m,p,eval_behavior);
                
                // Accept or not new point according to this probability.
                if(accep_scenarios_dep_gen[pr][(s*nb_saa_thinning+t)] <= h){
                    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i)
                        current_sample(i) = test(i);
                }
            }
        }
        obj_value += computeLeaderObjFollowerVars(current_sample.data());
    }
    obj_value = obj_value/(double)nb_saa_scenarios;
}

////////////////////////////////////////////////////////////////////////////////////////////////

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
        if(eval_behavior == Input::TypesDepGeneral::GenFragile || 
           eval_behavior == Input::TypesDepGeneral::GenProportional){
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
        double a = pi_c + m*(-F_c)/F_delta;
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
    problem = hp.attr("round")(problem,"simplify"_a=true);
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

    // Print NumPy array
    // py::print(samples);

    auto buf = samples.request();

    int B = buf.shape[0];
    int R = buf.shape[1];
    int C = buf.shape[2];

    double* ptr = static_cast<double*>(buf.ptr);

    double eval = 0.0; double feval = 0.0;
    for (int b = 0; b < B; ++b) {
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> M(ptr + b * R * C, R, C);
        double e = 0.0; double f = 0.0;
        for (int i = 0; i < M.rows(); ++i) {
            e += computeLeaderObjFollowerVars(M.row(i).data());
            f += computeFollowerObj(M.row(i).data());
        }
        e /= M.rows();
        f /= M.rows();
        eval += e;
        feval += f;
    }
    eval /= B;
    feval /= B;

    std::cout << " Python Eval value: " << eval << std::endl;
    std::cout << " Python FEval value: " << feval << std::endl;

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