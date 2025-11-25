#include "instance.hpp"


Instance::Instance(const Input & input): input(input), instance_file(input.getInstanceFile()),seed(input.getSeed()),
                                    rd_engine(input.getSeed()),uniform_eps(1e-12,1.0-1e-12),uniform(0.0,1.0),normal(0.0,1.0),coordinate_mcmc(true),
                                    nb_saa_problems(input.getNbProblemsSAA()),nb_saa_scenarios(input.getNbScenariosSAA()),
                                    nb_val_scenarios(input.getNbValidateScenarios()){
    defineName();
    defineCriticalValues(); 
    bilevelmodel = new BilevelModel(instance_file);
    verifyInstance();
    defineDepGenParams();
    uniform_int = std::uniform_int_distribution<int>(0, bilevelmodel->nb_follower_vars-1);
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

void Instance::verifyInstance(){
    if(bilevelmodel->leader_lb <= -bilevelmodel->inf || bilevelmodel->leader_ub >= bilevelmodel->inf)
        throw std::runtime_error("");
    if(bilevelmodel->follower_lb <= -bilevelmodel->inf || bilevelmodel->follower_ub >= bilevelmodel->inf)
        throw std::runtime_error("");
// TODO
// problems: unboundness for lower level (cannot use instance at all - include something alerting that)
// unboundness for upper level : not able to use general proportional and define upper bound too for strong fragile case (include a test in c++)
}

void Instance::defineDepGenParams(){
    // Define extra parameters for decision-dependent general case.
    gen_nb_intervals = input.getNbIntervalsGeneral();
    
    gen_intervals.resize(gen_nb_intervals+1);
    double interval = 1.0/(double)gen_nb_intervals;
    gen_intervals[0] = 0.0;
    for(int i = 1; i <= gen_nb_intervals; ++i) 
        gen_intervals[i] = gen_intervals[i-1] + interval;
    for(int i = 0; i < gen_nb_intervals; ++i)
        gen_coeff_intervals.push_back(0.5 - (gen_intervals[i]+gen_intervals[i+1])/2.0);

    double Fc = (bilevelmodel->leader_lb + bilevelmodel->leader_ub)/2.0;
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional){
        double delta = (bilevelmodel->leader_ub - bilevelmodel->leader_lb)/2.0;
        gen_ref_pt_probab = input.getScalingGeneral()*0.5*delta;
        gen_min_probab = gen_ref_pt_probab - input.getScalingGeneral()*0.5*(bilevelmodel->leader_ub-Fc);
        gen_max_probab = gen_ref_pt_probab + input.getScalingGeneral()*0.5*(bilevelmodel->leader_ub-Fc);
    }if(input.getTypeDepGeneral() == Input::TypesDepGeneral::StrongFragile){
        gen_ref_pt_probab = 1.0;
        gen_min_probab = std::exp(-input.getScalingGeneral()*0.5*(bilevelmodel->leader_ub-Fc)/(bilevelmodel->leader_ub-bilevelmodel->leader_lb));
        gen_max_probab = std::exp(input.getScalingGeneral()*0.5*(bilevelmodel->leader_ub-Fc)/(bilevelmodel->leader_ub-bilevelmodel->leader_lb));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::generateScenarios(){
    // Define SAA scenarios for Decision-Dependent Strong Weak using SAA method
    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepStrongWeak &&
        input.getSolverApproach() == Input::SolverApproach::TR_SAA){
        scenarios_dep_strwk.resize(nb_saa_problems);
        for(int pr = 0; pr < nb_saa_problems; ++pr){
            for(int s = 0; s < nb_saa_scenarios; ++s){
                scenarios_dep_strwk[pr].push_back(inverseDepStrongWeakScenario());
            }
            std::sort(scenarios_dep_strwk[pr].begin(),scenarios_dep_strwk[pr].end(),std::greater<double>());
        }
    }

    // Compute null space for equality constraints for formulation Decision-Dependent General
    Eigen::MatrixXd N;
    if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false) 
        nullSpaceFollowerEqConstrs(N);

    // Define SAA scenarios for Decision-Dependent General 
    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral){
        dir_scenarios_dep_gen.resize(nb_saa_problems);
        step_scenarios_dep_gen.resize(nb_saa_problems);
        accep_scenarios_dep_gen.resize(nb_saa_problems);
        for(int pr = 0; pr < nb_saa_problems; ++pr){
            dir_scenarios_dep_gen[pr].resize(nb_saa_scenarios);
            for(int s = 0; s < nb_saa_scenarios; ++s){
                dirDepGeneralScenario(N,dir_scenarios_dep_gen[pr][s]);
                step_scenarios_dep_gen[pr].push_back(uniform(rd_engine));
                accep_scenarios_dep_gen[pr].push_back(uniform(rd_engine));
            }
            defineGeneralStepInfo(pr,coeff_obj_alpha[pr],coeff_alpha[pr],constrs_alpha[pr],bdlbs_alpha[pr],bdubs_alpha[pr]);
        }
    }

    // Define evaluation scenarios for the Decision-Dependent General case.
    // Strong Weak cases have only 2 endogenous scenarios, that can be enumerated for evaluation
    dir_eval_dep_gen.resize(nb_val_scenarios);
    for(int s = 0; s < nb_val_scenarios; ++s){
        dirDepGeneralScenario(N,dir_eval_dep_gen[s]);
        step_eval_dep_gen.push_back(uniform(rd_engine));
        accep_eval_dep_gen.push_back(uniform(rd_engine));
    }
    defineGeneralStepInfo(nb_saa_problems,eval_coeff_obj_alpha,eval_coeff_alpha,eval_constrs_alpha,eval_bdlbs_alpha,eval_bdubs_alpha);
}

double Instance::inverseDepStrongWeakScenario(){
    if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional)
        return bilevelmodel->follower_ub - uniform_eps(rd_engine)*(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Threshold)
        return (bilevelmodel->follower_ub + bilevelmodel->follower_lb)/2.0 + (1.0/input.getParamThreshold())*std::log((1.0/uniform_eps(rd_engine) - 1.0));
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Strong)
        return bilevelmodel->follower_lb + (bilevelmodel->follower_ub - bilevelmodel->follower_lb)*(1.0 + (1.0/input.getParamStrong())*std::log((1.0 - uniform_eps(rd_engine)*(1.0 - std::exp(-input.getParamStrong())))));
    else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Fragile)
        return bilevelmodel->follower_lb - ((bilevelmodel->follower_ub - bilevelmodel->follower_lb)/input.getParamFragile())*std::log((std::exp(-input.getParamFragile()) + uniform_eps(rd_engine)*(1.0 - std::exp(-input.getParamFragile()))));
    else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior."));
}

void Instance::nullSpaceFollowerEqConstrs(Eigen::MatrixXd & N){
    // Number of equality constraints follower problem + optimality constraint when follower is optimal
    int nb_total_eq = bilevelmodel->nb_follower_eq_constrs;
    if(input.isFollowerNearOpt() == false) nb_total_eq += 1;

    // Defining Eigen Matrix of EqConstrs
    Eigen::MatrixXd EqConstrs(nb_total_eq, bilevelmodel->nb_follower_vars);
    int constr = 0;
    for(int i = 0; i < bilevelmodel->nb_follower_constrs; ++i){
        if(bilevelmodel->follower_constrs[i].sense == '='){
            for(int j = 0; j < bilevelmodel->nb_follower_vars; ++j){
                BilevelVariable var = bilevelmodel->follower_vars[j];
                EqConstrs(constr,j) = bilevelmodel->follower_constrs[i].coeffs[var.id];
            }
            ++constr;
        }
    }
    if(input.isFollowerNearOpt() == false){
        for(int id = 0; id < bilevelmodel->nb_follower_vars; ++id){
            EqConstrs(constr,id) = bilevelmodel->follower_vars[id].obj_follower;
        }
    }

    // Compute null space of EqConstrs
    Eigen::FullPivLU<Eigen::MatrixXd> lu(EqConstrs);
    N = lu.kernel();
}

void Instance::dirDepGeneralScenario(Eigen::MatrixXd & N, std::vector<double> & dir){
    if(coordinate_mcmc == true && bilevelmodel->nb_follower_eq_constrs == 0){
        double chosen = uniform_int(rd_engine);
        double direct = uniform(rd_engine);

        dir.resize(bilevelmodel->nb_follower_vars, 0.0);
        if(direct <= 0.5) dir[chosen] = 1.0;
        else dir[chosen] = -1.0;

        dir[chosen] = 0.0;
        dir[3] = 0.0765625;
        dir[8] = -0.107188;
        dir[9] = -0.06125;
    }else{
        if(bilevelmodel->nb_follower_eq_constrs > 0 || input.isFollowerNearOpt() == false){
            if(N.cols() == 0) {
                dir.resize(bilevelmodel->nb_follower_vars,0.0); 
            }else{
                // Sample a random vector in the null space
                Eigen::VectorXd aux = N * (Eigen::VectorXd::NullaryExpr(N.cols(), [&](){ return normal(rd_engine); })); 
                // Normalize to unit length
                aux.normalize();
                // Pass information to vector 
                dir.assign(aux.data(), aux.data() + aux.size());
            }
        }else{
            // Sample a normal random vector 
            Eigen::VectorXd aux(bilevelmodel->nb_follower_vars);
            for (int id = 0; id < bilevelmodel->nb_follower_vars; ++id) aux[id] = normal(rd_engine);
            // Normalize to unit length
            aux.normalize();
            // Pass information to vector 
            dir.assign(aux.data(), aux.data() + aux.size());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::getStrongWeakProbab(double opt_val_follower) const{
    if(input.getFollowerBehavior() == Input::FollowerBehavior::FixStrongWeak) 
        return input.getFixCoopLevel();
    else if(input.getFollowerBehavior() == Input::FollowerBehavior::DepStrongWeak){
        if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional)
            return (bilevelmodel->follower_ub - opt_val_follower)/(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Threshold){
            double phi_c = (bilevelmodel->follower_lb + bilevelmodel->follower_ub) / 2.0;
            return 1.0 / (1.0 + std::exp(input.getParamThreshold() * (opt_val_follower - phi_c)));
        }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Strong){
            double t = (opt_val_follower - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
            return (1.0 - std::exp(-input.getParamStrong() * (1.0 - t))) / (1.0 - std::exp(-input.getParamStrong()));
        }else if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Fragile){
            double t = (opt_val_follower - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
            return 1.0 - (1.0 - std::exp(-input.getParamFragile() * t)) / (1.0 - std::exp(-input.getParamFragile()));
        }else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior."));
    }
    else throw std::runtime_error(std::string("Incorrect choice of Behavior."));
}

double Instance::getEvalDepStrongWeakProbab(double opt_val_follower, Input::TypesDepStrongWeak eval_behavior, double param) const {
    if(eval_behavior == Input::TypesDepStrongWeak::Proportional)
        return (bilevelmodel->follower_ub - opt_val_follower)/(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    else if(eval_behavior== Input::TypesDepStrongWeak::Threshold){
        double phi_c = (bilevelmodel->follower_lb + bilevelmodel->follower_ub) / 2.0;
        return 1.0 / (1.0 + std::exp(param * (opt_val_follower - phi_c)));
    }else if(eval_behavior == Input::TypesDepStrongWeak::Strong){
        double t = (opt_val_follower - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return (1.0 - std::exp(-param * (1.0 - t))) / (1.0 - std::exp(-param));
    }else if(eval_behavior == Input::TypesDepStrongWeak::Fragile){
        double t = (opt_val_follower - bilevelmodel->follower_lb) / (bilevelmodel->follower_ub - bilevelmodel->follower_lb);
        return 1.0 - (1.0 - std::exp(-param * t)) / (1.0 - std::exp(-param));
    }else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent Strong Weak Behavior."));
}

////////////////////////////////////////////////////////////////////////////////////////////////

double Instance::getDepGeneralProbab(double opt_val_follower, const double * x_, std::vector<double> & y_) const{
    // Compute leader value
    double leader_value = 0.0;
    for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
        leader_value += bilevelmodel->leader_vars[i].obj_leader*x_[i];
    }
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        leader_value += bilevelmodel->follower_vars[i].obj_leader*y_[i];
    }

    // Central point leader value
    double Fc = (bilevelmodel->leader_lb + bilevelmodel->leader_ub)/2.0;

    // Select interval coefficient
    double follower_value = (bilevelmodel->follower_ub - opt_val_follower)/(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    double m = 0.0;
    for(int k = 0; k < gen_nb_intervals; ++k){
        if(follower_value >= gen_intervals[k] && follower_value <= gen_intervals[k+1]){
            m = input.getScalingGeneral()*gen_coeff_intervals[k];
            break;
        }
    }

    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional) return (gen_ref_pt_probab + m*(leader_value - Fc));
    else if(input.getTypeDepGeneral() == Input::TypesDepGeneral::StrongFragile) return (gen_ref_pt_probab*std::exp(m*(leader_value-Fc)));
    else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent General Behavior."));
}

double Instance::getEvalDepGeneralProb(double opt_val_follower, const double * x_, std::vector<double> & y_, Input::TypesDepGeneral eval_behavior, int nb_intervals, double scaling) const {
    // Compute leader value
    double leader_value = 0.0;
    for(int i = 0; i < bilevelmodel->nb_leader_vars; ++i){
        leader_value += bilevelmodel->leader_vars[i].obj_leader*x_[i];
    }
    double only_leader_value = 0.0;
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        leader_value += bilevelmodel->follower_vars[i].obj_leader*y_[i];
        only_leader_value += bilevelmodel->follower_vars[i].obj_leader*y_[i];
    }
    std::cout << "leader value: " << leader_value << " " << only_leader_value << std::endl;

    // Central point leader value
    double Fc = (bilevelmodel->leader_lb + bilevelmodel->leader_ub)/2.0;

    // Select interval coefficient
    double follower_value = (bilevelmodel->follower_ub - opt_val_follower)/(bilevelmodel->follower_ub - bilevelmodel->follower_lb);
    double interval = 1.0/(double)nb_intervals;
    double m = 0.0;
    for(int i = 0; i < nb_intervals; ++i){
        if(i*interval <= follower_value && follower_value <= (i+1)*interval){
            m = scaling*(0.5 - ((double)i*interval + (double)(i+1.0)*interval)/2.0);
            break;
        }
    }

    double ref_probab = 1.0; 
    if(eval_behavior == Input::TypesDepGeneral::GenProportional){
        double delta = (bilevelmodel->leader_ub - bilevelmodel->leader_lb)/2.0;
        ref_probab = scaling*0.5*delta;
        std::cout << "proportional " << (ref_probab + m*(leader_value - Fc)) << std::endl;
    }

    if(eval_behavior == Input::TypesDepGeneral::StrongFragile){
        std::cout << "strongfragile " << (ref_probab*std::exp(m*(leader_value-Fc)/(bilevelmodel->leader_ub - bilevelmodel->leader_lb) )) << std::endl;
    }

    if(eval_behavior == Input::TypesDepGeneral::GenProportional) return (ref_probab + m*(leader_value - Fc));
    else if(eval_behavior == Input::TypesDepGeneral::StrongFragile) return (ref_probab*std::exp(m*(leader_value-Fc)/(bilevelmodel->leader_ub - bilevelmodel->leader_lb) ));
    else throw std::runtime_error(std::string("Incorrect choice of Decision-Dependent General Behavior."));
}

////////////////////////////////////////////////////////////////////////////////////////////////

void Instance::defineGeneralStepInfo(int pr, std::vector<double> & coeff_obj,
                                            std::vector<std::vector<double>> & coeff,
                                            std::vector<std::vector<std::vector<int>>> & constrs, 
                                            std::vector<std::vector<std::vector<int>>> & bdlbs, 
                                            std::vector<std::vector<std::vector<int>>> & bdubs) {
    int num_scenarios = nb_saa_scenarios;
    if(pr == nb_saa_problems) num_scenarios = nb_val_scenarios;
    
    // Coefficient step variables (alpha_min and alpha_max).
    coeff.resize(num_scenarios);
    coeff_obj.resize(num_scenarios);
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
                if(pr == nb_saa_problems) coeff[s][c] += constr.coeffs[var.id]*getGeneralDirectionEval(s,i);
                else coeff[s][c] += constr.coeffs[var.id]*getGeneralDirection(pr,s,i);
            }
            if(constr.sense == '<'){
                if(coeff[s][c] >= 1e-6) constrs[s][1].push_back(c);
                else if(coeff[s][c] <= -1e-6) constrs[s][0].push_back(c);
            }
            if(constr.sense == '>'){
                if(coeff[s][c] >= 1e-6) constrs[s][0].push_back(c);
                else if(coeff[s][c] <= -1e-6) constrs[s][1].push_back(c);
            }
        }
    }

    // Define coefficients considering follower objective.
    for(int s = 0; s < num_scenarios; ++s){
        coeff_obj[s] = 0.0;
        for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
            BilevelVariable var = bilevelmodel->follower_vars[i];
            
            if(pr == nb_saa_problems) coeff_obj[s] += var.obj_follower*getGeneralDirectionEval(s,i);
            else coeff_obj[s] += var.obj_follower*getGeneralDirection(pr,s,i);
        }
        // Also consider near optimality constraint in this case
        if(input.isFollowerNearOpt() == true){ 
            // This constraints is <=.
            if(coeff_obj[s] >= 1e-6) constrs[s][1].push_back(bilevelmodel->nb_follower_constrs);
            else if(coeff_obj[s] <= -1e-6) constrs[s][0].push_back(bilevelmodel->nb_follower_constrs);
        }
    }

    // Define bound constraints for step problem (to define alpha variables).
    for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        for(int s = 0; s < num_scenarios; ++s){
            double dir = 0.0;
            if(pr == nb_saa_problems) dir = getGeneralDirectionEval(s,i);
            else dir = getGeneralDirection(pr,s,i);

            if(!(var.lb <= -bilevelmodel->inf)){
                if(dir >= 1e-6) bdlbs[s][0].push_back(i);
                else if(dir <= -1e-6) bdlbs[s][1].push_back(i);
            }
            if(!(var.ub >=  bilevelmodel->inf)){
                if(dir >= 1e-6) bdubs[s][1].push_back(i);
                else if(dir <= -1e-6) bdubs[s][0].push_back(i);
            }
        }
    }
}

double Instance::defineEvalMinStep(int s, const double * x_, std::vector<double> & y_, double fs_) const {
    // Lower bound on step (alpha_min).
    double alpha_min = -getStepSizeInterval();

    // Constraints defining alpha_min.
    for(int c : eval_constrs_alpha[s][0]){
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

            if(constr.sense == '>' && (constr.rhs - value) >= 1e-6){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }
            if(constr.sense == '<' && (constr.rhs - value) <= -1e-6){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }

            // Avoid problems with numeric instability.
            double diff = (constr.rhs - value);
            if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha value.
            alpha_min = std::max(alpha_min, diff/eval_coeff_alpha[s][c]);

            std::cout << "min \t" << constr.rhs << " " << value << " " << eval_coeff_alpha[s][c] << " " << diff/eval_coeff_alpha[s][c] << std::endl;
        }
        // Near optimality constraint.
        else{
            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                BilevelVariable var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-6){
                throw std::runtime_error("Constraint on near optimality was not respected by current y.");
            }
            
            // Avoid problems with numerical instability.
            double diff = ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value);
            if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha min.
            alpha_min = std::max(alpha_min, diff/eval_coeff_obj_alpha[s]);

            std::cout << "min obj\t" << (1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub << " " << value << " " << eval_coeff_obj_alpha[s] << " " << ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value)/eval_coeff_obj_alpha[s] << std::endl;
        }
    }

    // Lower bound constraints.
    for(int i: eval_bdlbs_alpha[s][0]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.lb - y_[i]) >= 1e-6){
            throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y.");
        }
        
        // Update alpha min.
        alpha_min = std::max(alpha_min, (var.lb - y_[i])/getGeneralDirectionEval(s,i));
        
        std::cout << "min lb\t" << var.lb << " " << y_[i] << " " << getGeneralDirectionEval(s,i) << " " << (var.lb - y_[i])/getGeneralDirectionEval(s,i) << std::endl;
    }

    // Upper bound constraints.
    for(int i: eval_bdubs_alpha[s][0]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.ub - y_[i]) <= -1e-6){
            throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Update alpha min.
        alpha_min = std::max(alpha_min, (var.ub - y_[i])/getGeneralDirectionEval(s,i));
        
        std::cout << "min ub\t" << var.ub << " " << y_[i] << " " << getGeneralDirectionEval(s,i) << " " << (var.ub - y_[i])/getGeneralDirectionEval(s,i) << std::endl;
    }

    if(alpha_min >= 1e-6) throw std::runtime_error(std::string("Incorrect min step value: Infeasible current point."));

    return alpha_min;
}

double Instance::defineEvalMaxStep(int s, const double * x_, std::vector<double> & y_, double fs_) const {
    // Upper bound on step (alpha_max).
    double alpha_max = getStepSizeInterval();

    // Constraints defining alpha_max.
    for(int c : eval_constrs_alpha[s][1]){
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

            if(constr.sense == '>' && (constr.rhs - value) >= 1e-6){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }
            if(constr.sense == '<' && (constr.rhs - value) <= -1e-6){
                throw std::runtime_error("Constraint " + std::to_string(c) + " was not respected by current y.");
            }

            // Avoid problems with numeric instability.
            double diff = (constr.rhs - value);
            if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha max.
            alpha_max = std::min(alpha_max, diff/eval_coeff_alpha[s][c]);

            std::cout << "max \t" << constr.rhs << " " << value << " " << eval_coeff_alpha[s][c] << " " << diff/eval_coeff_alpha[s][c] << std::endl;
        }
        // Near optimality constraint.
        else{
            double value = 0.0;
            for(int i = 0; i < bilevelmodel->nb_follower_vars; ++i){
                BilevelVariable var = bilevelmodel->follower_vars[i];
                value += var.obj_follower*y_[i];
            }

            if(((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value) <= -1e-6){
                throw std::runtime_error("Constraint on near optimality was not respected by current y.");
            }

            // Avoid problems with numerical instability.
            double diff = ((1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub - value);
            if(diff >= -1e-6 && diff <= 1e-6) diff = 0.0;

            // Update alpha_max value.
            alpha_max = std::min(alpha_max, diff/eval_coeff_obj_alpha[s]);

            std::cout << "max obj\t" << (1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*bilevelmodel->follower_ub << " " << value << " " << eval_coeff_obj_alpha[s] << " " << diff/eval_coeff_obj_alpha[s] << std::endl;
        }
    }

    // Lower bound constraints.
    for(int i: eval_bdlbs_alpha[s][1]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.lb - y_[i]) >= 1e-6){
            throw std::runtime_error("Lower bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Update alpha max.
        alpha_max = std::min(alpha_max, (var.lb - y_[i])/getGeneralDirectionEval(s,i));

        std::cout << "max lb\t" << var.lb << " " << y_[i] << " " << getGeneralDirectionEval(s,i) << " " << (var.lb - y_[i])/getGeneralDirectionEval(s,i) << std::endl;
    }

    // Upper bound constraints.
    for(int i: eval_bdubs_alpha[s][1]){
        BilevelVariable var = bilevelmodel->follower_vars[i];

        if((var.ub - y_[i]) <= -1e-6){
            throw std::runtime_error("Upper bound constraint " + std::to_string(i) + " was not respected by current y.");
        }

        // Define alpha max.
        alpha_max = std::min(alpha_max, (var.ub - y_[i])/getGeneralDirectionEval(s,i));

        std::cout << "max ub\t" << var.ub << " " << y_[i] << " " << getGeneralDirectionEval(s,i) << " " << (var.ub - y_[i])/getGeneralDirectionEval(s,i) << std::endl;
    }

    if(alpha_max <= -1e-6) throw std::runtime_error(std::string("Incorrect max step value: Infeasible current point."));

    return alpha_max;
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
