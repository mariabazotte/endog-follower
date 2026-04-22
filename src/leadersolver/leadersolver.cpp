#include "leadersolver.hpp"
#include "../followersolver/followersolver.hpp" 
#include "leadersolverfactory.hpp"

double inf_norm_diff(const double* a, const double* b, int n){
    double norm = 0.0;
    for (int i = 0; i < n; ++i)
        norm = std::max(norm, std::abs(a[i] - b[i]));
    return norm;
}

std::string AbstractLeaderSolver::writeHeadComp() const {
    std::string result = "FOLLOWER_BEHAVIOR;COOPERATION_TYPE;PARAMS;";
    result += "LEADER_OBJ_MEAN;LEADER_OBJ_VARIANCE;LEADER_VARS_ESS;LEADER_VARS_RHAT;";
    result += "FOLLOWER_OBJ_OPT;";
    if(input.isFollowerNearOpt() == true){
        result += "FOLLOWER_OBJ_NEAROPT_MEAN;FOLLOWER_OBJ_NEAROPT_VARIANCE;";
    }
    return result;
}

void AbstractLeaderSolver::writeCompFixStrongWeak(std::string & output, double cooperation, 
                                                double eval, double fs_, double f_eval) const {
    output += std::to_string(Input::FollowerBehavior::FixStrongWeak) + ";-;" + 
              std::to_string(cooperation) + ";" + 
              std::to_string(eval) + ";-;-;-;" + std::to_string(fs_) + ";";
    if(input.isFollowerNearOpt() == true) 
        output += std::to_string(f_eval) + ";-;";
    else output += "-;-;";
    output += "\n";
}

void AbstractLeaderSolver::writeCompDepStrongWeak(std::string & output, Input::TypesDepStrongWeak type, double param, 
                                                double eval, double fs_, double f_eval) const {
    output += std::to_string(Input::FollowerBehavior::DepStrongWeak) + ";" + 
              std::to_string(type) + ";" + std::to_string(param) + ";" + 
              std::to_string(eval) + ";-;-;-;" + std::to_string(fs_) + ";";
    if(input.isFollowerNearOpt() == true) 
        output += std::to_string(f_eval) + ";-;";
    else output += "-;-;";
    output += "\n";
}

void AbstractLeaderSolver::writeCompDepGeneral(std::string & output, Input::TypesDepGeneral type, int nb_int, double param, 
                                                double eval, double var_eval, double ess, double rhat,
                                                double fs_, double f_eval, double f_var_eval) const {
    output += std::to_string(Input::FollowerBehavior::DepGeneral) + ";" + 
              std::to_string(type) + ";" + std::to_string(nb_int) + "/" + std::to_string(param) + ";" +
              std::to_string(eval) + ";" + std::to_string(var_eval) + ";" +
              std::to_string(ess) + ";" + std::to_string(rhat) + ";" + std::to_string(fs_) + ";";
    if(input.isFollowerNearOpt() == true){  
        output += std::to_string(f_eval) + ";" + std::to_string(f_var_eval) + ";";
    } 
    output += "\n";
}

std::string AbstractLeaderSolver::writeComp() const {
    // Define output string.
    std::string output = "";

    if(nb_sol > 0){
        // Define leader objective involving only leader variables.
        double x_leader_obj = 0.0;
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
            x_leader_obj += instance.getModel()->leader_vars[i].obj_leader*getX_(i);
        
        // Auxiliar parameters.
        double Fs_ = getFollower()->getLeaderOptObj();
        double Fw_ = getFollower()->getLeaderPesObj();
        double fs_ = getFollower()->getFollowerOptimalObj();
        double fs_eps_ = (input.isFollowerNearOpt() == true) ? getFollower()->getFollowerNearOptimalOptObj() : 0.0;
        double fw_eps_ = (input.isFollowerNearOpt() == true) ? getFollower()->getFollowerNearOptimalPesObj() : 0.0;
        double eval, f_eval, var_eval, f_var_eval, ess, rhat;

        // Fixed Strong-Weak configurations.
        for(double beta: instance.getFixStrongWeakCoopConfigs()){
            instance.evaluateFixStrongWeak(eval,f_eval,Fs_,Fw_,fs_eps_,fw_eps_,beta);
            writeCompFixStrongWeak(output,beta,(x_leader_obj+eval),fs_,f_eval);
        }

        // Dependent Strong-Weak configurations.
        
        // Proportional.
        instance.evaluateDepStrongWeak(eval,f_eval,Fs_,Fw_,fs_,fs_eps_,fw_eps_,Input::TypesDepStrongWeak::Proportional,0.0);
        writeCompDepStrongWeak(output,Input::TypesDepStrongWeak::Proportional,0.0,(x_leader_obj+eval),fs_,f_eval);
        
        // Threshold.
        for(double param : instance.getDepStrongWeakThrConfigs()){
            instance.evaluateDepStrongWeak(eval,f_eval,Fs_,Fw_,fs_,fs_eps_,fw_eps_,Input::TypesDepStrongWeak::Threshold,param);
            writeCompDepStrongWeak(output,Input::TypesDepStrongWeak::Threshold,param,(x_leader_obj+eval),fs_,f_eval);
        }

        // Strong.
        for(double param : instance.getDepStrongWeakStrConfigs()){
            instance.evaluateDepStrongWeak(eval,f_eval,Fs_,Fw_,fs_,fs_eps_,fw_eps_,Input::TypesDepStrongWeak::Strong,param);
            writeCompDepStrongWeak(output,Input::TypesDepStrongWeak::Strong,param,(x_leader_obj+eval),fs_,f_eval);
        }

        // Fragile.
        for(double param : instance.getDepStrongWeakFrgConfigs()){
            instance.evaluateDepStrongWeak(eval,f_eval,Fs_,Fw_,fs_,fs_eps_,fw_eps_,Input::TypesDepStrongWeak::Fragile,param);
            writeCompDepStrongWeak(output,Input::TypesDepStrongWeak::Fragile,param,(x_leader_obj+eval),fs_,f_eval);
        }

        // Strong Power.
        for(double param : instance.getDepStrongWeakStrPowConfigs()){
            instance.evaluateDepStrongWeak(eval,f_eval,Fs_,Fw_,fs_,fs_eps_,fw_eps_,Input::TypesDepStrongWeak::StrongPower,param);
            writeCompDepStrongWeak(output,Input::TypesDepStrongWeak::StrongPower,param,(x_leader_obj+eval),fs_,f_eval);
        }

        // Fragile Power.
        for(double param : instance.getDepStrongWeakFrgPowConfigs()){
            instance.evaluateDepStrongWeak(eval,f_eval,Fs_,Fw_,fs_,fs_eps_,fw_eps_,Input::TypesDepStrongWeak::FragilePower,param);
            writeCompDepStrongWeak(output,Input::TypesDepStrongWeak::FragilePower,param,(x_leader_obj+eval),fs_,f_eval);
        }

        // Dependent General configurations.
        const double * x_ = getX_();
        const double * yi_ = getFollower()->getYi_();

        // Uniform.
        if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == true)){ eval = Fs_; var_eval = 0.0; f_eval = fs_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
        else if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == false) && (std::abs(fw_eps_ - fs_eps_) <= 1e-3)){ eval = Fs_; var_eval = 0.0; f_eval = fs_eps_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
        else instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::Neutral,0,0.0);
        writeCompDepGeneral(output,Input::TypesDepGeneral::Neutral,0,0.0,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);

        // Proportional.
        for(int nb_int: instance.getDepGeneralIntConfigs()){
            if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == true)){ eval = Fs_; var_eval = 0.0; f_eval = fs_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
            else if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == false) && (std::abs(fw_eps_ - fs_eps_) <= 1e-3)){ eval = Fs_; var_eval = 0.0; f_eval = fs_eps_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
            else instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenProportional,nb_int,1.0);
            writeCompDepGeneral(output,Input::TypesDepGeneral::GenProportional,nb_int,1.0,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
        }

        // Fragile. 
        for(int nb_int: instance.getDepGeneralIntConfigs()){
            for(double scal: instance.getDepGeneralScalFrgConfigs()){
                if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == true)){ eval = Fs_; var_eval = 0.0; f_eval = fs_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
                else if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == false) && (std::abs(fw_eps_ - fs_eps_) <= 1e-3)){ eval = Fs_; var_eval = 0.0; f_eval = fs_eps_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
                else instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenFragile,nb_int,scal);
                writeCompDepGeneral(output,Input::TypesDepGeneral::GenFragile,nb_int,scal,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
            }
        }

        // Fragile Power. 
        for(int nb_int: instance.getDepGeneralIntConfigs()){
            for(double scal: instance.getDepGeneralScalFrgPowConfigs()){
                if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == true)){ eval = Fs_; var_eval = 0.0; f_eval = fs_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
                else if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == false) && (std::abs(fw_eps_ - fs_eps_) <= 1e-3)){ eval = Fs_; var_eval = 0.0; f_eval = fs_eps_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
                else instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenFragilePower,nb_int,scal);
                writeCompDepGeneral(output,Input::TypesDepGeneral::GenFragilePower,nb_int,scal,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
            }
        }

        // Strong Power. 
        for(int nb_int: instance.getDepGeneralIntConfigs()){
            for(double scal: instance.getDepGeneralScalStrPowConfigs()){
                if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == true)){ eval = Fs_; var_eval = 0.0; f_eval = fs_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
                else if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == false) && (std::abs(fw_eps_ - fs_eps_) <= 1e-3)){ eval = Fs_; var_eval = 0.0; f_eval = fs_eps_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
                else instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenStrongPower,nb_int,scal);
                writeCompDepGeneral(output,Input::TypesDepGeneral::GenStrongPower,nb_int,scal,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
            }
        }
    }
    return output;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

LeaderSolver::LeaderSolver(const Input & input, Instance & instance, std::string name) : 
                                            LeaderSolver(input,instance,name,-1,new GRBEnv(),instance.getModel()->leader_lb,instance.getModel()->leader_ub,0.0){} 

LeaderSolver::LeaderSolver(const Input & input, Instance & instance, std::string name, int pr, GRBEnv *env, double lb, double ub, double tp) : 
                                AbstractLeaderSolver(input,instance,name,env), pr(pr), model(new GRBModel(*env)) {
    general_lb = lb;
    general_ub = ub;
    time_prep = tp;

    // Create leader problem
    params();
    create();
    
    // Create follower problem
    FollowerSolverFactory factory;
    follower = factory.createFollowerSolver(input,instance,pr,this);
    x_ = new double[instance.getModel()->nb_leader_vars];
}

LeaderSolver::~LeaderSolver(){
    delete follower;
    if(x_) delete[] x_;
    if(x) delete[] x;
    if(x_binary) delete[] x_binary;
    delete model;
    if(pr == -1) delete env;
}

void LeaderSolver::params(){
    model->set(GRB_IntParam_Threads, input.getNbThreads());
    model->set(GRB_DoubleParam_TimeLimit, std::max(0.0,input.getTimeLimit() - time_prep));
    model->set(GRB_DoubleParam_MIPGap, 1e-4);
    model->set(GRB_DoubleParam_IntFeasTol, 1e-5);
    model->set(GRB_DoubleParam_FeasibilityTol, 1e-6);
    // model->set(GRB_IntParam_Presolve, 0);
    if(input.getVerbose() < 1) model->set(GRB_IntParam_OutputFlag, 0);
    else model->set(GRB_IntParam_OutputFlag, 1);
}

void LeaderSolver::create(){
    x = model->addVars(instance.getModel()->nb_leader_vars,GRB_CONTINUOUS);
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
        BilevelVariable var = instance.getModel()->leader_vars[i];

        if(var.lb <= -instance.getModel()->inf) 
            x[i].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
        else x[i].set(GRB_DoubleAttr_LB,var.lb);
        if(var.ub >=  instance.getModel()->inf) 
            x[i].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        else x[i].set(GRB_DoubleAttr_UB,var.ub);

        x[i].set(GRB_StringAttr_VarName,var.name);
        x[i].set(GRB_DoubleAttr_Obj,var.obj_leader);
        x[i].set(GRB_CharAttr_VType,var.type);
    }

    // Non-coupling constraints
    for(int c = 0; c < instance.getModel()->nb_leader_constrs; ++c){
        BilevelConstraint constr = instance.getModel()->leader_constrs[c];
        
        GRBLinExpr expr = 0;
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
            BilevelVariable var = instance.getModel()->leader_vars[i];
            expr += constr.coeffs[var.id]*x[i];
        }

        if(constr.sense == '=') model->addConstr(expr == constr.rhs, constr.name);
        else if(constr.sense == '<') model->addConstr(expr <= constr.rhs, constr.name);
        else if(constr.sense == '>') model->addConstr(expr >= constr.rhs, constr.name);
    }
}

void LeaderSolver::write_inf(){
    model->computeIIS();
    model->write(name + ".ilp");
    model->write(name + ".lp");
}

void LeaderSolver::write_sol(){
    model->write(name + ".lp");
    model->write(name + ".sol");
}

void LeaderSolver::verify_stat(){
    if(status == Status::Unbounded){
        std::cerr << "The " << name << " is unbounded, this cannot happen according to the problem definition." << std::endl;
        throw std::runtime_error( name + " is unbounded.\n");
        exit(0);
    }
    if(status == Status::Infeasible_or_Unbounded){
        std::cerr << "The " << name << " is unbounded or infeasible, this cannot happen according to the problem definition." << std::endl;
        write_inf();
        throw std::runtime_error(name + " is unbounded or infeasible.\n");
        exit(0); 
    } 
    if(status == Status::Infeasible){
        std::cerr << "The " << name << " is infeasible, this cannot happen according to the problem definition." << std::endl;
        write_inf();
        throw std::runtime_error(name + " is infeasible.\n");
        exit(0); 
    }
}

void LeaderSolver::print(){
    std::string output;
    if(status == Status::Optimal){
        output = "OPTIMAL: The optimal solution of the " + name + " was " + std::to_string(ub) + ".\n";
    }else if(status == Status::Time_Limit){
        output = "TIME: The best bound found for " + name  + " was " + std::to_string(lb) + ".\n";
        output = "TIME: The best feasible solution found for " + name + " was " + std::to_string(ub) + ".\n";
    }
    else if(status == Status::Iteration_Limit){
        output = "ITERATION: The best bound found for " + name  + " was " + std::to_string(lb) + ".\n";
        output = "ITERATION: The best feasible solution found for " + name + " was " + std::to_string(ub) + ".\n";
    }
    std::cout << output; 
    if(input.getVerbose() == 2) write_sol(); 
}

void LeaderSolver::upd_solution(){
    lb = model->get(GRB_DoubleAttr_ObjBound);
    ub = model->get(GRB_DoubleAttr_ObjVal);
    nb_bnb = model->get(GRB_DoubleAttr_NodeCount);
    nb_sol = model->get(GRB_IntAttr_SolCount);
    gap = (ub-lb)/std::abs(ub);
    if(nb_sol > 0) get_variables();
}

void LeaderSolver::get_variables(){
    if(x_) delete[] x_;
    x_ = model->get(GRB_DoubleAttr_X,x,instance.getModel()->nb_leader_vars);
}

void LeaderSolver::get_variables(double & inf_norm){
    // Copy previous solution for comparison.
    double * prev_x_ = new double[instance.getModel()->nb_leader_vars];
    std::copy(x_, x_ + instance.getModel()->nb_leader_vars, prev_x_);

    // Get current solution.
    nb_sol = model->get(GRB_IntAttr_SolCount);
    if(nb_sol > 0) get_variables();

    // Compute infinity norm of the difference between current and previous solutions.
    inf_norm = inf_norm_diff(prev_x_, x_, instance.getModel()->nb_leader_vars);

    // Delete previous solution.
    delete[] prev_x_;
}

void LeaderSolver::solve(){
    follower->solve();
}

std::string LeaderSolver::write() const {
    double m = 0.0; 
    double fm = 0.0;
    double v, fv;
    if(nb_sol > 0) getFollower()->evaluate(m, v, fm, fv); // m is the expected value of leader objective, fm is the expected value of follower objective. Variance here should be zero.

    std::string output = std::to_string(status) + std::string(";") + 
                        std::to_string(lb) + std::string(";") + 
                        std::to_string(ub) + std::string(";") + 
                        std::to_string(gap) + std::string(";") + 
                        std::to_string(m) + std::string(";") +
                        std::to_string(follower->getFollowerOptimalObj()) + std::string(";");
    
    if(input.isFollowerNearOpt() == true){
        output += std::to_string(follower->getFollowerNearOptimalOptObj()) + std::string(";") +
                  std::to_string(follower->getFollowerNearOptimalPesObj()) + std::string(";") +
                  std::to_string(fm) + std::string(";");
    }
    
    output += std::to_string(time_) + std::string(";") +
              std::to_string(nb_bnb) + std::string(";") +
              std::to_string((nb_sol > 0)) + std::string(";");
    
    return output;
}

std::string LeaderSolver::writeHead() const {
    std::string output = "STATUS;LB;UB;GAP;EVAL;OBJ_FOLLOWER;";
    if(input.isFollowerNearOpt() == true)
        output += "OBJ_OPT_FOLLOWER;OBJ_PES_FOLLOWER;OBJ_MEAN_FOLLOWER;";
    output += "TIME;NB_BNB;NB_SOL;";
    return output;
}

void LeaderSolver::defineBinaryVariables() {
    total_nb_x_binary = 0;
    nb_x_binary.resize(instance.getModel()->nb_leader_vars,0);
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
        BilevelVariable var = instance.getModel()->leader_vars[i];

        if(var.type == 'C') 
            throw std::runtime_error("Callback method only defined for integer or binary leader variables. Continuous leader variables are not accepted.");

        int r = var.ub - var.lb + 1;   
        int k = 0; int p = 1;
        while(p < r) { p <<= 1; ++k; }

        total_nb_x_binary += k;
        nb_x_binary[i] = k;
    }

    int j = 0;
    x_binary = model->addVars(total_nb_x_binary,GRB_BINARY);
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
        const BilevelVariable & var = instance.getModel()->leader_vars[i];

        GRBLinExpr expr = var.lb;
        for(int k = 0; k < nb_x_binary[i]; ++k){
            expr += (1LL << k) * x_binary[j+k];
        }
        model->addConstr(x[i] == expr);
        j+= nb_x_binary[i];
    }
}

///////////////////////////////////////////////////////////////////////////////////

SAALeaderSolver::SAALeaderSolver(const Input & input, Instance & instance, std::string name) : 
                                    AbstractLeaderSolver(input,instance,name,new GRBEnv()) {
    var_lb = 0.0;
    var_ub = 0.0;
    var_gap = 0.0;
    var_nb_bnb = 0.0;

    stat_lb = 0.0;
    stat_ub = 0.0;
    stat_gap = 0.0;
    nb_notopt = 0;
    best_problem = -1;

    time_lb = 0.0;
    time_eval = 0.0;

    f_mean = 0.0;
    f_var = 0.0;
                            
    lb_vector.resize(instance.getNbProblems());
    ub_vector.resize(instance.getNbProblems());
    gap_vector.resize(instance.getNbProblems());
    nb_bnb_vector.resize(instance.getNbProblems());
    
    defineLeaderBounds();
    solvers.resize(instance.getNbProblems(),NULL);
    for(int i = 0; i < instance.getNbProblems(); ++i){
        solvers[i] = new LeaderSolver(input,instance,name+" Problem "+std::to_string(i+1),i,env,general_lb,general_ub,time_prep);
    }
} 

SAALeaderSolver::~SAALeaderSolver(){
    for(int i = 0; i < instance.getNbProblems(); ++i) delete solvers[i];
    delete env;
}

void SAALeaderSolver::solve(){
    double start_time = time(NULL); 
    for(int i = 0; i < instance.getNbProblems(); ++i){
       solvers[i]->solve();
       Status curr_status = solvers[i]->getStatus(); 
       status = std::max(status,curr_status);
       if(curr_status == Status::Optimal || curr_status == Status::Time_Limit || status == Status::Iteration_Limit){
          lb_vector[i] = solvers[i]->getLB();
          ub_vector[i] = solvers[i]->getUB();
          gap_vector[i] = solvers[i]->getGAP();
          nb_bnb_vector[i] = solvers[i]->getNbBnB();
          time_lb += solvers[i]->getTime();
          if(curr_status == Status::Time_Limit) ++nb_notopt;
       }
       if(solvers[i]->getNbSol() > 0) { nb_sol = 1; testProblem(i); }
    }
    estimators();
    time_ = time(NULL) - start_time + time_prep;
}
 
void SAALeaderSolver::estimators(){
    var_lb = 0.0;
    var_nb_bnb = 0.0;
 
    long long M = instance.getNbProblems();
    lb = std::accumulate(lb_vector.begin(),lb_vector.end(),0.0)/M;
    nb_bnb = std::accumulate(nb_bnb_vector.begin(),nb_bnb_vector.end(),0.0)/M;
 
    for(int i = 0; i < instance.getNbProblems(); ++i) {
       var_lb += std::pow((lb_vector[i]-lb),2);
       var_nb_bnb += std::pow((nb_bnb_vector[i]-nb_bnb),2);
    }
    if(M > 1) var_lb = var_lb/(M*(M-1));
    if(M > 1) var_nb_bnb = var_nb_bnb/(M*(M-1));
    var_gap = var_lb + var_ub;
    gap = (ub - lb)/std::abs(ub); 
 
    // Statistical estimators considering 95% accuracy (alpha = 5 %)
    stat_gap = (ub - lb) + instance.getCriticalNormal()*std::sqrt(var_gap);             
    stat_lb = lb - instance.getCriticalTStudent()*std::sqrt(var_lb);     
    stat_ub = ub + instance.getCriticalNormal()*std::sqrt(var_ub);
}
 
void SAALeaderSolver::testProblem(int i){
    double start_eval_time = time(NULL); 
    double mean, variance, follower_mean, follower_variance;
    solvers[i]->getFollower()->evaluate(mean,variance, follower_mean, follower_variance);
    time_eval += (time(NULL) - start_eval_time);
    
    if(mean <= ub - 1e-6){
        ub = std::min(ub, mean);
        best_problem = i;
        var_ub = variance;
        
        f_mean = follower_mean;
        f_var = follower_variance;
    }
}

std::string SAALeaderSolver::write() const {
    std::string output = "";
    for(int i = 0; i < instance.getNbProblems(); ++i){
        output += std::to_string(lb_vector[i]) + std::string(";") + 
                    std::to_string(ub_vector[i]) + std::string(";") + 
                    std::to_string(gap_vector[i]) + std::string(";") +
                    std::to_string(solvers[i]->getFollower()->getFollowerOptimalObj()) + std::string(";") +
                    std::to_string(nb_bnb_vector[i]) + std::string(";");
    }

    output += std::to_string(status) + std::string(";") + 
              std::to_string(lb) + std::string(";") + 
              std::to_string(ub) + std::string(";") + 
              std::to_string(gap) + std::string(";") + 
              std::to_string(getFollower()->getFollowerOptimalObj()) + std::string(";");

    if(input.isFollowerNearOpt() == true){
        output += std::to_string(getFollower()->getFollowerNearOptimalOptObj()) + std::string(";") +
                  std::to_string(getFollower()->getFollowerNearOptimalPesObj()) + std::string(";") +
                  std::to_string(f_mean) + std::string(";") +
                  std::to_string(f_var) + std::string(";");
    }

    output += std::to_string(time_) + std::string(";") +
              std::to_string(time_lb) + std::string(";") +
              std::to_string(time_eval) + std::string(";") +
              std::to_string(nb_bnb) + std::string(";") +
              std::to_string(var_nb_bnb) + std::string(";") + 
              std::to_string(stat_lb) + std::string(";") + 
              std::to_string(stat_ub) + std::string(";") + 
              std::to_string(stat_gap) + std::string(";") + 
              std::to_string(var_lb) + std::string(";") + 
              std::to_string(var_ub) + std::string(";") + 
              std::to_string(var_gap) + std::string(";") + 
              std::to_string(nb_notopt) + std::string(";") +
              std::to_string((nb_sol > 0)) + std::string(";");

    return output;
}

std::string SAALeaderSolver::writeHead() const {
    std::string output = "";
    for(int i = 0; i < instance.getNbProblems(); ++i){
        output += "LB_" + std::to_string(i+1) + ";" +
                  "UB_" + std::to_string(i+1) + ";" +
                  "GAP_" + std::to_string(i+1) + ";" +
                  "OBJ_FOLLOWER_" + std::to_string(i+1) + ";" +
                  "NB_BNB_" + std::to_string(i+1) + ";";
    }
    output += "STATUS;LB;UB;GAP;OBJ_FOLLOWER;";
    if(input.isFollowerNearOpt() == true)
        output += "OBJ_OPT_FOLLOWER;OBJ_PES_FOLLOWER;OBJ_MEAN_FOLLOWER;OBJ_VAR_FOLLOWER;";
    output += "TIME;TIME_LB;TIME_EVAL;NB_BNB;VAR_NB_BNB;STAT_LB;STAT_UB;STAT_GAP;VAR_LB;VAR_UB;VAR_GAP;NB_NOT_OPT;NB_SOL;";
    return output;
}

void SAALeaderSolver::defineLeaderBounds(){
    double time_start = time(NULL);

    general_lb = instance.getModel()->leader_lb;
    general_ub = instance.getModel()->leader_ub;

    if((input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral) && (input.getMethodDepGeneral() == Input::MethodDepGeneral::Callback)){
        Input input_opt(input,1.0);
        LeaderSolver * leader_solver_opt = new LeaderSolver(input_opt,instance,"Leader-opt",0,env,instance.getModel()->leader_lb,instance.getModel()->leader_ub,0.0);
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++ i){
            ((LeaderSolver*)leader_solver_opt)->getX(i).set(GRB_DoubleAttr_Obj,0.0);
        }
        leader_solver_opt->solve();

        if(leader_solver_opt->getStatus() == Status::Optimal) general_lb = ((LeaderSolver*)leader_solver_opt)->getGRBModel()->get(GRB_DoubleAttr_ObjVal);
        else general_lb = std::max(general_lb,((LeaderSolver*)leader_solver_opt)->getGRBModel()->get(GRB_DoubleAttr_ObjBound));

        delete leader_solver_opt;
    }

    bool improve_bounds = false;
    if((input.getFollowerBehavior() == Input::FollowerBehavior::DepStrongWeak) && improve_bounds == true){
        Input input_opt(input,1.0);
        LeaderSolver * leader_solver_opt = new LeaderSolver(input_opt,instance,"Leader-opt",0,env,instance.getModel()->leader_lb,instance.getModel()->leader_ub,0.0);
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++ i){
            ((LeaderSolver*)leader_solver_opt)->getX(i).set(GRB_DoubleAttr_Obj,0.0);
        }
        leader_solver_opt->solve();

        if(leader_solver_opt->getStatus() == Status::Optimal) general_lb = ((LeaderSolver*)leader_solver_opt)->getGRBModel()->get(GRB_DoubleAttr_ObjVal);
        else general_lb = std::max(general_lb,((LeaderSolver*)leader_solver_opt)->getGRBModel()->get(GRB_DoubleAttr_ObjBound));

        delete leader_solver_opt;

        Input input_pes(input,1.0);
        LeaderSolver * leader_solver_pes = new LeaderSolver(input_pes,instance,"Leader-pes",0,env,instance.getModel()->leader_lb,instance.getModel()->leader_ub,0.0);
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++ i){
            ((LeaderSolver*)leader_solver_pes)->getX(i).set(GRB_DoubleAttr_Obj,0.0);
        }
        ((LeaderSolver*)leader_solver_pes)->getGRBModel()->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        leader_solver_pes->solve();

        if(leader_solver_pes->getStatus() == Status::Optimal) general_ub = ((LeaderSolver*)leader_solver_pes)->getGRBModel()->get(GRB_DoubleAttr_ObjVal);
        else general_ub = std::min(general_ub,((LeaderSolver*)leader_solver_pes)->getGRBModel()->get(GRB_DoubleAttr_ObjBound));

        delete leader_solver_pes;
    }

    std::cout << "General LB: " << instance.getModel()->leader_lb << " " << general_lb << std::endl;
    std::cout << "General UB: " << instance.getModel()->leader_ub << " " << general_ub << std::endl;

    time_prep = time(NULL) - time_start;
}