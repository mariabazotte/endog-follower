#include "leadersolver.hpp"
#include "../followersolver/followersolver.hpp" 

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

    // Define output string.
    std::string output = "";

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
    instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::Neutral,0,0.0);
    writeCompDepGeneral(output,Input::TypesDepGeneral::Neutral,0,0.0,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);

    // Proportional.
    for(int nb_int: instance.getDepGeneralIntConfigs()){
        instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenProportional,nb_int,1.0);
        writeCompDepGeneral(output,Input::TypesDepGeneral::GenProportional,nb_int,1.0,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
    }

    // Fragile. 
    for(int nb_int: instance.getDepGeneralIntConfigs()){
        for(double scal: instance.getDepGeneralScalFrgConfigs()){
            instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenFragile,nb_int,scal);
            writeCompDepGeneral(output,Input::TypesDepGeneral::GenFragile,nb_int,scal,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
        }
    }

    // Fragile Power. 
    for(int nb_int: instance.getDepGeneralIntConfigs()){
        for(double scal: instance.getDepGeneralScalFrgPowConfigs()){
            instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenFragilePower,nb_int,scal);
            writeCompDepGeneral(output,Input::TypesDepGeneral::GenFragilePower,nb_int,scal,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
        }
    }

    // Strong Power. 
    for(int nb_int: instance.getDepGeneralIntConfigs()){
        for(double scal: instance.getDepGeneralScalStrPowConfigs()){
            instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,x_,yi_,Fs_,Fw_,fs_,Input::TypesDepGeneral::GenStrongPower,nb_int,scal);
            writeCompDepGeneral(output,Input::TypesDepGeneral::GenStrongPower,nb_int,scal,(x_leader_obj+eval),var_eval,ess,rhat,fs_,f_eval,f_var_eval);
        }
    }

    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

LeaderSolver::LeaderSolver(const Input & input, Instance & instance, std::string name) : 
                                            LeaderSolver(input,instance,name,-1,new GRBEnv()){} 

LeaderSolver::LeaderSolver(const Input & input, Instance & instance, std::string name, int pr, GRBEnv *env) : 
                                AbstractLeaderSolver(input,instance,name,env), pr(pr), model(new GRBModel(*env)) {
    // Create leader problem
    params();
    create();
    
    // Create follower problem
    FollowerSolverFactory factory;
    follower = factory.createFollowerSolver(input,instance,pr,this);
}

LeaderSolver::~LeaderSolver(){
    delete follower;
    if(x_) delete[] x_;
    if(x) delete[] x;
    delete model;
    if(pr == -1) delete env;
}

void LeaderSolver::params(){
    model->set(GRB_IntParam_Threads, input.getNbThreads());
    model->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
    model->set(GRB_DoubleParam_MIPGap, 1e-4);
    model->set(GRB_DoubleParam_IntFeasTol, 1e-5);
    model->set(GRB_DoubleParam_FeasibilityTol, 1e-6);
    // model->set(GRB_IntParam_Presolve, 0);
    if(input.getVerbose() < 1) model->set(GRB_IntParam_OutputFlag, 0);
    else model->set(GRB_IntParam_OutputFlag, 1);

    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral && 
       input.useLazyCallback() == true){
        model->set(GRB_IntParam_LazyConstraints, 1);
    }
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
    gap = (ub-lb)/ub;
    if(nb_sol > 0){
        if(x_) delete[] x_;
        x_ = model->get(GRB_DoubleAttr_X,x,instance.getModel()->nb_leader_vars);
    }
    follower->upd_solution();
}

void LeaderSolver::solve(){
    double start_time = time(NULL); 

    if(input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral && input.useLazyCallback() == true) 
        model->setCallback(((DepGeneralFollowerSolver*)getFollower()));

    model->update();
    model->optimize();

    status = statusFromGurobi(model->get(GRB_IntAttr_Status));
    if(status == Status::Optimal || status == Status::Time_Limit) upd_solution();
    else verify_stat(); 

    time_ = time(NULL) - start_time;
    if(input.getVerbose() >= 1) print();

    // Follower strong, weak, and interior solutions.
    getFollower()->computeStrongWeakInteriorSolutions();
}

std::string LeaderSolver::write() const {
    std::string output = std::to_string(status) + std::string(";") + 
                        std::to_string(lb) + std::string(";") + 
                        std::to_string(ub) + std::string(";") + 
                        std::to_string(gap) + std::string(";") + 
                        std::to_string(follower->getFollowerOptimalObj()) + std::string(";");
    
    if(input.isFollowerNearOpt() == true){
        double m, v, fm, fv;
        getFollower()->evaluate(m, v, fm, fv); // fm is the expected value of follower objective. Variance here should be zero.
        output += std::to_string(follower->getFollowerNearOptimalOptObj()) + std::string(";") +
                  std::to_string(follower->getFollowerNearOptimalPesObj()) + std::string(";") +
                  std::to_string(fm) + std::string(";");
    }
    
    output += std::to_string(time_) + std::string(";") +
              std::to_string(nb_bnb) + std::string(";");
    return output;
}

std::string LeaderSolver::writeHead() const {
    std::string output = "STATUS;LB;UB;GAP;OBJ_FOLLOWER;";
    if(input.isFollowerNearOpt() == true)
        output += "OBJ_OPT_FOLLOWER;OBJ_PES_FOLLOWER;OBJ_MEAN_FOLLOWER;";
    output += "TIME;NB_BNB;";
    return output;
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

    solvers.resize(instance.getNbProblems(),NULL);
    for(int i = 0; i < instance.getNbProblems(); ++i){
        solvers[i] = new LeaderSolver(input,instance,name+" Problem "+std::to_string(i+1),i,env);
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
       if(solvers[i]->getNbSol() > 0) testProblem(i);
    }
    estimators();
    time_ = time(NULL) - start_time;
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
    gap = (ub - lb)/ub;  
 
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
    
    if(mean < ub){
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
              std::to_string(nb_notopt) + std::string(";");

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
    output += "TIME;TIME_LB;TIME_EVAL;NB_BNB;VAR_NB_BNB;STAT_LB;STAT_UB;STAT_GAP;VAR_LB;VAR_UB;VAR_GAP;NB_NOT_OPT;";
    return output;
}