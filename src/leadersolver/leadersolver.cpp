#include "leadersolver.hpp"
#include "../followersolver/followersolver.hpp" 

std::string AbstractLeaderSolver::writeHeadComp() const {
    return "FOLLOWER_BEHAVIOR;COOPERATION_TYPE;PARAMS;VALUE;";
}

std::string AbstractLeaderSolver::writeComp() const {
    // Define part of objective involving only leader variables.
    double obj = 0.0;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        obj += instance.getModel()->leader_vars[i].obj_leader*getX_(i);

    // Compute follower strong and weak solutions for the fixed leader decision.
    getFollower()->computeStrongWeakInteriorSolutions();

    double leader_val = 0.0;
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i)
        leader_val += instance.getModel()->follower_vars[i].obj_leader*getFollower()->getYi_(i);

    std::cout << "leader: " << leader_val << std::endl;

    std::string output = "";

    // Fixed Strong-Weak configurations.
    for(double beta: instance.getFixStrongWeakCoopConfigs()){
        double val = beta*getFollower()->getLeaderOptObj() + (1.0-beta)*getFollower()->getLeaderPesObj();
        
        output += std::to_string(Input::FollowerBehavior::FixStrongWeak) + ";-;" + 
                  std::to_string(beta) + ";" + std::to_string((obj+val)) + ";\n";
    }

    // Dependent Strong Weak configuations.

    // Proportional
    double beta = instance.getEvalDepStrongWeakProbab(getFollower()->getFollowerOptimalObj(),
                                                Input::TypesDepStrongWeak::Proportional,0.0);
    double val = beta*getFollower()->getLeaderOptObj() + (1.0-beta)*getFollower()->getLeaderPesObj();

    output += std::to_string(Input::FollowerBehavior::DepStrongWeak) + ";" + 
              std::to_string(Input::TypesDepStrongWeak::Proportional) + ";-;" + 
              std::to_string((obj+val)) + ";\n";


    // Threshold
    for(double param : instance.getDepStrongWeakThrConfigs()){
        beta = instance.getEvalDepStrongWeakProbab(getFollower()->getFollowerOptimalObj(),
                                                Input::TypesDepStrongWeak::Threshold,param);
        val = beta*getFollower()->getLeaderOptObj() + (1.0-beta)*getFollower()->getLeaderPesObj();

        output += std::to_string(Input::FollowerBehavior::DepStrongWeak) + ";" + 
                std::to_string(Input::TypesDepStrongWeak::Threshold) + ";" +
                std::to_string(param) + ";" + std::to_string((obj+val)) + ";\n";
    }

    // Strong
    for(double param : instance.getDepStrongWeakStrConfigs()){
        beta = instance.getEvalDepStrongWeakProbab(getFollower()->getFollowerOptimalObj(),
                                                Input::TypesDepStrongWeak::Strong,param);
        val = beta*getFollower()->getLeaderOptObj() + (1.0-beta)*getFollower()->getLeaderPesObj();

        output += std::to_string(Input::FollowerBehavior::DepStrongWeak) + ";" + 
                std::to_string(Input::TypesDepStrongWeak::Strong) + ";" +
                std::to_string(param) + ";" + std::to_string((obj+val)) + ";\n";
    }

    // Fragile
    for(double param : instance.getDepStrongWeakFrgConfigs()){
        beta = instance.getEvalDepStrongWeakProbab(getFollower()->getFollowerOptimalObj(),
                                                Input::TypesDepStrongWeak::Fragile,param);
        val = beta*getFollower()->getLeaderOptObj() + (1.0-beta)*getFollower()->getLeaderPesObj();

        output += std::to_string(Input::FollowerBehavior::DepStrongWeak) + ";" + 
                std::to_string(Input::TypesDepStrongWeak::Fragile) + ";" +
                std::to_string(param) + ";" + std::to_string((obj+val)) + ";\n";
    }

    // Dependent General configuations.

    // Uniform
    double eval, var_eval, f_eval, f_var_eval;
    instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,
        getX_(),getFollower()->getYi_(),getFollower()->getFollowerOptimalObj(),
        Input::TypesDepGeneral::Neutral, 0, 0.0);

    output += std::to_string(Input::FollowerBehavior::DepGeneral) + ";" + 
                std::to_string(Input::TypesDepGeneral::Neutral) + ";-;" + 
                std::to_string((obj+val)) + ";\n";

    std::cout << "---------------------" << std::endl;

    exit(0);

    // Proportional 
    for(int nb_int: instance.getDepGeneralIntConfigs()){
        for(double scal: instance.getDepGeneralScalConfigs()){
            double eval, var_eval, f_eval, f_var_eval;
            instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,
                getX_(),getFollower()->getYi_(),getFollower()->getFollowerOptimalObj(),
                Input::TypesDepGeneral::GenProportional, nb_int, scal);

            output += std::to_string(Input::FollowerBehavior::DepGeneral) + ";" + 
                    std::to_string(Input::TypesDepGeneral::GenProportional) + ";" +
                    std::to_string(nb_int) + "/" + std::to_string(scal) + ";" + 
                    std::to_string((obj+val)) + ";\n";
            break;
        }
        break;
    }

    std::cout << "---------------------" << std::endl;

    // Strong Fragile 
    for(int nb_int: instance.getDepGeneralIntConfigs()){
        for(double scal: instance.getDepGeneralScalConfigs()){
            double eval, var_eval, f_eval, f_var_eval;
            instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,
                getX_(),getFollower()->getYi_(),getFollower()->getFollowerOptimalObj(),
                Input::TypesDepGeneral::StrongFragile, nb_int, scal);

            output += std::to_string(Input::FollowerBehavior::DepGeneral) + ";" + 
                    std::to_string(Input::TypesDepGeneral::StrongFragile) + ";" +
                    std::to_string(nb_int) + "/" + std::to_string(scal) + ";" + 
                    std::to_string((obj+val)) + ";\n";
            break;
        }
        break;
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
    if(x) delete[] x;
    if(x_) delete[] x_;
    delete model;
}

void LeaderSolver::params(){
    model->set(GRB_IntParam_Threads, input.getNbThreads());
    model->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
    model->set(GRB_DoubleParam_MIPGap, 1e-4);
    model->set(GRB_DoubleParam_IntFeasTol, 1e-5);
    model->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
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
    gap = (ub-lb)/ub;
    if(nb_sol > 0){
        if(x_) delete[] x_;
        x_ = model->get(GRB_DoubleAttr_X,x,instance.getModel()->nb_leader_vars);
    }
    follower->upd_solution();
}

void LeaderSolver::solve(){
    double start_time = time(NULL); 

    model->update();
    model->optimize();

    status = statusFromGurobi(model->get(GRB_IntAttr_Status));
    if(status == Status::Optimal || status == Status::Time_Limit) upd_solution();
    else verify_stat(); 

    time_ = time(NULL) - start_time;
    if(input.getVerbose() >= 1) print();
}

std::string LeaderSolver::write() const {
    std::string output = std::to_string(status) + std::string(";") + 
                        std::to_string(lb) + std::string(";") + 
                        std::to_string(ub) + std::string(";") + 
                        std::to_string(gap) + std::string(";") + 
                        std::to_string(follower->getFollowerOptimalObj()) + std::string(";");
    
    if(input.isFollowerNearOpt() == true){
        output += std::to_string(follower->getFollowerNearOptimalOptObj()) + std::string(";") +
                  std::to_string(follower->getFollowerNearOptimalPesObj()) + std::string(";");
    }
    
    output += std::to_string(time_) + std::string(";") +
              std::to_string(nb_bnb) + std::string(";");
    return output;
}

std::string LeaderSolver::writeHead() const {
    std::string output = "STATUS;LB;UB;GAP;OBJ_FOLLOWER;";
    if(input.isFollowerNearOpt() == true)
        output += "OBJ_OPT_FOLLOWER;OBJ_PES_FOLLOWER;";
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
    double eval = solvers[i]->getFollower()->evaluate();
    time_eval += (time(NULL) - start_eval_time);
    
    if(eval < ub){
        ub = std::min(ub, eval);
        best_problem = i;

        if(input.getFollowerBehavior() == Input::FollowerBehavior::DepGeneral){
            var_ub = 0.0;
            for(int s = 0; s < instance.getNbValidateScenarios(); ++s)
                var_ub += std::pow(solvers[i]->getFollower()->getDiffEvalScenario(s),2);
            long long N = instance.getNbValidateScenarios();
            if(N > 1) var_ub = var_ub/(N*(N-1));
        }
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
                  std::to_string(getFollower()->getFollowerNearOptimalPesObj()) + std::string(";");
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
        output += "OBJ_OPT_FOLLOWER;OBJ_PES_FOLLOWER;";
    output += "TIME;TIME_LB;TIME_EVAL;NB_BNB;VAR_NB_BNB;STAT_LB;STAT_UB;STAT_GAP;VAR_LB;VAR_UB;VAR_GAP;NB_NOT_OPT;";
    return output;
}