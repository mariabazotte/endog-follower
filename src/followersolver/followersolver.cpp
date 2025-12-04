#include "followersolver.hpp"
#include "../leadersolver/leadersolver.hpp"

void AbstractFollowerSolver::definePrimalVars(GRBVar * & y, std::string v){
    y = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        
        if(var.lb <= -instance.getModel()->inf) 
            y[i].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
        else y[i].set(GRB_DoubleAttr_LB,var.lb);
        if(var.ub >=  instance.getModel()->inf) 
            y[i].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        else y[i].set(GRB_DoubleAttr_UB,var.ub);

        y[i].set(GRB_StringAttr_VarName,var.name + v);
        y[i].set(GRB_CharAttr_VType,var.type);
    }
}

void AbstractFollowerSolver::definePrimalConstrs(GRBVar * & y, GRBVar * & slacks, GRBVar * & slack_lbs, GRBVar * & slack_ubs, bool int_problem, std::string v){
    // Define bound slack variables.
    slack_lbs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    slack_ubs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];

        if(!(var.lb <= -instance.getModel()->inf)){
            if(int_problem == false) leader->getGRBModel()->addConstr(slack_lbs[i] == y[i] - var.lb);
            else leader->getGRBModel()->addConstr(slack_lbs[i] + si == y[i] - var.lb);
        }if(!(var.ub >=  instance.getModel()->inf)){
            if(int_problem == false) leader->getGRBModel()->addConstr(slack_ubs[i] == var.ub - y[i]);
            else leader->getGRBModel()->addConstr(slack_ubs[i] + si == var.ub - y[i]);
        }
    }

    // Primal constraints and their slack variables.
    slacks =  leader->getGRBModel()->addVars(instance.getModel()->nb_follower_constrs,GRB_CONTINUOUS);
    for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
        BilevelConstraint constr = instance.getModel()->follower_constrs[c];
        
        GRBLinExpr expr = 0;
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
            BilevelVariable var = instance.getModel()->leader_vars[i];
            expr += constr.coeffs[var.id]*leader->getX(i); 
        }
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            BilevelVariable var = instance.getModel()->follower_vars[i];
            expr += constr.coeffs[var.id]*y[i];
        }

        // Constraints and Slacks.
        if(constr.sense == '=') 
            leader->getGRBModel()->addConstr(expr == constr.rhs, constr.name + v);
        else if(constr.sense == '<'){
            if(int_problem == false) leader->getGRBModel()->addConstr(expr + slacks[c] == constr.rhs, constr.name + v);
            else leader->getGRBModel()->addConstr(expr + slacks[c] + si == constr.rhs, constr.name + v);
        }else if(constr.sense == '>'){
            if(int_problem == false) leader->getGRBModel()->addConstr(expr - slacks[c] == constr.rhs, constr.name + v);
            else leader->getGRBModel()->addConstr(expr - slacks[c] - si == constr.rhs, constr.name + v);
        }
    }
}

void AbstractFollowerSolver::definePrimalConstrs(GRBVar * & y, std::string v){
    // Define primal constraints.
    for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
        BilevelConstraint constr = instance.getModel()->follower_constrs[c];
        
        GRBLinExpr expr = 0;
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
            BilevelVariable var = instance.getModel()->leader_vars[i];
            expr += constr.coeffs[var.id]*leader->getX(i); 
        }
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            BilevelVariable var = instance.getModel()->follower_vars[i];
            expr += constr.coeffs[var.id]*y[i];
        }

        // Constraints.
        if(constr.sense == '=') 
            leader->getGRBModel()->addConstr(expr == constr.rhs, constr.name + v);
        else if(constr.sense == '<'){
            leader->getGRBModel()->addConstr(expr <= constr.rhs, constr.name + v);
        }else if(constr.sense == '>'){
            leader->getGRBModel()->addConstr(expr >= constr.rhs, constr.name + v);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void AbstractFollowerSolver::defineDualVars(GRBVar * & ds, GRBVar * & dlbs, GRBVar * & dubs, bool min_problem, bool int_problem){
    // Dual variables.
    int nb_constrs = instance.getModel()->nb_follower_constrs;
    if(min_problem == false) ++nb_constrs;  // Optimality or near-optimality constraint.
    if(int_problem == true) ++nb_constrs;   // Non-negativity of slack variable.
    
    ds = leader->getGRBModel()->addVars(nb_constrs,GRB_CONTINUOUS);
    for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
        BilevelConstraint constr = instance.getModel()->follower_constrs[c];

        if(constr.sense == '=') {
            ds[c].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
            ds[c].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        }else if(constr.sense == '<'){
            if(min_problem == true){
                ds[c].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
                ds[c].set(GRB_DoubleAttr_UB,0.0);
            }else{
                ds[c].set(GRB_DoubleAttr_LB,0.0);
                ds[c].set(GRB_DoubleAttr_UB,GRB_INFINITY);
            }
        }else if(constr.sense == '>' ){
            if(min_problem == true){
                ds[c].set(GRB_DoubleAttr_LB,0.0);
                ds[c].set(GRB_DoubleAttr_UB,GRB_INFINITY);
            }else{
                ds[c].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
                ds[c].set(GRB_DoubleAttr_UB,0.0);
            }
        }
    }
    if(min_problem == false){ 
        //  Maximization problem.
        if(input.isFollowerNearOpt() == true){
            // Optimality constraint is <=.
            ds[instance.getModel()->nb_follower_constrs].set(GRB_DoubleAttr_LB,0.0);
            ds[instance.getModel()->nb_follower_constrs].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        }else{
            // Optimality constraint is ==. 
            ds[instance.getModel()->nb_follower_constrs].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
            ds[instance.getModel()->nb_follower_constrs].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        }
    }
    if(min_problem == false && int_problem == true){
        // Maximization problem. Slack non-negativity constraint (>=).
        ds[instance.getModel()->nb_follower_constrs+1].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
        ds[instance.getModel()->nb_follower_constrs+1].set(GRB_DoubleAttr_UB,0.0);
    } 

    // Dual variables for bound constraints.
    dlbs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    dubs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
}

void AbstractFollowerSolver::defineDualConstrs(GRBVar * & ds, GRBVar * & dlbs, GRBVar * & dubs, bool min_problem, bool int_problem){
    // Dual constraints for primal variables.
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];

        // Consider coefficient all primal constraints.
        GRBLinExpr expr = 0;
        for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
            BilevelConstraint constr = instance.getModel()->follower_constrs[c];
            expr += constr.coeffs[var.id]*ds[c];
        }

        // Pessimistic and interior cases we should also consider optimality constraint.
        if(min_problem == false) expr += var.obj_follower*ds[instance.getModel()->nb_follower_constrs];

        // Optimistic problem: guarantee min objective follower.
        if(min_problem == true){
            if(var.lb <= -instance.getModel()->inf && var.ub >= instance.getModel()->inf){
                leader->getGRBModel()->addConstr(expr == var.obj_follower);
            }else if(var.lb <= -instance.getModel()->inf){ 
                leader->getGRBModel()->addConstr(expr - dubs[i] == var.obj_follower);
            }else if(var.ub >= instance.getModel()->inf){
                leader->getGRBModel()->addConstr(expr + dlbs[i] == var.obj_follower);
            }else{
                leader->getGRBModel()->addConstr(expr + dlbs[i] - dubs[i] == var.obj_follower);
            }
        }
        // Pessimistic problem: guarantee max objective leader.
        if(min_problem == false && int_problem == false){
            if(var.lb <= -instance.getModel()->inf && var.ub >= instance.getModel()->inf){
                leader->getGRBModel()->addConstr(expr == var.obj_leader);
            }else if(var.lb <= -instance.getModel()->inf){ 
                leader->getGRBModel()->addConstr(expr + dubs[i] == var.obj_leader);
            }else if(var.ub >= instance.getModel()->inf){
                leader->getGRBModel()->addConstr(expr - dlbs[i] == var.obj_leader);
            }else{
                leader->getGRBModel()->addConstr(expr - dlbs[i] + dubs[i] == var.obj_leader);
            } 
        }
        // Interior problem: guarantee max slack.
        if(min_problem == false && int_problem == true){
            if(var.lb <= -instance.getModel()->inf && var.ub >= instance.getModel()->inf){ 
                leader->getGRBModel()->addConstr(expr == 0.0);
            }else if(var.lb <= -instance.getModel()->inf){ 
                leader->getGRBModel()->addConstr(expr + dubs[i] == 0.0);
            }else if(var.ub >= instance.getModel()->inf){
                leader->getGRBModel()->addConstr(expr - dlbs[i] == 0.0);
            }else{
                leader->getGRBModel()->addConstr(expr - dlbs[i] + dubs[i] == 0.0);
            } 
        }
    }
    // Dual constraint for slack variable (si >= 0) for interior solution.
    if(min_problem == false && int_problem == true){
        GRBLinExpr expr = 0;
        for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
            BilevelConstraint constr = instance.getModel()->follower_constrs[c];
            if(constr.sense == '<') expr += ds[c];
            else if(constr.sense == '>') expr += -ds[c];
        }
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            BilevelVariable var = instance.getModel()->follower_vars[i];
            if(!(var.lb <= -instance.getModel()->inf)) expr += dlbs[i];
            if(!(var.ub >=  instance.getModel()->inf)) expr += dubs[i];
        }
        expr += ds[instance.getModel()->nb_follower_constrs+1];
        leader->getGRBModel()->addConstr(expr == 1.0);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void AbstractFollowerSolver::defineCompConstrs(GRBVar * & slacks, GRBVar * & slack_lbs, GRBVar * & slack_ubs, 
                                                GRBVar * & ds, GRBVar * & dlbs, GRBVar * & dubs){
    
    // Complementary constraints for primal constraints. 
    for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
        BilevelConstraint constr = instance.getModel()->follower_constrs[c];

        if(constr.sense != '='){
            GRBVar v[2] = {slacks[c], ds[c]};
            double w[2] = {1.0, 2.0};

            // Add SOS1 constraint.
            leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);
        }
    }

    // Complementary constraints for primal bound constraints. 
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];

        if(!(var.lb <= -instance.getModel()->inf)){
            GRBVar vlb[2] = {slack_lbs[i], dlbs[i]};
            double wlb[2] = {1.0, 2.0};

            // Add SOS1 constraints.
            leader->getGRBModel()->addSOS(vlb, wlb, 2, GRB_SOS_TYPE1);
        }

        if(!(var.ub >= instance.getModel()->inf)){
            GRBVar vub[2] = {slack_ubs[i], dubs[i]};
            double wub[2] = {1.0, 2.0};

            // Add SOS1 constraints.
            leader->getGRBModel()->addSOS(vub, wub, 2, GRB_SOS_TYPE1);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void AbstractFollowerSolver::defineFollowerOptimalObj(GRBVar * & y, std::string v){
    // Define follower optimal objective value.
    if(fs_not_init == true){
        fs = leader->getGRBModel()->addVar(-GRB_INFINITY,GRB_INFINITY,0.0,GRB_CONTINUOUS,"fs");
        fs_not_init = false;
    }
    GRBLinExpr expr = 0;
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        expr += var.obj_follower*y[i];
    }
    leader->getGRBModel()->addConstr(expr == fs,"Follower_Optimal_Objective" + v);
}

void AbstractFollowerSolver::defineFollowerNearOptimalObj(GRBVar & f_eps, GRBVar * & y_eps, GRBVar & slack_eps, std::string v){
    slack_eps = leader->getGRBModel()->addVar(0.0,GRB_INFINITY,0.0,GRB_CONTINUOUS);
    f_eps = leader->getGRBModel()->addVar(-GRB_INFINITY,GRB_INFINITY,0.0,GRB_CONTINUOUS,"f_eps" + v);
    GRBLinExpr expr = 0;
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        expr += var.obj_follower*y_eps[i];
    }
    leader->getGRBModel()->addConstr(f_eps == expr, "Follower_Eps_Objective" + v);
    if(input.isNearOptMult() == true)
        leader->getGRBModel()->addConstr(expr + slack_eps == ((1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub) ,"Constr_Follower_Eps_Objective" + v);
    else
        leader->getGRBModel()->addConstr(expr + slack_eps == fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) ,"Constr_Follower_Eps_Objective" + v); 
}

void AbstractFollowerSolver::defineLeaderObj(GRBVar & F, GRBVar * & y, std::string v){
    // Include variable with leader obj value related to only follower variables. 
    // (Define it according to optimality or near optimality)
    F = leader->getGRBModel()->addVar(-GRB_INFINITY,GRB_INFINITY,0.0,GRB_CONTINUOUS,"F" + v);

    GRBLinExpr expr = 0;
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        expr += var.obj_leader*y[i];
    }
    leader->getGRBModel()->addConstr(F == expr,"Leader_Objective" + v);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void AbstractFollowerSolver::defineOptimalFollower(){
    // Define primal variables.
    definePrimalVars(ys,"_s");

    // Define primal constraints and slacks.
    GRBVar * slacks = NULL;    
    GRBVar * slack_lbs = NULL;
    GRBVar * slack_ubs = NULL;
    definePrimalConstrs(ys,slacks,slack_lbs,slack_ubs,false,"_s");
    
    // Define dual variables. (primal is a minimization problem)
    GRBVar * ds = NULL;
    GRBVar * dlbs = NULL;
    GRBVar * dubs = NULL;
    defineDualVars(ds,dlbs,dubs,true,false);

    // Define dual constraints. (primal is a minimization problem)
    defineDualConstrs(ds,dlbs,dubs,true,false);

    // Define complementarity constraints.
    defineCompConstrs(slacks,slack_lbs,slack_ubs,ds,dlbs,dubs);

    // Define follower optimal value.
    defineFollowerOptimalObj(ys,"_s");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void AbstractFollowerSolver::defineOptimisticFollower(){
    // Define follower's optimality.
    defineOptimalFollower();

    // Define auxiliar primal constraints in the case of near optimality.
    if(input.isFollowerNearOpt() == true){
        // Define primal eps near optimal variables.
        definePrimalVars(ys_eps,"_s_eps");

        // Define primal constraints.
        definePrimalConstrs(ys_eps,"_s_eps");

        // Define near optimality constraint.
        GRBVar slack_eps;
        defineFollowerNearOptimalObj(fs_eps,ys_eps,slack_eps,"_s_eps");
    }

    // Define leader optimal value.
    if(input.isFollowerNearOpt() == true) 
        defineLeaderObj(Fs,ys_eps,"_s_eps");
    else defineLeaderObj(Fs,ys,"_s");
}

void AbstractFollowerSolver::definePessimisticFollower(){
    // Define primal variables.
    definePrimalVars(yw,"_w");

    // Define primal constraints and slacks.
    GRBVar * slacks =  NULL;
    GRBVar * slack_lbs = NULL;
    GRBVar * slack_ubs = NULL;
    definePrimalConstrs(yw,slacks,slack_lbs,slack_ubs,false,"_w");

    // Define primal constraint and slack corresponding to optimality, specifically for pessimistic problem.
    // Define it according to optimality or near optimality.
    GRBVar slack_eps;
    if(input.isFollowerNearOpt() == true)
        defineFollowerNearOptimalObj(fw_eps,yw,slack_eps,"_w_eps");
    else defineFollowerOptimalObj(yw,"_w");

    // Define dual variables. (primal is a maximization problem).
    GRBVar * ds = NULL;
    GRBVar * dlbs = NULL;
    GRBVar * dubs = NULL;
    defineDualVars(ds,dlbs,dubs,false,false);

    // Define dual constraints. (primal is a maximization problem).
    defineDualConstrs(ds,dlbs,dubs,false,false);
    
    // Complementarity constraints.
    defineCompConstrs(slacks,slack_lbs,slack_ubs,ds,dlbs,dubs);

    // Complementary constraints for optimality constraint.
    if(input.isFollowerNearOpt() == true){
        GRBVar v[2] = {slack_eps, ds[instance.getModel()->nb_follower_constrs]};
        double w[2] = {1.0, 2.0};
        // Add SOS1 constraint.
        leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);
    }

    // Leader objective value.
    defineLeaderObj(Fw,yw,"_pes");    
}

void AbstractFollowerSolver::defineInteriorFollower(){
    // Define follower's optimality.
    defineOptimalFollower();

    // Define primal interior variables.
    definePrimalVars(yi,"_i");
    si = leader->getGRBModel()->addVar(0.0,GRB_INFINITY,0.0,GRB_CONTINUOUS,"s_i");

    // Define primal constraints and slacks.
    GRBVar * slacks =  NULL;
    GRBVar * slack_lbs = NULL;
    GRBVar * slack_ubs = NULL;
    definePrimalConstrs(yi,slacks,slack_lbs,slack_ubs,true,"_i");

    // Define primal constraint and slack corresponding to optimality.
    // Define it according to optimality or near optimality.
    GRBVar slack_eps;
    if(input.isFollowerNearOpt() == true)
        defineFollowerNearOptimalObj(fi_eps,yi,slack_eps,"_i_eps");
    else defineFollowerOptimalObj(yi,"_i");

    // Define dual variables. (primal is a maximization problem).
    GRBVar * ds = NULL;
    GRBVar * dlbs = NULL;
    GRBVar * dubs = NULL;
    defineDualVars(ds,dlbs,dubs,false,true);

    // Define dual constraints. (primal is a maximization problem).
    defineDualConstrs(ds,dlbs,dubs,false,true);

    // Complementarity constraints.
    defineCompConstrs(slacks,slack_lbs,slack_ubs,ds,dlbs,dubs);

    // Complementary constraints for optimality constraint.
    if(input.isFollowerNearOpt() == true){
        GRBVar v[2] = {slack_eps, ds[instance.getModel()->nb_follower_constrs]};
        double w[2] = {1.0, 2.0};
        // Add SOS1 constraint.
        leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);
    }

    // Complementary constraint for slack variable (si >= 0).
    GRBVar v[2] = {si, ds[instance.getModel()->nb_follower_constrs+1]};
    double w[2] = {1.0, 2.0};
    // Add SOS1 constraint.
    leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void AbstractFollowerSolver::upd_solution(){
    if(leader->getNbSol() > 0){
        if(ys_) delete[] ys_;
        if(yw_) delete[] yw_;
        if(yi_) delete[] yi_;
        if(ys_eps_) delete[] ys_eps_;

        if(ys) ys_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,ys,instance.getModel()->nb_follower_vars);
        if(yw) yw_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,yw,instance.getModel()->nb_follower_vars);
        if(yi) yi_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,yi,instance.getModel()->nb_follower_vars);
        if(ys_eps) ys_eps_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,ys_eps,instance.getModel()->nb_follower_vars);

        Fs_ = Fs.get(GRB_DoubleAttr_X);
        Fw_ = Fw.get(GRB_DoubleAttr_X);

        std::cout << "Fs: " << Fs_ << " Fw: " << Fw_ << std::endl;

        fs_ = fs.get(GRB_DoubleAttr_X);
        std::cout << "fs: " << fs_ << std::endl;
        if(input.isFollowerNearOpt() == true) {
            fs_eps_ = fs_eps.get(GRB_DoubleAttr_X);
            fw_eps_ = fw_eps.get(GRB_DoubleAttr_X);
            std::cout << "fs: " << fs_eps_ << " fw: " << fw_eps_ << std::endl;

            long double obj_follower_pes = 0.0;
            long double obj_follower_opt = 0.0;
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                std::cout << ys_[i] << " " << ys_eps_[i] << " " << yw_[i] << std::endl;
                BilevelVariable var = instance.getModel()->follower_vars[i];
                
                obj_follower_opt += var.obj_follower*ys_eps_[i];
                obj_follower_pes += var.obj_follower*yw_[i];
            }
            std::cout << "fs: " << obj_follower_opt << " fw: " << obj_follower_pes << std::endl;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

// Compute strong, weak, and interior follower solutions for current leader decision if needed to evaluate other follower behaviors.
void AbstractFollowerSolver::computeStrongWeakInteriorSolutions(bool compute_str, bool compute_wk, bool compute_int){
    
    GRBModel * model_eval = new GRBModel(*leader->getGRBEnv());

    // Follower variables.
    GRBVar * y_eval = model_eval->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        
        if(var.lb <= -instance.getModel()->inf) 
            y_eval[i].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
        else y_eval[i].set(GRB_DoubleAttr_LB,var.lb);
        if(var.ub >=  instance.getModel()->inf) 
            y_eval[i].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        else y_eval[i].set(GRB_DoubleAttr_UB,var.ub);

        y_eval[i].set(GRB_StringAttr_VarName,var.name);
        y_eval[i].set(GRB_CharAttr_VType,var.type);

        // Define leader objective.
        y_eval[i].set(GRB_DoubleAttr_Obj,var.obj_leader);
    }

    // Follower constraints.
    std::vector<GRBConstr> grb_constrs;
    for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
        BilevelConstraint constr = instance.getModel()->follower_constrs[c];

        GRBLinExpr expr = 0;
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
            BilevelVariable var = instance.getModel()->leader_vars[i];
            expr += constr.coeffs[var.id]*leader->getX_(i); 
        }

        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            BilevelVariable var = instance.getModel()->follower_vars[i];
            expr += constr.coeffs[var.id]*y_eval[i];
        }

        if(constr.sense == '=') 
            grb_constrs.push_back(model_eval->addConstr(expr == constr.rhs, constr.name));
        else if(constr.sense == '<'){
            grb_constrs.push_back(model_eval->addConstr(expr <= constr.rhs, constr.name));
        }else if(constr.sense == '>'){
            grb_constrs.push_back(model_eval->addConstr(expr >= constr.rhs, constr.name));
        }
    }

    // Follower optimality constraint.
    GRBVar f_eval = model_eval->addVar(-GRB_INFINITY,GRB_INFINITY,0.0,GRB_CONTINUOUS);
    GRBLinExpr expr = 0;
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        expr += var.obj_follower*y_eval[i];
    }
    model_eval->addConstr(expr == f_eval);
    if(input.isFollowerNearOpt() == true) {  
        if(input.isNearOptMult() == true)
            grb_constrs.push_back(model_eval->addConstr(expr <= (1.0 - input.getEpsNearOpt())*fs_ + input.getEpsNearOpt()*instance.getModel()->follower_ub,"Follower_Objective"));
        else
            grb_constrs.push_back(model_eval->addConstr(expr <= fs_ + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) ,"Follower_Objective"));
    }else{
        model_eval->addConstr(expr == fs_,"Follower_Objective");
    }

    // Solve optimistic problem.
    if(compute_str == true){
        model_eval->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        model_eval->update();
        model_eval->optimize();

        Status status_eval = statusFromGurobi(model_eval->get(GRB_IntAttr_Status));
        if(status_eval == Status::Optimal){
            Fs_ = model_eval->get(GRB_DoubleAttr_ObjVal);

            if(input.isFollowerNearOpt() == true){
                if(ys_eps_) delete[] ys_eps_;
                ys_eps_ = model_eval->get(GRB_DoubleAttr_X,y_eval,instance.getModel()->nb_follower_vars);
                fs_eps_ = f_eval.get(GRB_DoubleAttr_X);
            }else{
                if(ys_) delete[] ys_;
                ys_ = model_eval->get(GRB_DoubleAttr_X,y_eval,instance.getModel()->nb_follower_vars);
            }
        }
    }

    // Solve pessimistic problem.
    if(compute_wk == true){
        model_eval->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

        model_eval->update();
        model_eval->optimize();

        Status status_eval = statusFromGurobi(model_eval->get(GRB_IntAttr_Status));
        if(status_eval == Status::Optimal){
            Fw_ = model_eval->get(GRB_DoubleAttr_ObjVal);

            if(yw_) delete[] yw_;
            yw_ = model_eval->get(GRB_DoubleAttr_X,y_eval,instance.getModel()->nb_follower_vars);

            if(input.isFollowerNearOpt() == true) fw_eps_ = f_eval.get(GRB_DoubleAttr_X);
        }
    }

    // Solve interior (slack) problem
    if(compute_int == true){
        model_eval->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

        // Define slack variable and objective.
        GRBVar slack = model_eval->addVar(0.0,GRB_INFINITY,1.0,GRB_CONTINUOUS);

        // Remove leader objective.
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            BilevelVariable var = instance.getModel()->follower_vars[i];
            y_eval[i].set(GRB_DoubleAttr_Obj,0.0);
        }

        // Include slack vars at follower constraints.
        for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
            BilevelConstraint constr = instance.getModel()->follower_constrs[c];

            double norm = 0.0;
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                BilevelVariable var = instance.getModel()->follower_vars[i];
                norm += constr.coeffs[var.id]*constr.coeffs[var.id];
            }
            norm = std::sqrt(norm);

            if(constr.sense == '<') model_eval->chgCoeff(grb_constrs[c], slack, 1.0*norm);
            else if(constr.sense == '>') model_eval->chgCoeff(grb_constrs[c], slack, -1.0*norm);
        }

        // Include at near opt constraint.
        if(input.isFollowerNearOpt() == true) {

            double norm = 0.0;
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                BilevelVariable var = instance.getModel()->follower_vars[i];
                norm += var.obj_follower*var.obj_follower;
            }
            norm = std::sqrt(norm);

            model_eval->chgCoeff(grb_constrs[instance.getModel()->nb_follower_constrs], slack, 1.0*norm);
        }

        // Include lb and ub constraints with slack.
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            BilevelVariable var = instance.getModel()->follower_vars[i];
            
            if(!(var.lb <= -instance.getModel()->inf))
                model_eval->addConstr(y_eval[i] - slack >= var.lb);
            if(!(var.ub >=  instance.getModel()->inf)) 
                model_eval->addConstr(y_eval[i] + slack <= var.ub);
        }

        model_eval->update();
        model_eval->optimize();

        // model_eval->write("slack.lp");

        Status status_eval = statusFromGurobi(model_eval->get(GRB_IntAttr_Status));
        if(status_eval == Status::Optimal){
            std::cout << "slack max: " << model_eval->get(GRB_DoubleAttr_ObjVal) << std::endl;

            if(yi_) delete[] yi_;
            yi_ = model_eval->get(GRB_DoubleAttr_X,y_eval,instance.getModel()->nb_follower_vars);

            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                std::cout << yi_[i] << std::endl;
            }

            // if(input.isFollowerNearOpt() == true) fw_eps_ = f_eval.get(GRB_DoubleAttr_X);
        }


    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////