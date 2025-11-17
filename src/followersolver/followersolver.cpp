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

void AbstractFollowerSolver::definePrimalConstrs(GRBVar * & y, GRBVar * & slacks, GRBVar * & slack_lbs, GRBVar * & slack_ubs, std::string v){
    // Define bound slack variables.
    slack_lbs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    slack_ubs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];

        if(!(var.lb <= -instance.getModel()->inf))
            leader->getGRBModel()->addConstr(slack_lbs[i] == y[i] - var.lb);
        if(!(var.ub >=  instance.getModel()->inf))
            leader->getGRBModel()->addConstr(slack_ubs[i] == var.ub - y[i]);
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
            leader->getGRBModel()->addConstr(expr + slacks[c] == constr.rhs, constr.name + v);
        }else if(constr.sense == '>'){
            leader->getGRBModel()->addConstr(expr - slacks[c] == constr.rhs, constr.name + v);
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

void AbstractFollowerSolver::defineDualVars(GRBVar * & ds, GRBVar * & dlbs, GRBVar * & dubs, bool min_problem){
    // Dual variables.
    int nb_constrs = instance.getModel()->nb_follower_constrs;
    if(min_problem == false) ++nb_constrs; // Optimality constraint
    
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

    // Dual variables for bound constraints.
    dlbs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    dubs = leader->getGRBModel()->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
}

void AbstractFollowerSolver::defineDualConstrs(GRBVar * & ds, GRBVar * & dlbs, GRBVar * & dubs, bool min_problem){
    // Dual constraints.
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];

        // Consider coefficient all primal constraints.
        GRBLinExpr expr = 0;
        for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
            BilevelConstraint constr = instance.getModel()->follower_constrs[c];
            expr += constr.coeffs[var.id]*ds[c];
        }

        // Pessimistic case we should also consider optimality constraint.
        if(min_problem == false) expr += var.obj_follower*ds[instance.getModel()->nb_follower_constrs];

        if(var.lb <= -instance.getModel()->inf && var.ub >= instance.getModel()->inf) 
            leader->getGRBModel()->addConstr(expr == var.obj_follower);
        else if(var.lb <= -instance.getModel()->inf){ 
            if(min_problem == true) 
                leader->getGRBModel()->addConstr(expr - dubs[i] == var.obj_follower);
            else leader->getGRBModel()->addConstr(expr + dubs[i] == var.obj_follower);
        }else if(var.ub >= instance.getModel()->inf){
            if(min_problem == true) 
                leader->getGRBModel()->addConstr(expr + dlbs[i] == var.obj_follower);
            else leader->getGRBModel()->addConstr(expr - dlbs[i] == var.obj_follower);
        }else{
            if(min_problem == true) 
                leader->getGRBModel()->addConstr(expr + dlbs[i] - dubs[i] == var.obj_follower);
            else leader->getGRBModel()->addConstr(expr - dlbs[i] + dubs[i] == var.obj_follower);
        } 
    }
}

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

void AbstractFollowerSolver::defineFollowerObj(GRBVar * & y, std::string v){
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
    leader->getGRBModel()->addConstr(fs == expr,"Follower_Optimal_Objective" + v);
}

void AbstractFollowerSolver::defineNearOptimalObj(GRBVar & f_eps, GRBVar * & y_eps, GRBVar & slack_eps, std::string v){
    slack_eps = leader->getGRBModel()->addVar(0.0,GRB_INFINITY,0.0,GRB_CONTINUOUS);
    f_eps = leader->getGRBModel()->addVar(-GRB_INFINITY,GRB_INFINITY,0.0,GRB_CONTINUOUS,"f_eps" + v);
    GRBLinExpr expr = 0;
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        expr += var.obj_follower*y_eps[i];
    }
    leader->getGRBModel()->addConstr(f_eps == expr, "Follower_Eps_Objective" + v);
    leader->getGRBModel()->addConstr(f_eps + slack_eps == (1.0 + input.getEpsNearOpt())*fs ,"Constr_Follower_Eps_Objective" + v);
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

void AbstractFollowerSolver::defineOptFollower(){
    // Define primal variables.
    definePrimalVars(ys,"_s");
    
    // Define primal constraints and slacks.
    GRBVar * slacks = NULL;    
    GRBVar * slack_lbs = NULL;
    GRBVar * slack_ubs = NULL;
    definePrimalConstrs(ys,slacks,slack_lbs,slack_ubs,"_s");
    
    // Define dual variables.
    GRBVar * ds = NULL;
    GRBVar * dlbs = NULL;
    GRBVar * dubs = NULL;
    defineDualVars(ds,dlbs,dubs,true);

    // Define dual constraints.
    defineDualConstrs(ds,dlbs,dubs,true);

    // Define complementarity constraints.
    defineCompConstrs(slacks,slack_lbs,slack_ubs,ds,dlbs,dubs);

    // Define follower optimal value.
    defineFollowerObj(ys,"_s");

    // Define auxiliar primal constraints in the case of near optimality.
    if(input.isFollowerNearOpt() == true){
        // Define primal eps near optimal variables.
        definePrimalVars(ys_eps,"_s_eps");

        // Define primal constraints.
        definePrimalConstrs(ys_eps,"_s_eps");

        // Define near optimality constraint.
        GRBVar slack_eps;
        defineNearOptimalObj(fs_eps,ys_eps,slack_eps,"_s_eps");
    }

    // Define leader optimal value.
    if(input.isFollowerNearOpt() == true) 
        defineLeaderObj(Fs,ys_eps,"_opt_eps");
    else defineLeaderObj(Fs,ys,"_opt");
}

void AbstractFollowerSolver::definePesFollower(){
    // Define primal variables.
    definePrimalVars(yw,"_w");

    // Define primal constraints and slacks.
    GRBVar * slacks =  NULL;
    GRBVar * slack_lbs = NULL;
    GRBVar * slack_ubs = NULL;
    definePrimalConstrs(yw,slacks,slack_lbs,slack_ubs,"_w");

    // Define primal constraint and slack corresponding to optimality, specifically for pessimistic problem.
    // Define it according to optimality or near optimality.
    GRBVar slack_eps;
    if(input.isFollowerNearOpt() == true)
        defineNearOptimalObj(fw_eps,yw,slack_eps,"_w_eps");
    else defineFollowerObj(yw,"_w");

    // Define dual variables (primal is a maximization problem)
    GRBVar * ds = NULL;
    GRBVar * dlbs = NULL;
    GRBVar * dubs = NULL;
    defineDualVars(ds,dlbs,dubs,false);

    // Define dual constraints.
    defineDualConstrs(ds,dlbs,dubs,false);
    
    // Complementarity constraints
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

void AbstractFollowerSolver::upd_solution(){
    if(leader->getNbSol() > 0){
        if(ys_) delete[] ys_;
        if(yw_) delete[] yw_;
        if(ys_eps_) delete[] ys_eps_;

        if(ys) ys_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,ys,instance.getModel()->nb_follower_vars);
        if(yw) yw_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,yw,instance.getModel()->nb_follower_vars);
        if(ys_eps) ys_eps_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,ys_eps,instance.getModel()->nb_follower_vars);

        Fs_ = Fs.get(GRB_DoubleAttr_X);
        Fw_ = Fw.get(GRB_DoubleAttr_X);

        fs_ = fs.get(GRB_DoubleAttr_X);
        if(input.isFollowerNearOpt() == true) {
            fs_eps_ = fs_eps.get(GRB_DoubleAttr_X);
            fw_eps_ = fw_eps.get(GRB_DoubleAttr_X);
        }
    }
}

void AbstractFollowerSolver::computeStrongWeakSolutions(bool compute_str, bool compute_wk){
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
            model_eval->addConstr(expr == constr.rhs, constr.name);
        else if(constr.sense == '<'){
            model_eval->addConstr(expr <= constr.rhs, constr.name);
        }else if(constr.sense == '>'){
            model_eval->addConstr(expr >= constr.rhs, constr.name);
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
    if(input.isFollowerNearOpt() == true)   
        model_eval->addConstr(expr <= (1.0 + input.getEpsNearOpt())*fs_,"Follower_Objective");
    else
        model_eval->addConstr(expr <= fs_,"Follower_Objective");

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
}
