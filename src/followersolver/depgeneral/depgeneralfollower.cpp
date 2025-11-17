#include "depgeneralfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

DepGeneralFollowerSolver::DepGeneralFollowerSolver(const Input & input,const Instance & instance, LeaderSolver *leader, int pr) :
                                    AbstractFollowerSolver(input,instance,leader), pr(pr) { 
    // Initialize pointers to NULL.
    y.resize(instance.getNbScenarios(),NULL);
}

void DepGeneralFollowerSolver::defineLeaderObj(){
    // Define variables.
    defineMCMCVars();

    // Define constraints.
    defineMCMCConstrs();

    // Define objective.
    defineMCMCObj();
}

void DepGeneralFollowerSolver::defineMCMCVars(){
    // Define scenario variables representing follower primal ones.
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        definePrimalVars(y[s],"_" + std::to_string(s+1));
    }

    // Define step and step maximum and minimum value variables (alpha, alpha_min and alpha_max).
    alpha = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_CONTINUOUS);
    alpha_min = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_CONTINUOUS);
    alpha_max = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_CONTINUOUS);
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        alpha[s].set(GRB_DoubleAttr_LB,-instance.getStepSizeInterval());
        alpha[s].set(GRB_DoubleAttr_UB,instance.getStepSizeInterval());

        alpha_min[s].set(GRB_DoubleAttr_LB,-instance.getStepSizeInterval());
        alpha_min[s].set(GRB_DoubleAttr_UB,0.0);

        alpha_max[s].set(GRB_DoubleAttr_LB,0.0);
        alpha_max[s].set(GRB_DoubleAttr_UB,instance.getStepSizeInterval());
    }

    if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
        // Acceptance variables.
        q = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_BINARY);

        // Interval variables.
        w = leader->getGRBModel()->addVars(instance.getGeneralNbIntervals(),GRB_BINARY);
    }
}

void DepGeneralFollowerSolver::defineMCMCConstrs(){
    // Define alpha variables from alpha_min and alpha_max.
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        double tau = instance.getGeneralStep(pr,s);
        leader->getGRBModel()->addConstr(alpha[s] == (1.0 - tau)*alpha_min[s] + tau*alpha_max[s]);
    }

    // Define maximum and minimum step variables (alpha_min and alpha_max).
    if(input.useStepExplicitGeneral() == true) defineExplicitStep();
    else defineKKTStep();

    // Define sequence of points y.
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::Neutral){
        // All new points are accepted.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                if(s == 0) leader->getGRBModel()->addConstr(y[s][i] == ys[i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
                else leader->getGRBModel()->addConstr(y[s][i] == y[s-1][i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
            }
        }
    }
    else{
        // Define interval for probability functions according to follower optimal value.
        GRBVar f_prop = leader->getGRBModel()->addVar(0.0,1.0,0.0,GRB_CONTINUOUS);
        leader->getGRBModel()->addConstr(f_prop == (instance.getModel()->follower_ub - fs)/(instance.getModel()->follower_ub - instance.getModel()->follower_lb)); 
        for(int k = 0; k < instance.getGeneralNbIntervals(); ++k){
            leader->getGRBModel()->addConstr(f_prop >= (instance.getGeneralIntervalValue(k) + input.getEpsBigM())*w[k]);
            leader->getGRBModel()->addConstr(f_prop <= 1.0 + (instance.getGeneralIntervalValue(k+1) - 1.0)*w[k]); 

            leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, f_prop >= instance.getGeneralIntervalValue(k) + input.getEpsBigM());
            leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, f_prop <= instance.getGeneralIntervalValue(k+1));
        }
        
        // Define acceptance constraints.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                if(s == 0){
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == ys[i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == ys[i]);
                }else{
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == y[s-1][i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == y[s-1][i]);
                }
            }
        }

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional){
            
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                double r = instance.getGeneralAccep(pr,s);

                GRBLinExpr expr = 0;
                for(int k = 0; k < instance.getGeneralNbIntervals(); ++k)
                    expr += (((r - 1.0)*instance.getGeneralRefProbab())/(input.getScalingGeneral()*instance.getGeneralIntervalCoeff(k)))*w[k];

                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr -= (1.0 - r)*var.obj_follower*ys[i];    
                    else expr -= (1.0 - r)*var.obj_follower*y[s-1][i];
                }
                expr -= instance.getGeneralCoeffAlphaObj(pr,s)*alpha[s];
                expr -= instance.getRefLeaderObj()*(r - 1.0);

                double bigm_ub = r*instance.getGeneralMaxProbab() - instance.getGeneralMinProbab();
                double bigm_lb = r*instance.getGeneralMinProbab() - instance.getGeneralMaxProbab();

                leader->getGRBModel()->addConstr(expr <= bigm_ub*(1.0 - q[s]));
                leader->getGRBModel()->addConstr(expr >= input.getEpsBigM() + bigm_lb*q[s]);
                
                leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
                leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.0 + input.getEpsBigM());
            }
        }

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::StrongFragile){
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                double r = instance.getGeneralAccep(pr,s);

                GRBLinExpr expr = 0;
                for(int k = 0; k < instance.getGeneralNbIntervals(); ++k)
                    expr += ((std::log(r))/(input.getScalingGeneral()*instance.getGeneralIntervalCoeff(k)))*w[k];
                expr -= instance.getGeneralCoeffAlphaObj(pr,s)*alpha[s];

                double bigm_ub = std::log(r) - std::log(instance.getGeneralMinProbab()) + std::log(instance.getGeneralMaxProbab());
                double bigm_lb = std::log(r) - std::log(instance.getGeneralMaxProbab()) + std::log(instance.getGeneralMinProbab());

                leader->getGRBModel()->addConstr(expr <= bigm_ub*(1.0 - q[s]));
                leader->getGRBModel()->addConstr(expr >= input.getEpsBigM() + bigm_lb*q[s]);
                
                leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
                leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.0 + input.getEpsBigM());
            }
        }
    }  
}

void DepGeneralFollowerSolver::defineMCMCObj(){
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            y[s][i].set(GRB_DoubleAttr_Obj,(instance.getModel()->follower_vars[i].obj_leader/(double)instance.getNbScenarios()));
        }
    }
}

void DepGeneralFollowerSolver::defineExplicitStep(){
    // Define constraints for alpha_min.
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        
        // Auxiliar binary variables.
        int nb_constrs_min = instance.getGeneralNbConstrsAlpha(pr,s,0) + instance.getGeneralNbBdlbsAlpha(pr,s,0) + instance.getGeneralNbBdubsAlpha(pr,s,0) + 1;
        GRBVar * b_alpha_min =  leader->getGRBModel()->addVars(nb_constrs_min,GRB_BINARY);

        // Only one possible lower bound can be equal to one.
        GRBLinExpr sum_expr = 0;
        for(int i = 0; i < nb_constrs_min; ++i) sum_expr += b_alpha_min[i];
        leader->getGRBModel()->addConstr(sum_expr == 1);

        // Constraints impacting alpha_min.
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,0,j);
            
            // Usual follower constraints.
            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];

                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                    BilevelVariable var = instance.getModel()->leader_vars[i];
                    expr += constr.coeffs[var.id]*leader->getX(i); 
                }

                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += constr.coeffs[var.id]*ys[i];    
                    else expr += constr.coeffs[var.id]*y[s-1][i];
                }

                if(constr.sense == '<'){
                    expr += instance.getGeneralCoeffAlpha(pr,s,c)*alpha_min[s];
                    leader->getGRBModel()->addConstr(expr <= constr.rhs);
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[j], 1, expr >= constr.rhs);
                }
                if(constr.sense == '>'){
                    expr += instance.getGeneralCoeffAlpha(pr,s,c)*alpha_min[s];
                    leader->getGRBModel()->addConstr(expr >= constr.rhs);
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[j], 1, expr <= constr.rhs);
                }
            }
            // Near optimality constraint.
            else{
                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += var.obj_follower*ys[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObj(pr,s)*alpha_min[s];
                leader->getGRBModel()->addConstr(expr <= (1.0 + input.getEpsNearOpt())*fs);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[j], 1, expr >= (1.0 + input.getEpsNearOpt())*fs);
            }
        }
        
        // Lower bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0){
                leader->getGRBModel()->addConstr(ys[i] + dir*alpha_min[s] >= var.lb);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+j], 1, ys[i] + dir*alpha_min[s] <= var.lb);
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*alpha_min[s] >= var.lb);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+j], 1, y[s-1][i] + dir*alpha_min[s] <= var.lb);
            }  
        }

        // Upper bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0){
                leader->getGRBModel()->addConstr(ys[i] + dir*alpha_min[s] <= var.ub);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,0)+j], 1, ys[i] + dir*alpha_min[s] >= var.ub);
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*alpha_min[s] <= var.ub);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,0)+j], 1, y[s-1][i] + dir*alpha_min[s] >= var.ub);
            }
        }

        // Lower bound on alpha_min.
        leader->getGRBModel()->addConstr(alpha_min[s] >= -instance.getStepSizeInterval());
        leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[nb_constrs_min-1], 1, alpha_min[s] <= -instance.getStepSizeInterval());
    }

    // Define constraints for alpha_max.
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        
        // Auxiliar binary variables.
        int nb_constrs_max = instance.getGeneralNbConstrsAlpha(pr,s,1) + instance.getGeneralNbBdlbsAlpha(pr,s,1) + instance.getGeneralNbBdubsAlpha(pr,s,1) + 1;
        GRBVar * b_alpha_max =  leader->getGRBModel()->addVars(nb_constrs_max,GRB_BINARY);

        // Only one possible upper bound can be equal to one.
        GRBLinExpr sum_expr = 0;
        for(int i = 0; i < nb_constrs_max; ++i) sum_expr += b_alpha_max[i];
        leader->getGRBModel()->addConstr(sum_expr == 1);

        // Constraints impacting alpha_max.
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,1); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,1,j);
            
            // Usual follower constraints.
            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];

                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                    BilevelVariable var = instance.getModel()->leader_vars[i];
                    expr += constr.coeffs[var.id]*leader->getX(i); 
                }

                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += constr.coeffs[var.id]*ys[i];    
                    else expr += constr.coeffs[var.id]*y[s-1][i];
                }

                if(constr.sense == '<'){
                    expr += instance.getGeneralCoeffAlpha(pr,s,c)*alpha_max[s];
                    leader->getGRBModel()->addConstr(expr <= constr.rhs);
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[j], 1, expr >= constr.rhs);
                }
                if(constr.sense == '>'){
                    expr += instance.getGeneralCoeffAlpha(pr,s,c)*alpha_max[s];
                    leader->getGRBModel()->addConstr(expr >= constr.rhs);
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[j], 1, expr <= constr.rhs);
                }
            }
            // Near optimality constraint.
            else{
                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += var.obj_follower*ys[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObj(pr,s)*alpha_max[s];
                leader->getGRBModel()->addConstr(expr <= (1.0 + input.getEpsNearOpt())*fs);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[j], 1, expr >= (1.0 + input.getEpsNearOpt())*fs);
            }
        }

        // Lower bound constraints defining alpha_max.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0){
                leader->getGRBModel()->addConstr(ys[i] + dir*alpha_max[s] >= var.lb);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+j], 1, ys[i] + dir*alpha_max[s] <= var.lb);
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*alpha_max[s] >= var.lb);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+j], 1, y[s-1][i] + dir*alpha_max[s] <= var.lb);
            }  
        }

        // Upper bound constraints defining alpha_max.
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0){
                leader->getGRBModel()->addConstr(ys[i] + dir*alpha_max[s] <= var.ub);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+instance.getGeneralNbBdlbsAlpha(pr,s,1)+j], 1, ys[i] + dir*alpha_max[s] >= var.ub);
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*alpha_max[s] <= var.ub);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+instance.getGeneralNbBdlbsAlpha(pr,s,1)+j], 1, y[s-1][i] + dir*alpha_max[s] >= var.ub);
            }
        }

        // Upper bound on alpha_max.
        leader->getGRBModel()->addConstr(alpha_max[s] <= instance.getStepSizeInterval());
        leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[nb_constrs_max-1], 1, alpha_max[s] >= instance.getStepSizeInterval());
    }  
}

void DepGeneralFollowerSolver::defineKKTStep(){
    // Define primal constraints and their slack variables.
    std::vector<GRBVar*> slacks;
    std::vector<GRBVar*> slack_lbs;
    std::vector<GRBVar*> slack_ubs;
    std::vector<GRBVar*> slack_min;
    std::vector<GRBVar*> slack_max;
    defineAlphaPrimalConstrs(slacks,slack_lbs,slack_ubs,slack_min,slack_max);

    // Define dual variables.
    std::vector<GRBVar*> duals;
    std::vector<GRBVar*> duals_lbs;
    std::vector<GRBVar*> duals_ubs;
    std::vector<GRBVar*> duals_min;
    std::vector<GRBVar*> duals_max;
    defineAlphaDualVars(duals,duals_lbs,duals_ubs,duals_min,duals_max);
    
    // Define dual constraints.
    defineAlphaDualConstrs(duals,duals_lbs,duals_ubs,duals_min,duals_max);

    // Define complementarity constraints.
    defineAlphaCompConstrs(slacks,slack_lbs,slack_ubs,slack_min,slack_max,
                           duals,duals_lbs,duals_ubs,duals_min,duals_max);
}

void DepGeneralFollowerSolver::defineAlphaPrimalConstrs(std::vector<GRBVar*> & slacks, 
                                        std::vector<GRBVar*> & slack_lbs, std::vector<GRBVar*> & slack_ubs, 
                                        std::vector<GRBVar*> & slack_min, std::vector<GRBVar*> & slack_max){
    slacks.resize(instance.getNbScenarios());
    slack_lbs.resize(instance.getNbScenarios());
    slack_ubs.resize(instance.getNbScenarios());
    slack_min.resize(instance.getNbScenarios());
    slack_max.resize(instance.getNbScenarios());

    for(int s = 0; s < instance.getNbScenarios(); ++s){
        slacks[s] = leader->getGRBModel()->addVars(instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbConstrsAlpha(pr,s,1),GRB_CONTINUOUS);
        slack_lbs[s] = leader->getGRBModel()->addVars(instance.getGeneralNbBdlbsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,1),GRB_CONTINUOUS);
        slack_ubs[s] = leader->getGRBModel()->addVars(instance.getGeneralNbBdubsAlpha(pr,s,0)+instance.getGeneralNbBdubsAlpha(pr,s,1),GRB_CONTINUOUS);
        slack_min[s] = leader->getGRBModel()->addVars(2,GRB_CONTINUOUS);
        slack_max[s] = leader->getGRBModel()->addVars(2,GRB_CONTINUOUS);
    }

    // Define primal constraints and their slack variables.
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        // alpha_min
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,0,j);

            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];

                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                    BilevelVariable var = instance.getModel()->leader_vars[i];
                    expr += constr.coeffs[var.id]*leader->getX(i); 
                }

                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += constr.coeffs[var.id]*ys[i];    
                    else expr += constr.coeffs[var.id]*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlpha(pr,s,c)*alpha_min[s];
                if(constr.sense == '<') leader->getGRBModel()->addConstr(expr + slacks[s][j] == constr.rhs);
                if(constr.sense == '>') leader->getGRBModel()->addConstr(expr - slacks[s][j] == constr.rhs);
            }
            else{
                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += var.obj_follower*ys[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObj(pr,s)*alpha_min[s];
                leader->getGRBModel()->addConstr(expr + slacks[s][j] == (1.0 + input.getEpsNearOpt())*fs);
            }
        }

        // alpha_max
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,1); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,1,j);

            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];

                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                    BilevelVariable var = instance.getModel()->leader_vars[i];
                    expr += constr.coeffs[var.id]*leader->getX(i); 
                }

                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += constr.coeffs[var.id]*ys[i];    
                    else expr += constr.coeffs[var.id]*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlpha(pr,s,c)*alpha_max[s];
                if(constr.sense == '<') leader->getGRBModel()->addConstr(expr + slacks[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j] == constr.rhs);
                if(constr.sense == '>') leader->getGRBModel()->addConstr(expr - slacks[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j] == constr.rhs);
            }
            else{
                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += var.obj_follower*ys[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObj(pr,s)*alpha_max[s];
                leader->getGRBModel()->addConstr(expr + slacks[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j] == (1.0 + input.getEpsNearOpt())*fs);
            }
        }
    }

    // Define primal bound constraints on y and their slacks.
    for(int s = 0; s < instance.getNbScenarios(); ++s){

        // alpha_min
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,0,j);

            BilevelVariable var = instance.getModel()->follower_vars[i];
            double dir = instance.getGeneralDirection(pr,s,i);
            leader->getGRBModel()->addConstr(y[s][i] + dir*alpha_min[s] - slack_lbs[s][j] == var.lb);
        }

        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,0,j);

            BilevelVariable var = instance.getModel()->follower_vars[i];
            double dir = instance.getGeneralDirection(pr,s,i);
            leader->getGRBModel()->addConstr(y[s][i] + dir*alpha_min[s] + slack_ubs[s][j] == var.ub);
        }

        // alpha_max
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            double dir = instance.getGeneralDirection(pr,s,i);
            leader->getGRBModel()->addConstr(y[s][i] + dir*alpha_max[s] - slack_lbs[s][instance.getGeneralNbBdlbsAlpha(pr,s,0)+j] == var.lb);
        }

        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            double dir = instance.getGeneralDirection(pr,s,i);
            leader->getGRBModel()->addConstr(y[s][i] + dir*alpha_max[s] + slack_ubs[s][instance.getGeneralNbBdubsAlpha(pr,s,0)+j] == var.ub);
        }

    }

    // Define primal constraints for alpha bounds.
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        leader->getGRBModel()->addConstr(alpha_min[s] + slack_min[s][1] == 0.0);
        leader->getGRBModel()->addConstr(alpha_min[s] - slack_min[s][0] == -instance.getStepSizeInterval());
        
        leader->getGRBModel()->addConstr(alpha_max[s] + slack_max[s][1] == instance.getStepSizeInterval());
        leader->getGRBModel()->addConstr(alpha_max[s] - slack_max[s][0] == 0.0);
    }
}

void DepGeneralFollowerSolver::defineAlphaDualVars(std::vector<GRBVar*> & duals, 
                                        std::vector<GRBVar*> & duals_lbs, std::vector<GRBVar*> & duals_ubs, 
                                        std::vector<GRBVar*> & duals_min, std::vector<GRBVar*> & duals_max){                               
    // Defining all dual constraints as >= 0 variables. There are no equalities in these problems.
    duals.resize(instance.getNbScenarios());
    duals_lbs.resize(instance.getNbScenarios());
    duals_ubs.resize(instance.getNbScenarios());
    duals_min.resize(instance.getNbScenarios());
    duals_max.resize(instance.getNbScenarios());

    for(int s = 0; s < instance.getNbScenarios(); ++s){
        duals[s] = leader->getGRBModel()->addVars(instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbConstrsAlpha(pr,s,1),GRB_CONTINUOUS);
        duals_lbs[s] = leader->getGRBModel()->addVars(instance.getGeneralNbBdlbsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,1),GRB_CONTINUOUS);
        duals_ubs[s] = leader->getGRBModel()->addVars(instance.getGeneralNbBdubsAlpha(pr,s,0)+instance.getGeneralNbBdubsAlpha(pr,s,1),GRB_CONTINUOUS);
        duals_min[s] = leader->getGRBModel()->addVars(2,GRB_CONTINUOUS);
        duals_max[s] = leader->getGRBModel()->addVars(2,GRB_CONTINUOUS);
    }
}

void DepGeneralFollowerSolver::defineAlphaDualConstrs(std::vector<GRBVar*> & duals, 
                                        std::vector<GRBVar*> & duals_lbs, std::vector<GRBVar*> & duals_ubs, 
                                        std::vector<GRBVar*> & duals_min, std::vector<GRBVar*> & duals_max){

    // Dual constraint for alpha_min. (Primal is a minimization problem).
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        GRBLinExpr expr = 0;
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,0,j);

            BilevelConstraint constr = instance.getModel()->follower_constrs[c];
            if(constr.sense == '<') expr += -duals[s][j]*instance.getGeneralCoeffAlpha(pr,s,c);
            else if(constr.sense == '>') expr += duals[s][j]*instance.getGeneralCoeffAlpha(pr,s,c);
        }
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j) {
            int i = instance.getGeneralBdlbsAlpha(pr,s,0,j);
            expr += duals_lbs[s][j]*instance.getGeneralDirection(pr,s,i);
        }
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j) {
            int i = instance.getGeneralBdubsAlpha(pr,s,0,j);
            expr += -duals_ubs[s][j]*instance.getGeneralDirection(pr,s,i);
        } 
        expr += duals_min[s][0] - duals_min[s][1];
        leader->getGRBModel()->addConstr(expr == 1.0);
    }
    
    // Dual constraint for alpha_max. (Primal is a maximization problem).
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        GRBLinExpr expr = 0;
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,1); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,1,j);
            BilevelConstraint constr = instance.getModel()->follower_constrs[c];
            if(constr.sense == '<') expr += duals[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j]*instance.getGeneralCoeffAlpha(pr,s,c);
            else if(constr.sense == '>') expr += -duals[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j]*instance.getGeneralCoeffAlpha(pr,s,c);
        }
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j) {
            int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
            expr += -duals_lbs[s][instance.getGeneralNbBdlbsAlpha(pr,s,0)+j]*instance.getGeneralDirection(pr,s,i);
        }
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j) {
            int i = instance.getGeneralBdubsAlpha(pr,s,1,j) ;
            expr += duals_lbs[s][instance.getGeneralNbBdubsAlpha(pr,s,0)+j]*instance.getGeneralDirection(pr,s,i);
        }
        expr += -duals_max[s][0] + duals_max[s][1];
        leader->getGRBModel()->addConstr(expr == 1.0);
    }
}

void DepGeneralFollowerSolver::defineAlphaCompConstrs(std::vector<GRBVar*> & slacks, 
                                        std::vector<GRBVar*> & slack_lbs, std::vector<GRBVar*> & slack_ubs, 
                                        std::vector<GRBVar*> & slack_min, std::vector<GRBVar*> & slack_max,
                                        std::vector<GRBVar*> & duals, 
                                        std::vector<GRBVar*> & duals_lbs, std::vector<GRBVar*> & duals_ubs, 
                                        std::vector<GRBVar*> & duals_min, std::vector<GRBVar*> & duals_max){
    
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        // alpha_min
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
            GRBVar v[2] = {slacks[s][j], duals[s][j]};
            double w[2] = {1.0, 2.0};
            leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);
        }
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j){
            GRBVar vlb[2] = {slack_lbs[s][j], duals_lbs[s][j]};
            double wlb[2] = {1.0, 2.0};
            leader->getGRBModel()->addSOS(vlb, wlb, 2, GRB_SOS_TYPE1);
        }
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j){
            GRBVar vub[2] = {slack_ubs[s][j], duals_ubs[s][j]};
            double wub[2] = {1.0, 2.0};
            leader->getGRBModel()->addSOS(vub, wub, 2, GRB_SOS_TYPE1);
        }
        GRBVar vminlb[2] = {slack_min[s][0],duals_min[s][0]};
        double wminlb[2] = {1.0, 2.0};
        leader->getGRBModel()->addSOS(vminlb, wminlb, 2, GRB_SOS_TYPE1);

        GRBVar vminub[2] = {slack_min[s][1],duals_min[s][1]};
        double wminub[2] = {1.0, 2.0};
        leader->getGRBModel()->addSOS(vminub, wminub, 2, GRB_SOS_TYPE1);

        // alpha_max
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,1); ++j){
            GRBVar v[2] = {slacks[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j], duals[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j]};
            double w[2] = {1.0, 2.0};
            leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);
        }
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j){
            GRBVar vlb[2] = {slack_lbs[s][instance.getGeneralNbBdlbsAlpha(pr,s,0)+j], duals_lbs[s][instance.getGeneralNbBdlbsAlpha(pr,s,0)+j]};
            double wlb[2] = {1.0, 2.0};
            leader->getGRBModel()->addSOS(vlb, wlb, 2, GRB_SOS_TYPE1);
        }
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j){
            GRBVar vub[2] = {slack_ubs[s][instance.getGeneralNbBdubsAlpha(pr,s,0)+j], duals_ubs[s][instance.getGeneralNbBdubsAlpha(pr,s,0)+j]};
            double wub[2] = {1.0, 2.0};
            leader->getGRBModel()->addSOS(vub, wub, 2, GRB_SOS_TYPE1);
        }
        GRBVar vmaxlb[2] = {slack_max[s][0],duals_max[s][0]};
        double wmaxlb[2] = {1.0, 2.0};
        leader->getGRBModel()->addSOS(vmaxlb, wmaxlb, 2, GRB_SOS_TYPE1);

        GRBVar vmaxub[2] = {slack_max[s][1],duals_max[s][1]};
        double wmaxub[2] = {1.0, 2.0};
        leader->getGRBModel()->addSOS(vmaxub, wmaxub, 2, GRB_SOS_TYPE1);
    }
}

double DepGeneralFollowerSolver::evaluate(){
    // Initialize y_eval and y_test.
    std::vector<double> y_eval(ys_, ys_ + instance.getModel()->nb_follower_vars);
    std::vector<double> y_test(ys_, ys_ + instance.getModel()->nb_follower_vars);

    // Initialize evaluation values.
    eval_avg = 0.0;
    eval_sce.resize(instance.getNbValidateScenarios(),0.0);

    for(int s = 0; s < instance.getNbValidateScenarios(); ++s){
        // Computing alpha_min.
        double eval_alpha_min = instance.defineEvalMinStep(s,leader->getX_(),y_eval,fs_);
        
        // Computing alpha_max.
        double eval_alpha_max = instance.defineEvalMaxStep(s,leader->getX_(),y_eval,fs_);
        
        // Computing alpha.
        double tau = instance.getGeneralStepEval(s);
        double eval_alpha = (1.0-tau)*eval_alpha_min + tau*eval_alpha_max;

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::Neutral){
            // New point always accepted.
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i)
                y_eval[i] = y_eval[i] + eval_alpha*instance.getGeneralDirectionEval(s,i);
        }
        else {
            // Compute possible new point.
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i)
                y_test[i] = y_eval[i] + eval_alpha*instance.getGeneralDirectionEval(s,i);

            // Compute probability to accept new point.
            double h = instance.getDepGeneralProbab(fs_,leader->getX_(),y_test)/instance.getDepGeneralProbab(fs_,leader->getX_(),y_eval);

            // Accept new point under this condition.
            if(instance.getGeneralAccepEval(s) <= h){
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i)
                    y_eval[i] = y_test[i];
            }
        }

        // Compute leader objective from this scenario.
        eval_sce[s] = 0.0;
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i)
            eval_sce[s] += instance.getModel()->follower_vars[i].obj_leader*y_eval[i];
        eval_avg += eval_sce[s]/(double)instance.getNbValidateScenarios();
    }

    double eval = eval_avg;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        eval += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);

    return eval;
}

void DepGeneralFollowerSolver::computeStrongWeakSolutions(){
    // Compute both strong and weak solutions.
    AbstractFollowerSolver::computeStrongWeakSolutions(true,true);
}

