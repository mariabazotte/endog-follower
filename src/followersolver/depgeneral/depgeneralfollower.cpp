#include "depgeneralfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

DepGeneralFollowerSolver::DepGeneralFollowerSolver(const Input & input, Instance & instance, LeaderSolver *leader, int pr) :
                                    AbstractFollowerSolver(input,instance,leader), pr(pr) { 
    // Initialize pointers to NULL.
    y.resize(instance.getNbScenarios(),NULL);
    
    // Create problem
    create();
}

void DepGeneralFollowerSolver::defineLeaderObj(){
    // Define variables.
    defineMCMCVars();

    // Define constraints.
    defineMCMCConstrs();

    // Define objective.
    defineMCMCObj();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        w = leader->getGRBModel()->addVars(input.getNbIntervalsGeneral(),GRB_BINARY);

        // Auxiliar variables for bilinear term w_k*(Fw-Fs).
        w_delta = leader->getGRBModel()->addVars(input.getNbIntervalsGeneral(),GRB_CONTINUOUS);
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

    // Define sequence of points y. Start sequence with chebyshev center, defined by yi.
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::Neutral){
        // All new points are accepted.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                if(s == 0) leader->getGRBModel()->addConstr(y[s][i] == yi[i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
                else leader->getGRBModel()->addConstr(y[s][i] == y[s-1][i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
            }
        }
    }
    else{
        bool use_bigm = false;
        // Define interval for probability distribution according to follower optimal value.
        GRBVar fsp = leader->getGRBModel()->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
        leader->getGRBModel()->addConstr(fsp == (instance.getModel()->follower_ub - fs) / (instance.getModel()->follower_ub - instance.getModel()->follower_lb)); 
        
        GRBLinExpr sum_expr = 0;
        for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
            sum_expr += w[k];

            leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, fsp >= input.getIntValueGeneral(k) + input.getEpsBigM());
            leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, fsp <= input.getIntValueGeneral(k+1));

            leader->getGRBModel()->addConstr(fsp >= (input.getIntValueGeneral(k) + input.getEpsBigM())*w[k]);
            leader->getGRBModel()->addConstr(fsp <= 1.0 + (input.getIntValueGeneral(k+1) - 1.0)*w[k]);   
        }
        leader->getGRBModel()->addConstr(sum_expr == 1.0);
        
        // Define w_delta_k as multiplication of w_k*(Fw - Fs).
        for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
            leader->getGRBModel()->addGenConstrIndicator(w[k], 0.0, w_delta[k] == 0.0);
            leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, w_delta[k] == (Fw -Fs));

            if(use_bigm == true){
                // Large bigM values.
                leader->getGRBModel()->addConstr(w_delta[k] <= (Fw -Fs));
                leader->getGRBModel()->addConstr(w_delta[k] >= (Fw -Fs) - (instance.getModel()->leader_ub - instance.getModel()->leader_lb)*(1.0 - w[k]));
                leader->getGRBModel()->addConstr(w_delta[k] <= (instance.getModel()->leader_ub - instance.getModel()->leader_lb)*w[k]);
            }
        }

        // Define constraints for sequence of points: new points are accepted according to probability.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                if(s == 0){
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == yi[i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == yi[i]);
                }else{
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == y[s-1][i] + alpha[s]*instance.getGeneralDirection(pr,s,i));
                    leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == y[s-1][i]);
                }
            }
        }

        // Definition of acceptance variables.
        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional){
            // Defining acceptance variables in proportional case.
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                double r = instance.getGeneralAccep(pr,s);

                GRBLinExpr expr = 0;
                for(int k = 0; k < input.getNbIntervalsGeneral(); ++k)
                    expr += (((r - 1.0) * input.getMaxIntCoeffGeneral()) / (2.0 * input.getIntCoeffGeneral(k)))*w_delta[k];

                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    
                    if(s == 0) expr += (r - 1.0)*var.obj_leader*yi[i];    
                    else expr += (r - 1.0)*var.obj_leader*y[s-1][i];
                }
                expr -= (r - 1.0) * (Fw + Fs) / 2.0;
                expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*alpha[s];

                leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
                leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.0 + input.getEpsBigM());
                
                if(use_bigm == true){
                    // Large bigM values.
                    double bigm_ub = r*instance.getGeneralMaxProbab() - instance.getGeneralMinProbab();
                    double bigm_lb = r*instance.getGeneralMinProbab() - instance.getGeneralMaxProbab();

                    leader->getGRBModel()->addConstr(expr <= bigm_ub*(1.0 - q[s]));
                    leader->getGRBModel()->addConstr(expr >= input.getEpsBigM() + bigm_lb*q[s]); 
                } 
            }
        }

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragile){
            // Defining acceptance variables in strong-fragile case.
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                double r = instance.getGeneralAccep(pr,s);

                GRBLinExpr expr = 0;
                for(int k = 0; k < input.getNbIntervalsGeneral(); ++k)
                    expr += (std::log(r) / input.getIntCoeffGeneral(k))*w_delta[k];

                expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*alpha[s];

                leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
                leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.0 + input.getEpsBigM());

                if(use_bigm == true){
                    // Large bigM values.
                    double bigm_ub = std::log(r) - std::log(instance.getGeneralMinProbab()) + std::log(instance.getGeneralMaxProbab());
                    double bigm_lb = std::log(r) - std::log(instance.getGeneralMaxProbab()) + std::log(instance.getGeneralMinProbab());

                    leader->getGRBModel()->addConstr(expr <= bigm_ub*(1.0 - q[s]));
                    leader->getGRBModel()->addConstr(expr >= input.getEpsBigM() + bigm_lb*q[s]); 
                }
            }
        }
    }  
}

void DepGeneralFollowerSolver::defineMCMCObj(){
    double nb_scenarios = instance.getNbScenarios()/instance.getNbThinning();
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            if(instance.chosenScenario(s) == true){
                y[s][i].set(GRB_DoubleAttr_Obj,var.obj_leader/nb_scenarios);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
                    if(s == 0) expr += constr.coeffs[var.id]*yi[i];    
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
                    if(s == 0) expr += var.obj_follower*yi[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*alpha_min[s];

                if(input.isNearOptMult() == true){
                    leader->getGRBModel()->addConstr(expr <= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[j], 1, expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                }else{
                    leader->getGRBModel()->addConstr(expr <= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[j], 1, expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                }
            }
        }
        
        // Lower bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0){
                leader->getGRBModel()->addConstr(yi[i] + dir*alpha_min[s] >= var.lb);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+j], 1, yi[i] + dir*alpha_min[s] <= var.lb);
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
                leader->getGRBModel()->addConstr(yi[i] + dir*alpha_min[s] <= var.ub);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,0)+j], 1, yi[i] + dir*alpha_min[s] >= var.ub);
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
                    if(s == 0) expr += constr.coeffs[var.id]*yi[i];    
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
                    if(s == 0) expr += var.obj_follower*yi[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*alpha_max[s];

                if(input.isNearOptMult() == true){
                    leader->getGRBModel()->addConstr(expr <= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[j], 1, expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                }else{
                    leader->getGRBModel()->addConstr(expr <= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[j], 1, expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                }
            }
        }

        // Lower bound constraints defining alpha_max.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0){
                leader->getGRBModel()->addConstr(yi[i] + dir*alpha_max[s] >= var.lb);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+j], 1, yi[i] + dir*alpha_max[s] <= var.lb);
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
                leader->getGRBModel()->addConstr(yi[i] + dir*alpha_max[s] <= var.ub);
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+instance.getGeneralNbBdlbsAlpha(pr,s,1)+j], 1, yi[i] + dir*alpha_max[s] >= var.ub);
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
                    if(s == 0) expr += constr.coeffs[var.id]*yi[i];    
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
                    if(s == 0) expr += var.obj_follower*yi[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*alpha_min[s];
                if(input.isNearOptMult() == true){
                    leader->getGRBModel()->addConstr(expr + slacks[s][j] == (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                }else{
                    leader->getGRBModel()->addConstr(expr + slacks[s][j] == fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                }
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
                    if(s == 0) expr += constr.coeffs[var.id]*yi[i];    
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
                    if(s == 0) expr += var.obj_follower*yi[i];    
                    else expr += var.obj_follower*y[s-1][i];
                }
                expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*alpha_max[s];
                if(input.isNearOptMult() == true){
                    leader->getGRBModel()->addConstr(expr + slacks[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j] == (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                }else{
                    leader->getGRBModel()->addConstr(expr + slacks[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j] == fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                }
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DepGeneralFollowerSolver::evaluate(double & mean, double & variance, double & f_mean, double & f_variance){
    double x_obj_leader = 0.0;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);

    double eval, var_eval, f_eval, f_var_eval, ess, rhat;
    instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,leader->getX_(),yi_,Fs_,Fw_,fs_,input.getTypeDepGeneral(),input.getNbIntervalsGeneral(),input.getScalingGeneral());
    
    mean = x_obj_leader + eval;
    variance = var_eval;

    f_mean = f_eval;
    f_variance = f_var_eval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DepGeneralFollowerSolver::computeStrongWeakInteriorSolutions(){
    // Compute both strong and weak solutions.
    AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,true,false);
}

