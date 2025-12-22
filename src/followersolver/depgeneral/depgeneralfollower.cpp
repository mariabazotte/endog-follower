#include "depgeneralfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

DepGeneralFollowerSolver::DepGeneralFollowerSolver(const Input & input, Instance & instance, LeaderSolver *leader, int pr) :
                                    AbstractFollowerSolver(input,instance,leader), pr(pr), 
                                    metropolis_hastings(input.metropolisHastings()) { 
    // Initialize pointers to NULL.
    y.resize(instance.getNbScenarios(),NULL);
    b_alpha_min.resize(instance.getNbScenarios(),NULL);
    b_alpha_max.resize(instance.getNbScenarios(),NULL);
    bound_alpha_min.resize(instance.getNbScenarios(),NULL);
    bound_alpha_max.resize(instance.getNbScenarios(),NULL);
    
    alpha_ = new double[instance.getNbScenarios()];
    alpha_min_ = new double[instance.getNbScenarios()];
    alpha_max_ = new double[instance.getNbScenarios()];
    b_alpha_min_ = new int[instance.getNbScenarios()];
    b_alpha_max_ = new int[instance.getNbScenarios()];
    y_.resize(instance.getNbScenarios(),NULL);
    bound_alpha_min_.resize(instance.getNbScenarios(),NULL);
    bound_alpha_max_.resize(instance.getNbScenarios(),NULL);
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        y_[s] = new double[instance.getModel()->nb_follower_vars];
    }
    
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
        alpha[s].set(GRB_DoubleAttr_LB,-instance.getStepSizeInterval()*numeric);
        alpha[s].set(GRB_DoubleAttr_UB,instance.getStepSizeInterval()*numeric);
        alpha[s].set(GRB_StringAttr_VarName,"alpha_"+std::to_string(s));

        alpha_min[s].set(GRB_DoubleAttr_LB,-instance.getStepSizeInterval()*numeric);
        alpha_min[s].set(GRB_DoubleAttr_UB,0.0);
        alpha_min[s].set(GRB_StringAttr_VarName,"alpha_min_"+std::to_string(s));

        alpha_max[s].set(GRB_DoubleAttr_LB,0.0);
        alpha_max[s].set(GRB_DoubleAttr_UB,instance.getStepSizeInterval()*numeric);
        alpha_max[s].set(GRB_StringAttr_VarName,"alpha_max_"+std::to_string(s));
    }

    if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
        if(metropolis_hastings == true){
            // Acceptance variables.
            q = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_BINARY);
        }

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

    // Define sequence.
    defineSequence();

    if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
        // Define constraints for interval coefficients.
        defineIntCoeffConstrs();

        // Define sequence according to probability.
        if(metropolis_hastings == true) 
            defineSequenceWithMetropolisHastings();
        else defineSequenceWithSliceSampling();
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

void DepGeneralFollowerSolver::defineIntCoeffConstrs(){
    // Define interval for probability distribution according to follower optimal value.
    GRBVar fsp = leader->getGRBModel()->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    leader->getGRBModel()->addConstr(fsp == (instance.getModel()->follower_ub - fs) / (instance.getModel()->follower_ub - instance.getModel()->follower_lb)); 
    
    // Define w_k.
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

        leader->getGRBModel()->addConstr(w_delta[k] <= (Fw -Fs));
        leader->getGRBModel()->addConstr(w_delta[k] >= (Fw -Fs) - (instance.getModel()->leader_ub - instance.getModel()->leader_lb)*(1.0 - w[k]));
        leader->getGRBModel()->addConstr(w_delta[k] <= (instance.getModel()->leader_ub - instance.getModel()->leader_lb)*w[k]);
    }
}

// Function to define sequence of points y according to the chosen sampling method.
// The sequence starts with the chebyshev center defined by variables yi.
void DepGeneralFollowerSolver::defineSequence(){
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::Neutral){
        // Accept all new points.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                if(s == 0) leader->getGRBModel()->addConstr(y[s][i] == yi[i] + (alpha[s]/numeric)*instance.getGeneralDirection(pr,s,i));
                else leader->getGRBModel()->addConstr(y[s][i] == y[s-1][i] + (alpha[s]/numeric)*instance.getGeneralDirection(pr,s,i));
            }
        }
    }else{
        if(metropolis_hastings == false){
            // Accept all new points.
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    if(s == 0) leader->getGRBModel()->addConstr(y[s][i] == yi[i] + (alpha[s]/numeric)*instance.getGeneralDirection(pr,s,i));
                    else leader->getGRBModel()->addConstr(y[s][i] == y[s-1][i] + (alpha[s]/numeric)*instance.getGeneralDirection(pr,s,i));
                }
            }
        }else{
            // Accept new points according to metropolis hastings rule.
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    if(s == 0){
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == yi[i] + (alpha[s]/numeric)*instance.getGeneralDirection(pr,s,i));
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == yi[i]);
                    }else{
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == y[s-1][i] + (alpha[s]/numeric)*instance.getGeneralDirection(pr,s,i));
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == y[s-1][i]);
                    }
                }
            }
        }
    }
}

void DepGeneralFollowerSolver::defineSequenceWithMetropolisHastings(){
    // Define metropolis hastings acceptance variables according to distribution.
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional){
        // Define acceptance variables for proportional case.
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
            expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*(alpha[s]/numeric);

            leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
            leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.1);
        }
    }

    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragile){
        // Define acceptance variables in fragile case.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            double r = instance.getGeneralAccep(pr,s);

            GRBLinExpr expr = 0;
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k)
                expr += (std::log(r) / input.getIntCoeffGeneral(k))*w_delta[k];

            expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*(alpha[s]/numeric);

            leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
            leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.0 + input.getEpsBigM());
        }
    }

    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragilePower || 
       input.getTypeDepGeneral() == Input::TypesDepGeneral::GenStrongPower){
        // Define acceptance variables in fragile and strong power cases.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            
            // In this case we need more auxiliar variables.
            GRBLinExpr leader_obj_s = 0.0;
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                BilevelVariable var = instance.getModel()->follower_vars[i];
                if(s == 0) leader_obj_s += var.obj_leader*yi[i];    
                else leader_obj_s += var.obj_leader*y[s-1][i];
            }

            GRBVar * w_s = leader->getGRBModel()->addVars(input.getNbIntervalsGeneral(),GRB_CONTINUOUS);
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                leader->getGRBModel()->addGenConstrIndicator(w[k], 0.0, w_s[k] == 0.0);
                leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, w_s[k] == leader_obj_s - (Fw + Fs)/2.0);
            }

            // Now defining constraint.
            double r = instance.getGeneralAccep(pr,s);

            GRBLinExpr expr = 0;
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragilePower){
                    expr += ((std::pow(r, 1.0/input.getIntCoeffPowerFrgGeneral(k)) - 1.0) / (2.0 * input.getIntCoeffSignalGeneral(k))) * w_delta[k];
                    expr += (std::pow(r, 1.0/input.getIntCoeffPowerFrgGeneral(k)) - 1.0) * w_s[k];
                }
                if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenStrongPower){
                    expr += ((std::pow(r, 1.0/input.getIntCoeffPowerStrGeneral(k)) - 1.0) / (2.0 * input.getIntCoeffSignalGeneral(k))) * w_delta[k];
                    expr += (std::pow(r, 1.0/input.getIntCoeffPowerStrGeneral(k)) - 1.0) * w_s[k];
                }
            }
            expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*(alpha[s]/numeric);
            
            leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
            leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.0 + input.getEpsBigM());
        }
    }
}

void DepGeneralFollowerSolver::defineSequenceWithSliceSampling(){
    for(int s = 0; s < instance.getNbScenarios(); ++s){  
        int index_min = instance.getGeneralNbConstrsAlpha(pr,s,0) + instance.getGeneralNbBdlbsAlpha(pr,s,0) + instance.getGeneralNbBdubsAlpha(pr,s,0) + 1;
        int index_max = instance.getGeneralNbConstrsAlpha(pr,s,1) + instance.getGeneralNbBdlbsAlpha(pr,s,1) + instance.getGeneralNbBdubsAlpha(pr,s,1) + 1;

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional){
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                double r = instance.getGeneralAccep(pr,s);

                GRBLinExpr expr = w_delta[k]*((r - 1.0) * input.getMaxIntCoeffGeneral()) / (2.0 * input.getIntCoeffGeneral(k) * instance.getGeneralCoeffAlphaObjLeader(pr,s));
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr += ((r - 1.0)*var.obj_leader*yi[i])/instance.getGeneralCoeffAlphaObjLeader(pr,s);    
                    else expr += ((r - 1.0)*var.obj_leader*y[s-1][i])/instance.getGeneralCoeffAlphaObjLeader(pr,s);
                }
                expr -= ((r - 1.0) * (Fw + Fs)) / (2.0 * instance.getGeneralCoeffAlphaObjLeader(pr,s));

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) >= (1e-12)){
                    // Constraint defines alpha_min.
                    if(relaxation == false) leader->getGRBModel()->addConstr(b_alpha_max[s][index_max] <= 1.0 - w[k]);
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_min[s]/numeric) >= expr);
                    if(relaxation == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_min], 1.0, (alpha_min[s]/numeric) <= expr);
                }

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k) <= (-1e-12))){
                    // Constraint defines alpha_max.   
                    if(relaxation == false) leader->getGRBModel()->addConstr(b_alpha_min[s][index_min] <= 1.0 - w[k]); 
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_max[s]/numeric) <= expr);
                    if(relaxation == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_max], 1.0, (alpha_max[s]/numeric) >= expr);
                }
            }
        }

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragile){
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                GRBLinExpr expr = (w_delta[k]*std::log(instance.getGeneralAccep(pr,s))) / 
                                  (input.getIntCoeffGeneral(k)*instance.getGeneralCoeffAlphaObjLeader(pr,s));

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) >= (1e-12)){
                    if(relaxation == false) leader->getGRBModel()->addConstr(b_alpha_max[s][index_max] <= 1.0 - w[k]);
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_min[s]/numeric) >= expr);
                    if(relaxation == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_min], 1.0, (alpha_min[s]/numeric) <= expr); 
                }

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k) <= (-1e-12))){
                    if(relaxation == false) leader->getGRBModel()->addConstr(b_alpha_min[s][index_min] <= 1.0 - w[k]);
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_max[s]/numeric) <= expr);
                    if(relaxation == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_max], 1.0, (alpha_max[s]/numeric) >= expr);
                }
            }
        }

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragilePower || 
           input.getTypeDepGeneral() == Input::TypesDepGeneral::GenStrongPower){

            GRBLinExpr leader_obj_s = 0.0;
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                BilevelVariable var = instance.getModel()->follower_vars[i];
                if(s == 0) leader_obj_s += var.obj_leader*yi[i];    
                else leader_obj_s += var.obj_leader*y[s-1][i];
            }

            GRBVar * w_s = leader->getGRBModel()->addVars(input.getNbIntervalsGeneral(),GRB_CONTINUOUS);
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                leader->getGRBModel()->addGenConstrIndicator(w[k], 0.0, w_s[k] == 0.0);
                leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, w_s[k] == leader_obj_s - (Fw + Fs)/2.0);
            }

            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                double u = 0.0;
                if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragilePower)
                    u = (std::pow(instance.getGeneralAccep(pr,s), 1.0/input.getIntCoeffPowerFrgGeneral(k)) - 1.0);
                if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenStrongPower)
                    u = (std::pow(instance.getGeneralAccep(pr,s), 1.0/input.getIntCoeffPowerStrGeneral(k)) - 1.0);

                GRBLinExpr expr = (w_delta[k]*0.5*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) +
                                  (w_s[k]*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s));
                
                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) >= (1e-12)){
                    if(relaxation == false) leader->getGRBModel()->addConstr(b_alpha_max[s][index_max] <= 1.0 - w[k]);
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_min[s]/numeric) >= expr);
                    if(relaxation == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_min], 1.0, (alpha_min[s]/numeric) <= expr); 
                }

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k) <= (-1e-12))){
                    if(relaxation == false) leader->getGRBModel()->addConstr(b_alpha_min[s][index_min] <= 1.0 - w[k]);
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_max[s]/numeric) <= expr);
                    if(relaxation == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_max], 1.0, (alpha_max[s]/numeric) >= expr);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DepGeneralFollowerSolver::defineNewExplicitStep(){
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        // Define primal constraints and slacks.
        GRBVar * slacks = NULL;    
        GRBVar * slack_lbs = NULL;
        GRBVar * slack_ubs = NULL;
        if(s == 0) definePrimalConstrs(yi,slacks,slack_lbs,slack_ubs,false,"_"+std::to_string(s+1));
        else definePrimalConstrs(y[s-1],slacks,slack_lbs,slack_ubs,false,"_"+std::to_string(s+1));

        GRBVar slack_eps, f_eps;
        if(input.isFollowerNearOpt() == true){
            if(s == 0) defineFollowerNearOptimalObj(f_eps,yi,slack_eps,"_" + std::to_string(s+1));
            else defineFollowerNearOptimalObj(f_eps,y[s-1],slack_eps,"_" + std::to_string(s+1));
        } else {
            if(s == 0) defineFollowerOptimalObj(yi,"_" + std::to_string(s+1));
            else defineFollowerOptimalObj(y[s-1],"_" + std::to_string(s+1));
        }
 
        if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12 && relaxation == false){
            // Binary variable defining best lower bound on alpha_min.
            b_alpha_min[s] =  leader->getGRBModel()->addVars(getNbConstrsMin(s),GRB_BINARY);

            GRBLinExpr sum_expr = 0;
            for(int i = 0; i < getNbConstrsMin(s); ++i) sum_expr += b_alpha_min[s][i];
            leader->getGRBModel()->addConstr(sum_expr == 1);
        }
        if(instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12 && relaxation == false){
            // Binary variable defining best upper bound on alpha_max.
            b_alpha_max[s] =  leader->getGRBModel()->addVars(getNbConstrsMax(s),GRB_BINARY);

            GRBLinExpr sum_expr_max = 0;
            for(int i = 0; i < getNbConstrsMax(s); ++i) sum_expr_max += b_alpha_max[s][i];
            leader->getGRBModel()->addConstr(sum_expr_max == 1);
        }

        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,0,j);
            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];
                if(constr.sense == '<') {
                    leader->getGRBModel()->addConstr(instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_min[s]/numeric) <= slacks[c]);
                    if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12 && relaxation == false){
                        leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_min[s]/numeric) >= slacks[c] - 1e-6);
                    }
                }else if(constr.sense == '>'){ 
                    leader->getGRBModel()->addConstr(instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_min[s]/numeric) >= -slacks[c]);
                    if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12 && relaxation == false){
                        leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_min[s]/numeric) <= -slacks[c] + 1e-6);
                    }
                }
            }else{
                leader->getGRBModel()->addConstr(instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_min[s]/numeric) <= slack_eps);
                if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12 && relaxation == false){
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_min[s]/numeric) >= slack_eps - 1e-6);
                }
            }
        }
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,0,j);
            leader->getGRBModel()->addConstr(instance.getGeneralDirection(pr,s,i)*(alpha_min[s]/numeric) >= -slack_lbs[i]);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12 && relaxation == false){
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j], 1, instance.getGeneralDirection(pr,s,i)*(alpha_min[s]/numeric) <= -slack_lbs[i] + 1e-6);
            }
        }
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,0,j);
            leader->getGRBModel()->addConstr(instance.getGeneralDirection(pr,s,i)*(alpha_min[s]/numeric) <= slack_ubs[i]);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12 && relaxation == false){
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,0)+j], 1, instance.getGeneralDirection(pr,s,i)*(alpha_min[s]/numeric) >= slack_ubs[i] - 1e-6);
            }
        }
        if((instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12) && relaxation == false){
            int index = instance.getGeneralNbConstrsAlpha(pr,s,0) + 
                        instance.getGeneralNbBdlbsAlpha(pr,s,0) + 
                        instance.getGeneralNbBdubsAlpha(pr,s,0);
            leader->getGRBModel()->addConstr((alpha_min[s]/numeric) <= -instance.getStepSizeInterval()*b_alpha_min[s][index]);
        }
                
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,1); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,1,j);
            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];
                if(constr.sense == '<') {
                    leader->getGRBModel()->addConstr(instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_max[s]/numeric) <= slacks[c]);
                    if(instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12 && relaxation == false){
                        leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_max[s]/numeric) >= slacks[c] - 1e-6);
                    }
                }else if(constr.sense == '>'){ 
                    leader->getGRBModel()->addConstr(instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_max[s]/numeric) >= -slacks[c]);
                    if((instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12) && relaxation == false){
                        leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_max[s]/numeric) <= -slacks[c] + 1e-6);
                    }
                }
            }else{
                leader->getGRBModel()->addConstr(instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_max[s]/numeric) <= slack_eps);
                if((instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12) && relaxation == false){
                    leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_max[s]/numeric) >= slack_eps - 1e-6);
                }
            }
        }
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
            leader->getGRBModel()->addConstr(instance.getGeneralDirection(pr,s,i)*(alpha_max[s]/numeric) >= -slack_lbs[i]);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            if((instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12) && relaxation == false){
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][instance.getGeneralNbConstrsAlpha(pr,s,1)+j], 1, instance.getGeneralDirection(pr,s,i)*(alpha_max[s]/numeric) <= -slack_lbs[i] + 1e-6);
            }
        }
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,1,j);
            leader->getGRBModel()->addConstr(instance.getGeneralDirection(pr,s,i)*(alpha_max[s]/numeric) <= slack_ubs[i]);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            if((instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12) && relaxation == false){
                leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][instance.getGeneralNbConstrsAlpha(pr,s,1)+instance.getGeneralNbBdlbsAlpha(pr,s,1)+j], 1, instance.getGeneralDirection(pr,s,i)*(alpha_max[s]/numeric) >= slack_ubs[i] - 1e-6);
            }
        }
        if((instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12) && relaxation == false){
            int index_max = instance.getGeneralNbConstrsAlpha(pr,s,1) + 
                        instance.getGeneralNbBdlbsAlpha(pr,s,1) + 
                        instance.getGeneralNbBdubsAlpha(pr,s,1);
            leader->getGRBModel()->addConstr((alpha_max[s]/numeric) >= instance.getStepSizeInterval()*b_alpha_max[s][index_max]);
        }
    }
}

void DepGeneralFollowerSolver::defineExplicitStep(){
    // Define constraints for alpha_min.
    for(int s = 0; s < instance.getNbScenarios(); ++s) {

        int nb_constrs_min = getNbConstrsMin(s);

        // Binary variable defining best lower bound on alpha_min.
        if(relaxation == false) b_alpha_min[s] = leader->getGRBModel()->addVars(nb_constrs_min,GRB_BINARY);

        // Continuous variable defining bounds.
        if(relaxation == false && lazy_callback == true){
            bound_alpha_min[s] = leader->getGRBModel()->addVars(nb_constrs_min,GRB_CONTINUOUS);
            for(int i = 0; i < nb_constrs_min; ++i){
                bound_alpha_min[s][i].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
                bound_alpha_min[s][i].set(GRB_DoubleAttr_UB,0.0);
            }
        }

        if(relaxation == false){
            GRBLinExpr sum_expr = 0;
            for(int i = 0; i < nb_constrs_min; ++i) sum_expr += b_alpha_min[s][i];
            leader->getGRBModel()->addConstr(sum_expr == 1);

            std::vector<GRBVar> v(nb_constrs_min);
            std::vector<double> w(nb_constrs_min);
            for (int i = 0; i < nb_constrs_min; ++i) {
                v[i] = b_alpha_min[s][i];
                w[i] = static_cast<double>(i + 1);
            }
            leader->getGRBModel()->addSOS(v.data(), w.data(), nb_constrs_min, GRB_SOS_TYPE1);
        }

        // GRBVar * slack_alpha_min = leader->getGRBModel()->addVars(nb_constrs_min,GRB_CONTINUOUS);
        // for(int i = 0; i < nb_constrs_min; ++i){
        //     GRBVar v[2] = {slack_alpha_min[i], b_alpha_min[s][i]};
        //     double w[2] = {1.0, 2.0};
        //     leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);  
        // }

        // Constraints impacting alpha_min.
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,0,j);
            
            // Follower constraints.
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
                GRBLinExpr bound_expr = (constr.rhs - expr)/instance.getGeneralCoeffAlpha(pr,s,c);
                expr += instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_min[s]/numeric);
                
                if(constr.sense == '<'){
                    leader->getGRBModel()->addConstr(expr <= constr.rhs);
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr >= constr.rhs - 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr >= constr.rhs - (constr.rhs - constr.bound)*(1.0 - b_alpha_min[s][j])); 
                    // leader->getGRBModel()->addConstr(expr + slack_alpha_min[j] == constr.rhs);
                }
                if(constr.sense == '>'){
                    leader->getGRBModel()->addConstr(expr >= constr.rhs);
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr <= constr.rhs + 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr <= constr.rhs + (constr.bound - constr.rhs)*(1.0 - b_alpha_min[s][j]));
                    // leader->getGRBModel()->addConstr(expr - slack_alpha_min[j] == constr.rhs);
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
                GRBLinExpr bound_expr = -(expr)/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_min[s]/numeric);

                if(input.isNearOptMult() == true){
                    bound_expr += ((1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub)/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr(expr <= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - instance.getModel()->follower_ub*(1.0 - b_alpha_min[s][j]));
                    // leader->getGRBModel()->addConstr(expr + slack_alpha_min[j] == (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                }else{
                    bound_expr += (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr(expr <= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - (instance.getModel()->follower_ub + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))*(1.0 - b_alpha_min[s][j]));
                    // leader->getGRBModel()->addConstr(expr + slack_alpha_min[j] == fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                }
            }
        }
        // Lower bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            int index_lb = instance.getGeneralNbConstrsAlpha(pr,s,0);
            if(s == 0){
                leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) >= var.lb);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_min[s][index_lb+j]/numeric) == var.lb);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_lb+j], 1, yi[i] + dir*(alpha_min[s]/numeric) <= var.lb + 1e-6);
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_min[s][index_lb+j]));
                // leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) - slack_alpha_min[index_lb+j] == var.lb);
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) >= var.lb);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_min[s][index_lb+j]/numeric) == var.lb);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_lb+j], 1, y[s-1][i] + dir*(alpha_min[s]/numeric) <= var.lb + 1e-6);
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_min[s][index_lb+j]));
                // leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) - slack_alpha_min[index_lb+j] == var.lb);
            }  
        }
        // Upper bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            
            double dir = instance.getGeneralDirection(pr,s,i);
            int index_ub = instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,0);
            if(s == 0){
                leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) <= var.ub);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_min[s][index_ub+j]/numeric) == var.ub);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_ub+j], 1, yi[i] + dir*(alpha_min[s]/numeric) >= var.ub - 1e-6);
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) >= var.ub - (var.ub-var.lb)*(1.0 - b_alpha_min[s][index_ub+j]));
                // leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) + slack_alpha_min[index_ub+j] == var.ub); 
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) <= var.ub);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_min[s][index_ub+j]/numeric) == var.ub);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_ub+j], 1, y[s-1][i] + dir*(alpha_min[s]/numeric) >= var.ub - 1e-6);
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) >= var.ub - (var.ub - var.lb)*(1.0 - b_alpha_min[s][index_ub+j]));
                // leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) + slack_alpha_min[index_ub+j] == var.ub);
            }
        }
        // Lower bound on alpha_min.
        int index = instance.getGeneralNbConstrsAlpha(pr,s,0) + 
                    instance.getGeneralNbBdlbsAlpha(pr,s,0) + 
                    instance.getGeneralNbBdubsAlpha(pr,s,0);
        leader->getGRBModel()->addConstr((alpha_min[s]/numeric) >= -instance.getStepSizeInterval());
        if(relaxation == false) leader->getGRBModel()->addConstr((alpha_min[s]/numeric) <= -instance.getStepSizeInterval()*b_alpha_min[s][index]);
        if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][index]/numeric) == -instance.getStepSizeInterval());
        // leader->getGRBModel()->addConstr((alpha_min[s]/numeric) - slack_alpha_min[index] == -instance.getStepSizeInterval());
    }

    // Define constraints for alpha_max.
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        int nb_constrs_max = getNbConstrsMax(s);

        if(relaxation == false) b_alpha_max[s] =  leader->getGRBModel()->addVars(nb_constrs_max,GRB_BINARY);

        // Continuous variable defining bounds.
        if(relaxation == false && lazy_callback == true){
            bound_alpha_max[s] = leader->getGRBModel()->addVars(nb_constrs_max,GRB_CONTINUOUS);
            for(int i = 0; i < nb_constrs_max; ++i){
                bound_alpha_max[s][i].set(GRB_DoubleAttr_LB,0.0);
                bound_alpha_max[s][i].set(GRB_DoubleAttr_UB,GRB_INFINITY);
            }
        }

        if(relaxation == false){
            GRBLinExpr sum_expr = 0;
            for(int i = 0; i < nb_constrs_max; ++i) sum_expr += b_alpha_max[s][i];
            leader->getGRBModel()->addConstr(sum_expr == 1);

            std::vector<GRBVar> v(nb_constrs_max);
            std::vector<double> w(nb_constrs_max);
            for (int i = 0; i < nb_constrs_max; ++i) {
                v[i] = b_alpha_max[s][i];
                w[i] = static_cast<double>(i + 1);
            }
            leader->getGRBModel()->addSOS(v.data(), w.data(), nb_constrs_max, GRB_SOS_TYPE1);
        }

        // GRBVar * slack_alpha_max = leader->getGRBModel()->addVars(nb_constrs_max,GRB_CONTINUOUS);
        // for(int i = 0; i < nb_constrs_max; ++i){
        //     GRBVar v[2] = {slack_alpha_max[i], b_alpha_max[s][i]};
        //     double w[2] = {1.0, 2.0};
        //     leader->getGRBModel()->addSOS(v, w, 2, GRB_SOS_TYPE1);  
        // }

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
                GRBLinExpr bound_expr = (constr.rhs - expr)/instance.getGeneralCoeffAlpha(pr,s,c);
                expr += instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_max[s]/numeric);
                
                if(constr.sense == '<'){
                    leader->getGRBModel()->addConstr(expr <= constr.rhs);
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr >= constr.rhs - 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr >= constr.rhs - (constr.rhs - constr.bound)*(1.0 - b_alpha_max[s][j]));
                    // leader->getGRBModel()->addConstr(expr + slack_alpha_max[j] == constr.rhs);
                }
                if(constr.sense == '>'){
                    leader->getGRBModel()->addConstr(expr >= constr.rhs);
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr <= constr.rhs + 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr <= constr.rhs + (constr.bound - constr.rhs)*(1.0 - b_alpha_max[s][j]));
                    // leader->getGRBModel()->addConstr(expr - slack_alpha_max[j] == constr.rhs);
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
                GRBLinExpr bound_expr = -(expr)/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_max[s]/numeric);

                if(input.isNearOptMult() == true){
                    bound_expr += ((1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub)/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr(expr <= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub);
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - instance.getModel()->follower_ub*(1.0 - b_alpha_max[s][j]));
                }else{
                    bound_expr += (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr(expr <= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                    if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - 1e-6);
                    if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - (instance.getModel()->follower_ub + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))*(1.0 - b_alpha_max[s][j]));
                    // leader->getGRBModel()->addConstr(expr + slack_alpha_max[j] == fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                }
            }
        }
        // Lower bound constraints defining alpha_max.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            int index_lb = instance.getGeneralNbConstrsAlpha(pr,s,1);
            if(s == 0){
                leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) >= var.lb);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_max[s][index_lb+j]/numeric) == var.lb);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_lb+j], 1, yi[i] + dir*(alpha_max[s]/numeric) <= var.lb + 1e-6);
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_lb+j]));
                // leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) - slack_alpha_max[index_lb+j] == var.lb);                 
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) >= var.lb);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_max[s][index_lb+j]/numeric) == var.lb);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_lb+j], 1, y[s-1][i] + dir*(alpha_max[s]/numeric) <= var.lb + 1e-6);
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_lb+j]));
                // leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) - slack_alpha_max[index_lb+j] == var.lb);
            }  
        }
        // Upper bound constraints defining alpha_max.
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            int index_ub = instance.getGeneralNbConstrsAlpha(pr,s,1)+instance.getGeneralNbBdlbsAlpha(pr,s,1);
            if(s == 0){
                leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) <= var.ub);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_max[s][index_ub+j]/numeric) == var.ub);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_ub+j], 1, yi[i] + dir*(alpha_max[s]/numeric) >= var.ub - 1e-6);    
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) >= var.ub - (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_ub+j]));
                // leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) + slack_alpha_max[index_ub+j] == var.ub);
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) <= var.ub);
                if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_max[s][index_ub+j]/numeric) == var.ub);
                if(relaxation == false && lazy_callback == false && bigm == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_ub+j], 1, y[s-1][i] + dir*(alpha_max[s]/numeric) >= var.ub - 1e-6);
                if(relaxation == false && lazy_callback == false && bigm == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) >= var.ub - (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_ub+j]));
                // leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) + slack_alpha_max[index_ub+j] == var.ub);
            }
        }
        // Upper bound on alpha_max.
        int index = instance.getGeneralNbConstrsAlpha(pr,s,1) + 
                    instance.getGeneralNbBdlbsAlpha(pr,s,1) + 
                    instance.getGeneralNbBdubsAlpha(pr,s,1);
        leader->getGRBModel()->addConstr((alpha_max[s]/numeric) <= instance.getStepSizeInterval());
        if(relaxation == false) leader->getGRBModel()->addConstr((alpha_max[s]/numeric) >= instance.getStepSizeInterval()*b_alpha_max[s][index]);
        if(relaxation == false && lazy_callback == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][index]/numeric) == instance.getStepSizeInterval());
        // leader->getGRBModel()->addConstr((alpha_max[s]/numeric) + slack_alpha_max[index] == instance.getStepSizeInterval());
    }  
}

void DepGeneralFollowerSolver::defineExplicitStepMaxMinConstr(){
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        // Alpha_min.
        int nb_constrs_min = instance.getGeneralNbConstrsAlpha(pr,s,0) + 
                             instance.getGeneralNbBdlbsAlpha(pr,s,0) + 
                             instance.getGeneralNbBdubsAlpha(pr,s,0) + 1;
        if((metropolis_hastings == false) && (input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral)) 
            nb_constrs_min += 1;

        GRBVar * gen_alpha_min = leader->getGRBModel()->addVars(nb_constrs_min,GRB_CONTINUOUS);
        for(int i = 0; i < nb_constrs_min; ++i){
            gen_alpha_min[i].set(GRB_DoubleAttr_UB,0.0);
            gen_alpha_min[i].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
        }

        // Constraints impacting alpha_min.
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,0,j);

            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];

                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                    BilevelVariable var = instance.getModel()->leader_vars[i];  
                    expr -= constr.coeffs[var.id]*leader->getX(i)/instance.getGeneralCoeffAlpha(pr,s,c);   
                }
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr -= constr.coeffs[var.id]*yi[i]/instance.getGeneralCoeffAlpha(pr,s,c);    
                    else expr -= constr.coeffs[var.id]*y[s-1][i]/instance.getGeneralCoeffAlpha(pr,s,c); 
                }
                expr += constr.rhs/instance.getGeneralCoeffAlpha(pr,s,c);
                leader->getGRBModel()->addConstr((gen_alpha_min[j]/numeric) == expr);   
            }else{
                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr -= var.obj_follower*yi[i]/instance.getGeneralCoeffAlphaObjFollower(pr,s);    
                    else expr -= var.obj_follower*y[s-1][i]/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                }
                if(input.isNearOptMult() == true){
                    expr += ((1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub)/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr((gen_alpha_min[j]/numeric) == expr);
                }
                else{
                    expr += (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr((gen_alpha_min[j]/numeric) == expr);
                }
            }
        }
        // Lower bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0) leader->getGRBModel()->addConstr((gen_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+j]/numeric) == (var.lb - yi[i])/dir);
            else leader->getGRBModel()->addConstr((gen_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+j]/numeric) == (var.lb - y[s-1][i])/dir); 
        }
        // Upper bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0) leader->getGRBModel()->addConstr((gen_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,0)+j]/numeric) == (var.ub - yi[i])/dir);    
            else leader->getGRBModel()->addConstr((gen_alpha_min[instance.getGeneralNbConstrsAlpha(pr,s,0)+instance.getGeneralNbBdlbsAlpha(pr,s,0)+j]/numeric) == (var.ub - y[s-1][i])/dir);
        }
        leader->getGRBModel()->addGenConstrMax(alpha_min[s], gen_alpha_min, nb_constrs_min, -instance.getStepSizeInterval()*numeric);
        
        // Alpha_max.
        int nb_constrs_max = instance.getGeneralNbConstrsAlpha(pr,s,1) + 
                                 instance.getGeneralNbBdlbsAlpha(pr,s,1) + 
                                 instance.getGeneralNbBdubsAlpha(pr,s,1) + 1;
            if((metropolis_hastings == false) && (input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral)) 
                nb_constrs_max += 1;
            
        GRBVar * gen_alpha_max = leader->getGRBModel()->addVars(nb_constrs_max,GRB_CONTINUOUS);
        for(int i = 0; i < nb_constrs_max; ++i){
            gen_alpha_max[i].set(GRB_DoubleAttr_LB,0.0);
            gen_alpha_max[i].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        }

        // Constraints impacting alpha_max.
        for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,1); ++j){
            int c = instance.getGeneralConstrsAlpha(pr,s,1,j);
            
            // Usual follower constraints.
            if(c != instance.getModel()->nb_follower_constrs){
                BilevelConstraint constr = instance.getModel()->follower_constrs[c];

                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                    BilevelVariable var = instance.getModel()->leader_vars[i]; 
                    expr -= constr.coeffs[var.id]*leader->getX(i)/instance.getGeneralCoeffAlpha(pr,s,c);  
                }
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr -= constr.coeffs[var.id]*yi[i]/instance.getGeneralCoeffAlpha(pr,s,c);    
                    else expr -= constr.coeffs[var.id]*y[s-1][i]/instance.getGeneralCoeffAlpha(pr,s,c);  
                }  
                expr += constr.rhs/instance.getGeneralCoeffAlpha(pr,s,c);
                leader->getGRBModel()->addConstr((gen_alpha_max[j]/numeric) == expr);
            }
            // Near optimality constraint.
            else{
                GRBLinExpr expr = 0;
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    BilevelVariable var = instance.getModel()->follower_vars[i];
                    if(s == 0) expr -= var.obj_follower*yi[i]/instance.getGeneralCoeffAlphaObjFollower(pr,s);    
                    else expr -= var.obj_follower*y[s-1][i]/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                }
                if(input.isNearOptMult() == true){
                    expr += ((1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub)/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr((gen_alpha_max[j]/numeric) == expr);
                }else{
                    expr += (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr((gen_alpha_max[j]/numeric) == expr);
                }
            }
        }
        // Lower bound constraints defining alpha_max.
        for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0) leader->getGRBModel()->addConstr((gen_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+j]/numeric) == (var.lb -yi[i])/dir);             
            else leader->getGRBModel()->addConstr((gen_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+j]/numeric) == (var.lb -y[s-1][i])/dir); 
        }
        // Upper bound constraints defining alpha_max.
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,1,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];

            double dir = instance.getGeneralDirection(pr,s,i);
            if(s == 0) leader->getGRBModel()->addConstr((gen_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+instance.getGeneralNbBdlbsAlpha(pr,s,1)+j]/numeric) == (var.ub - yi[i])/dir);
            else leader->getGRBModel()->addConstr((gen_alpha_max[instance.getGeneralNbConstrsAlpha(pr,s,1)+instance.getGeneralNbBdlbsAlpha(pr,s,1)+j]/numeric) == (var.ub - y[s-1][i])/dir);
        }
        leader->getGRBModel()->addGenConstrMin(alpha_max[s], gen_alpha_max, nb_constrs_max, instance.getStepSizeInterval());
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
        if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12){
            GRBLinExpr expr = 0;
            for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,0); ++j){
                int c = instance.getGeneralConstrsAlpha(pr,s,0,j);

                if(c != instance.getModel()->nb_follower_constrs){
                    BilevelConstraint constr = instance.getModel()->follower_constrs[c];
                    if(constr.sense == '<') expr += -duals[s][j]*instance.getGeneralCoeffAlpha(pr,s,c);
                    else if(constr.sense == '>') expr += duals[s][j]*instance.getGeneralCoeffAlpha(pr,s,c);
                }else{
                    expr += -duals[s][j]*instance.getGeneralCoeffAlphaObjFollower(pr,s);
                }
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
    }
    
    // Dual constraint for alpha_max. (Primal is a maximization problem).
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        if(instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12){
            GRBLinExpr expr = 0;
            for(int j = 0; j < instance.getGeneralNbConstrsAlpha(pr,s,1); ++j){
                int c = instance.getGeneralConstrsAlpha(pr,s,1,j);

                if(c != instance.getModel()->nb_follower_constrs){
                    BilevelConstraint constr = instance.getModel()->follower_constrs[c];
                    if(constr.sense == '<') expr += duals[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j]*instance.getGeneralCoeffAlpha(pr,s,c);
                    else if(constr.sense == '>') expr += -duals[s][instance.getGeneralNbConstrsAlpha(pr,s,0)+j]*instance.getGeneralCoeffAlpha(pr,s,c);
                }else{
                    expr += duals[s][j]*instance.getGeneralCoeffAlphaObjFollower(pr,s);
                }
            }
            for(int j = 0; j < instance.getGeneralNbBdlbsAlpha(pr,s,1); ++j) {
                int i = instance.getGeneralBdlbsAlpha(pr,s,1,j);
                expr += -duals_lbs[s][instance.getGeneralNbBdlbsAlpha(pr,s,0)+j]*instance.getGeneralDirection(pr,s,i);
            }
            for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,1); ++j) {
                int i = instance.getGeneralBdubsAlpha(pr,s,1,j) ;
                expr += duals_ubs[s][instance.getGeneralNbBdubsAlpha(pr,s,0)+j]*instance.getGeneralDirection(pr,s,i);
            }
            expr += -duals_max[s][0] + duals_max[s][1];
            leader->getGRBModel()->addConstr(expr == 1.0);
        }
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
        if(instance.getGeneralCoeffAlphaObjLeader(pr,s) <= -1e-12){
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
        }

        // alpha_max
        if(instance.getGeneralCoeffAlphaObjLeader(pr,s) >= 1e-12){
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DepGeneralFollowerSolver::upd_solution(){
    AbstractFollowerSolver::upd_solution();
    

    double x_obj_leader = 0.0;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);



    double y_obj_leader = 0.0; 
    instance.evaluateSAADepGeneral(y_obj_leader,alpha_,alpha_min_,b_alpha_min_,alpha_max_,b_alpha_max_,y_,
                                    pr,leader->getX_(),yi_,Fs_,Fw_,fs_,
                                    input.getTypeDepGeneral(),
                                    input.getNbIntervalsGeneral(),
                                    input.getScalingGeneral());
    
    std::vector<double*> sample(instance.getNbScenarios(),NULL);
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        std::cout << "scenario " << s << std::endl;
        std::cout << "alpha_min: " << alpha_min_[s] << "|| " << alpha_min[s].get(GRB_DoubleAttr_X)/numeric << std::endl;
        std::cout << "alpha_max: " << alpha_max_[s] << "|| " << alpha_max[s].get(GRB_DoubleAttr_X)/numeric << " "  << std::endl;
        std::cout << "alpha: " << alpha_[s] << "||" << alpha[s].get(GRB_DoubleAttr_X)/numeric << " "  << std::endl;
        for(int j = 0; j < getNbConstrsMin(s); ++j){
            std::cout << bound_alpha_min[s][j].get(GRB_DoubleAttr_X) << " ";
        }
        std::cout << "\n";
        for(int j = 0; j < getNbConstrsMax(s); ++j){
            std::cout << bound_alpha_max[s][j].get(GRB_DoubleAttr_X) << " ";
        }
        std::cout << "\n";
        sample[s] = leader->getGRBModel()->get(GRB_DoubleAttr_X,y[s],instance.getModel()->nb_follower_vars);
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            if(s > 0){
                std::cout << y_[s][i] << " " << (y_[s-1][i] + instance.getGeneralDirection(pr,s,i)*alpha_[s]) << "||" << sample[s][i] << " " << (sample[s-1][i] + instance.getGeneralDirection(pr,s,i)*(alpha[s].get(GRB_DoubleAttr_X)/numeric))<< std::endl; 
            }
        }
        
        // std::cout << "coeff:" << instance.getGeneralCoeffAlphaObjLeader(pr,s) << std::endl;

        // if(alpha_min[s].get(GRB_DoubleAttr_X) <= alpha_min_[s] - 1e-12) std::cout << "Error alpha_min" << std::endl;
        // if(alpha_max[s].get(GRB_DoubleAttr_X) >= alpha_max_[s] + 1e-12) std::cout << "Error alpha_max" << std::endl;
    }

    std::cout << "real objective: " << x_obj_leader+y_obj_leader << std::endl;
    exit(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DepGeneralFollowerSolver::callback(){
    try{
        if(where == GRB_CB_MIPSOL){
            std::cout << "Running callback\n";
            // Verifying any violated constraints.
            if(alpha_min_) delete[] alpha_min_;
            if(alpha_max_) delete[] alpha_max_;
            alpha_min_ = getSolution(alpha_min,instance.getNbScenarios());
            alpha_max_ = getSolution(alpha_min,instance.getNbScenarios());

            std::cout << "Solution found\n";
            
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                if(bound_alpha_min_[s]) delete[] bound_alpha_min_[s];
                if(bound_alpha_max_[s]) delete[] bound_alpha_max_[s];

                int nb_constrs_min = getNbConstrsMin(s);
                int nb_constrs_max = getNbConstrsMax(s);

                bound_alpha_min_[s] = getSolution(bound_alpha_min[s], nb_constrs_min);
                bound_alpha_max_[s] = getSolution(bound_alpha_max[s], nb_constrs_max);

                auto min_ = std::max_element(bound_alpha_min_[s], bound_alpha_min_[s] + nb_constrs_min);
                auto max_ = std::min_element(bound_alpha_max_[s], bound_alpha_max_[s] + nb_constrs_max);

                int min_it = std::distance(bound_alpha_min_[s], min_);
                int max_it = std::distance(bound_alpha_max_[s], max_);

                double min_val = *min_;
                double max_val = *max_;

                std::cout << alpha_min_[s] << " " << min_val << " " << bound_alpha_min_[s][min_it] << std::endl;
                std::cout << alpha_max_[s] << " " << max_val << " " << bound_alpha_max_[s][max_it] << std::endl;

                bool included = false;
                if(alpha_min_[s] >= min_val + 1e-12){
                    included = true;
                    GRBLinExpr expr = 0;
                    getExpression(expr, s, 0, min_it);
                    addLazy(expr >= 0);
                }
                if(alpha_max_[s] <= max_val - 1e-12){
                    included = true;
                    GRBLinExpr expr = 0;
                    getExpression(expr, s, 1, max_it);
                    addLazy(expr >= 0);
                }
                if(included == true) break;
            }

            if(x_) delete[] x_; 
            if(yi_) delete[] yi_;
            
            x_ = getSolution(leader->getX(),instance.getModel()->nb_leader_vars);
            yi_ = getSolution(yi,instance.getModel()->nb_follower_vars);
            fs_ = getSolution(fs);
            if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
                Fs_ = getSolution(Fs);
                Fw_ = getSolution(Fw);
            }

            double x_obj_leader = 0.0;
            for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
                x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*x_[i];

            std::fill(alpha_, alpha_ + instance.getNbScenarios(), -1.0);
            std::fill(alpha_min_, alpha_min_ + instance.getNbScenarios(), -1.0);
            std::fill(alpha_max_, alpha_max_ + instance.getNbScenarios(), -1.0);
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                std::fill(y_[s], y_[s] + instance.getModel()->nb_follower_vars, -1.0);
            }

            double y_obj_leader = 0.0;
            instance.evaluateSAADepGeneral(y_obj_leader,alpha_,alpha_min_,b_alpha_min_,alpha_max_,b_alpha_max_,y_,
                                          pr,x_,yi_,Fs_,Fw_,fs_,
                                          input.getTypeDepGeneral(),
                                          input.getNbIntervalsGeneral(),
                                          input.getScalingGeneral());

            std::cout << (x_obj_leader + y_obj_leader) <<" " << getDoubleInfo(GRB_CB_MIPSOL_OBJBST) << std::endl;
            if((x_obj_leader + y_obj_leader) <= getDoubleInfo(GRB_CB_MIPSOL_OBJBST)){
                std::cout << "Solution found: " << (x_obj_leader + y_obj_leader)  << std::endl;
                
                setSolution(leader->getX(), x_, instance.getModel()->nb_leader_vars);
                setSolution(yi, yi_, instance.getModel()->nb_follower_vars);
                setSolution(alpha, alpha_, instance.getNbScenarios());
                setSolution(alpha_min, alpha_min_, instance.getNbScenarios());
                setSolution(alpha_max, alpha_max_, instance.getNbScenarios());
                for(int s = 0; s < instance.getNbScenarios(); ++s){
                    setSolution(y[s], y_[s], instance.getModel()->nb_follower_vars);
                    for(int j = 0; j < getNbConstrsMin(s); ++j){
                        if(j == b_alpha_min_[s]) setSolution(b_alpha_min[s][j],1.0);
                        else setSolution(b_alpha_min[s][j],0.0);
                    }
                    for(int j = 0; j < getNbConstrsMax(s); ++j){
                        if(j == b_alpha_max_[s]) setSolution(b_alpha_max[s][j],1.0);
                        else setSolution(b_alpha_max[s][j],0.0);
                    }
                }

                for(int s = 0; s < instance.getNbScenarios(); ++s){
                    std::cout << "scenario " << s << std::endl;
                    std::cout << "alpha_min: " << alpha_min_[s] << std::endl;
                    std::cout << "alpha_max: " << alpha_max_[s] << std::endl;
                    std::cout << "alpha: " << alpha_[s] << "||" << std::endl;
                    std::cout << "\n";
                    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                        if(s > 0){
                            std::cout << y_[s][i] << " " << (y_[s-1][i] + instance.getGeneralDirection(pr,s,i)*alpha_[s]) << std::endl; 
                        }
                    }
                }
                
                // useSolution();
            }
            std::cout << "---------"<< std::endl;
        }
    }
    catch (GRBException e){
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...){
        std::cout << "Error during callback" << std::endl;
    }
}

void DepGeneralFollowerSolver::getExpression(GRBLinExpr & expr, int s, int min_or_max, int it){
    expr = 0;
    if(it < instance.getGeneralNbConstrsAlpha(pr, s, min_or_max)){
        int c = instance.getGeneralConstrsAlpha(pr, s, min_or_max, it);
    
        if(c != instance.getModel()->nb_follower_constrs){
            BilevelConstraint constr = instance.getModel()->follower_constrs[c];

            expr = 0;
            for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                BilevelVariable var = instance.getModel()->leader_vars[i];
                expr += constr.coeffs[var.id]*leader->getX(i);     
            }
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                BilevelVariable var = instance.getModel()->follower_vars[i];  
                if(s == 0) expr += constr.coeffs[var.id]*yi[i];    
                else expr += constr.coeffs[var.id]*y[s-1][i];
            }
            if(min_or_max == 0) expr += instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_min[s]/numeric);
            else if(min_or_max == 1) expr += instance.getGeneralCoeffAlpha(pr,s,c)*(alpha_max[s]/numeric);

            if(min_or_max == 0){
                if(constr.sense == '<') {
                    expr -= (constr.rhs - (constr.rhs - constr.bound)*(1.0 - b_alpha_min[s][it])); 
                }
                else if(constr.sense == '>') expr = (constr.rhs + (constr.bound - constr.rhs)*(1.0 - b_alpha_min[s][it]) - expr);
            }
            else if(min_or_max == 1){
                if(constr.sense == '<') expr -= (constr.rhs - (constr.rhs - constr.bound)*(1.0 - b_alpha_max[s][it]));
                else if(constr.sense == '>') expr = (constr.rhs + (constr.bound - constr.rhs)*(1.0 - b_alpha_max[s][it]) - expr);
            }
        }else{
            expr = 0;
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                BilevelVariable var = instance.getModel()->follower_vars[i];
                if(s == 0) expr += var.obj_follower*yi[i];    
                else expr += var.obj_follower*y[s-1][i];
            }
            if(min_or_max == 0) expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_min[s]/numeric);
            else if(min_or_max == 1) expr += instance.getGeneralCoeffAlphaObjFollower(pr,s)*(alpha_max[s]/numeric);

            if(input.isNearOptMult() == true){
                if(min_or_max == 0) expr -= ((1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - instance.getModel()->follower_ub*(1.0 - b_alpha_min[s][it]));
                else if(min_or_max == 1) expr -= ((1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - instance.getModel()->follower_ub*(1.0 - b_alpha_max[s][it]));
            }else{
                if(min_or_max == 0) expr -= (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - (instance.getModel()->follower_ub + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))*(1.0 - b_alpha_min[s][it]));
                else if(min_or_max == 1) expr -= (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - (instance.getModel()->follower_ub + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))*(1.0 - b_alpha_max[s][it]));
            }
        }
    }
    else if(it < (instance.getGeneralNbConstrsAlpha(pr, s, min_or_max) + instance.getGeneralNbBdlbsAlpha(pr, s, min_or_max))){
        int i = instance.getGeneralBdlbsAlpha(pr, s, min_or_max, it - instance.getGeneralNbConstrsAlpha(pr, s, min_or_max));
        double dir = instance.getGeneralDirection(pr,s,i);
        BilevelVariable var = instance.getModel()->follower_vars[i];

        expr = 0;
        if(s == 0){
            if(min_or_max == 0) expr = (var.lb + (var.ub - var.lb)*(1.0 - b_alpha_min[s][it])) - (yi[i] + dir*(alpha_min[s]/numeric));
            else if(min_or_max == 1) expr = (var.lb + (var.ub - var.lb)*(1.0 - b_alpha_max[s][it])) - (yi[i] + dir*(alpha_max[s]/numeric));
        }else{
            if(min_or_max == 0) expr = (var.lb + (var.ub - var.lb)*(1.0 - b_alpha_min[s][it])) - (y[s-1][i] + dir*(alpha_min[s]/numeric));
            else if(min_or_max == 1) expr = (var.lb + (var.ub - var.lb)*(1.0 - b_alpha_max[s][it])) - (y[s-1][i] + dir*(alpha_max[s]/numeric));
        }
    }
    else if(it < (instance.getGeneralNbConstrsAlpha(pr, s, min_or_max) + instance.getGeneralNbBdlbsAlpha(pr, s, min_or_max) + instance.getGeneralNbBdubsAlpha(pr, s, min_or_max))){
        int i = instance.getGeneralBdlbsAlpha(pr, s, min_or_max, it - instance.getGeneralNbConstrsAlpha(pr, s, min_or_max) - instance.getGeneralNbBdlbsAlpha(pr, s, min_or_max));
        double dir = instance.getGeneralDirection(pr,s,i);
        BilevelVariable var = instance.getModel()->follower_vars[i];

        expr = 0;
        if(s == 0){
            if(min_or_max == 0) expr = (yi[i] + dir*(alpha_min[s]/numeric)) - (var.ub - (var.ub - var.lb)*(1.0 - b_alpha_min[s][it]));
            else if(min_or_max == 1) expr = (yi[i] + dir*(alpha_max[s]/numeric)) - (var.ub - (var.ub - var.lb)*(1.0 - b_alpha_max[s][it]));
        }else{
            if(min_or_max == 0) expr = (y[s-1][i] + dir*(alpha_min[s]/numeric)) - (var.ub - (var.ub - var.lb)*(1.0 - b_alpha_min[s][it]));
            else if(min_or_max == 1) expr = (y[s-1][i] + dir*(alpha_max[s]/numeric)) - (var.ub - (var.ub - var.lb)*(1.0 - b_alpha_max[s][it]));
        }
    }
}