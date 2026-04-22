#include "depgeneralfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

DepGeneralFollowerSolver::DepGeneralFollowerSolver(const Input & input, Instance & instance, LeaderSolver *leader, int pr) :
                                    AbstractFollowerSolver(input,instance,leader), pr(pr), metropolis_hastings(input.metropolisHastings()), 
                                    grb_minmax(input.useStepExplicitGRBMinMax()), big_m(false), numeric(1.0), cuts_mipnode(true) { 
    // Initialize pointers to NULL.
    y.resize(instance.getNbScenarios(),NULL);
    y_.resize(instance.getNbScenarios(),NULL);
    
    b_alpha_min.resize(instance.getNbScenarios(),NULL);
    b_alpha_max.resize(instance.getNbScenarios(),NULL);

    bound_alpha_min.resize(instance.getNbScenarios(),NULL);
    bound_alpha_max.resize(instance.getNbScenarios(),NULL);

    // alpha_ = new double[instance.getNbScenarios()];
    alpha_min_ = new double[instance.getNbScenarios()];
    alpha_max_ = new double[instance.getNbScenarios()];

    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty){
        b_alpha_min_ = new int[instance.getNbScenarios()];
        b_alpha_max_ = new int[instance.getNbScenarios()];
        std::fill(b_alpha_min_, b_alpha_min_+instance.getNbScenarios(), -1);
        std::fill(b_alpha_max_, b_alpha_max_+instance.getNbScenarios(), -1);

        bound_alpha_min_ = new double[instance.getNbScenarios()];
        bound_alpha_max_ = new double[instance.getNbScenarios()];
        
        for(int s = 0; s < instance.getNbScenarios(); ++s)
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
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        definePrimalVars(y[s],"_" + std::to_string(s+1));
        // definePrimalConstrs(y[s],"_" + std::to_string(s+1));
    }

    // Define step and step maximum and minimum value variables (alpha, alpha_min and alpha_max).
    // alpha = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_CONTINUOUS);
    alpha_min = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_CONTINUOUS);
    alpha_max = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_CONTINUOUS);
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        // alpha[s].set(GRB_DoubleAttr_LB,-instance.getStepSizeInterval()*numeric);
        // alpha[s].set(GRB_DoubleAttr_UB,instance.getStepSizeInterval()*numeric);
        // alpha[s].set(GRB_StringAttr_VarName,"alpha_"+std::to_string(s));

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
    // for(int s = 0; s < instance.getNbScenarios(); ++s){
    //     double tau = instance.getGeneralStep(pr,s);
    //     leader->getGRBModel()->addConstr(alpha[s] == (1.0 - tau)*alpha_min[s] + tau*alpha_max[s]);
    // }

    // Define maximum and minimum step variables (alpha_min and alpha_max).
    if((input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel) && 
       (input.useStepExplicitGeneral() == false)) defineKKTStep();  
    else defineExplicitStep();

    // Define sequence.
    defineSequence();

    if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
        // Define constraints for interval coefficients.
        defineIntCoeffConstrs();

        // Define sequence according to probability.
        if(metropolis_hastings == true) 
            defineSequenceWithMetropolisHastings();
        else{
            defineSequenceWithSliceSampling();
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

void DepGeneralFollowerSolver::defineIntCoeffConstrs(){
    // Define interval for probability distribution according to follower optimal value.
    GRBVar fsp = leader->getGRBModel()->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
    leader->getGRBModel()->addConstr(fsp == (instance.getModel()->follower_ub - fs) / (instance.getModel()->follower_ub - instance.getModel()->follower_lb)); 
    
    // Define w_k.
    GRBLinExpr sum_expr = 0;
    for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
        sum_expr += w[k];

        // leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, fsp >= input.getIntValueGeneral(k) + input.getEpsBigM());
        // leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, fsp <= input.getIntValueGeneral(k+1));

        leader->getGRBModel()->addConstr(fsp >= (input.getIntValueGeneral(k) + input.getEpsBigM())*w[k]);
        leader->getGRBModel()->addConstr(fsp <= 1.0 + (input.getIntValueGeneral(k+1) - 1.0)*w[k]);   
    }
    leader->getGRBModel()->addConstr(sum_expr == 1.0);
        
    // Define w_delta_k as multiplication of w_k*(Fw - Fs).
    double bigm = (leader->getGeneralUB() - leader->getGeneralLB());
    // double bigm = (instance.getModel()->leader_ub - instance.getModel()->leader_lb);

    for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
        // leader->getGRBModel()->addGenConstrIndicator(w[k], 0.0, w_delta[k] == 0.0);
        // leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, w_delta[k] == (Fw -Fs));

        leader->getGRBModel()->addConstr(w_delta[k] <= (Fw - Fs));
        leader->getGRBModel()->addConstr(w_delta[k] >= (Fw - Fs) - bigm*(1.0 - w[k]));
        leader->getGRBModel()->addConstr(w_delta[k] <= bigm*w[k]);
    }
}

// Function to define sequence of points y according to the chosen sampling method.
// The sequence starts with the chebyshev center defined by variables yi.
void DepGeneralFollowerSolver::defineSequence(){
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::Neutral){
        // Accept all new points.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                double tau = instance.getGeneralStep(pr,s);
                if(s == 0) leader->getGRBModel()->addConstr(y[s][i] == yi[i] + ((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric))*instance.getGeneralDirection(pr,s,i));
                else leader->getGRBModel()->addConstr(y[s][i] == y[s-1][i] + ((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric))*instance.getGeneralDirection(pr,s,i));
            }
        }
    }else{
        if(metropolis_hastings == false){
            // Accept all new points.
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    double tau = instance.getGeneralStep(pr,s);
                    if(s == 0) leader->getGRBModel()->addConstr(y[s][i] == yi[i] + ((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric))*instance.getGeneralDirection(pr,s,i));
                    else leader->getGRBModel()->addConstr(y[s][i] == y[s-1][i] + ((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric))*instance.getGeneralDirection(pr,s,i));
                }
            }
        }else{
            // Accept new points according to metropolis hastings rule.
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                    double tau = instance.getGeneralStep(pr,s);
                    if(s == 0){
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == yi[i] + ((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric))*instance.getGeneralDirection(pr,s,i));
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == yi[i]);
                    }else{
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, y[s][i] == y[s-1][i] + ((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric))*instance.getGeneralDirection(pr,s,i));
                        leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, y[s][i] == y[s-1][i]);
                    }
                }
            }
        }
    }
}

// Function to define acceptance variables for sequence in the case of metropolis hastings.
void DepGeneralFollowerSolver::defineSequenceWithMetropolisHastings(){
    // Define metropolis hastings acceptance variables according to distribution.
    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenProportional){
        // Define acceptance variables for proportional case.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            double r = instance.getGeneralAccep(pr,s);
            double tau = instance.getGeneralStep(pr,s);

            GRBLinExpr expr = 0;
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k)
                expr += (((r - 1.0) * input.getMaxIntCoeffGeneral()) / (2.0 * input.getIntCoeffGeneral(k)))*w_delta[k];

            for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
                BilevelVariable var = instance.getModel()->follower_vars[i];
                
                if(s == 0) expr += (r - 1.0)*var.obj_leader*yi[i];    
                else expr += (r - 1.0)*var.obj_leader*y[s-1][i];
            }
            expr -= (r - 1.0) * (Fw + Fs) / 2.0;
            expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric));

            leader->getGRBModel()->addGenConstrIndicator(q[s], 1.0, expr <= 0.0);
            leader->getGRBModel()->addGenConstrIndicator(q[s], 0.0, expr >= 0.0 + input.getEpsBigM());
        }
    }

    if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragile){
        // Define acceptance variables in fragile case.
        for(int s = 0; s < instance.getNbScenarios(); ++s){
            double r = instance.getGeneralAccep(pr,s);
            double tau = instance.getGeneralStep(pr,s);

            GRBLinExpr expr = 0;
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k)
                expr += (std::log(r) / input.getIntCoeffGeneral(k))*w_delta[k];

            expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric));

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
            double tau = instance.getGeneralStep(pr,s);

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
            expr -= instance.getGeneralCoeffAlphaObjLeader(pr,s)*((1.0 - tau)*(alpha_min[s]/numeric) + tau*(alpha_max[s]/numeric));
            
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
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_min[s]/numeric) >= expr);

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == instance.getStepSizeInterval());
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == instance.getStepSizeInterval());

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addConstr(b_alpha_max[s][index_max] <= 1.0 - w[k]);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_min], 1.0, (alpha_min[s]/numeric) <= expr);
                }

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k) <= (-1e-12))){
                    // Constraint defines alpha_max.   
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_max[s]/numeric) <= expr);

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == -instance.getStepSizeInterval());
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == -instance.getStepSizeInterval());

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addConstr(b_alpha_min[s][index_min] <= 1.0 - w[k]); 
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_max], 1.0, (alpha_max[s]/numeric) >= expr);
                }
            }
        }

        if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragile){
            GRBLinExpr expr_alpha_min = 0;
            GRBLinExpr expr_alpha_max = 0;
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                GRBLinExpr expr = (w_delta[k]*std::log(instance.getGeneralAccep(pr,s))) / 
                                  (input.getIntCoeffGeneral(k)*instance.getGeneralCoeffAlphaObjLeader(pr,s));

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) >= (1e-12)){
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_min[s]/numeric) >= expr);

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == instance.getStepSizeInterval());
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == instance.getStepSizeInterval());
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addConstr(b_alpha_max[s][index_max] <= 1.0 - w[k]);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_min], 1.0, (alpha_min[s]/numeric) <= expr); 
                    
                    expr_alpha_min += (w_delta[k]*std::log(instance.getGeneralAccep(pr,s))) / 
                                      (input.getIntCoeffGeneral(k)*instance.getGeneralCoeffAlphaObjLeader(pr,s));
                    expr_alpha_max += w[k]*instance.getStepSizeInterval();
                }

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k) <= (-1e-12))){
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_max[s]/numeric) <= expr);

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == -instance.getStepSizeInterval());

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == -instance.getStepSizeInterval());

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addConstr(b_alpha_min[s][index_min] <= 1.0 - w[k]);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_max], 1.0, (alpha_max[s]/numeric) >= expr);
                    
                    expr_alpha_min += w[k]*(-instance.getStepSizeInterval());
                    expr_alpha_max += (w_delta[k]*std::log(instance.getGeneralAccep(pr,s))) / 
                                      (input.getIntCoeffGeneral(k)*instance.getGeneralCoeffAlphaObjLeader(pr,s));
                }
            }

            // leader->getGRBModel()->addConstr((alpha_min[s]/numeric) >= expr_alpha_min);
            // leader->getGRBModel()->addConstr((alpha_max[s]/numeric) <= expr_alpha_max);

            // if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_min[s][index_min]/numeric) == expr_alpha_min);
            // if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_max[s][index_max]/numeric) == expr_alpha_max);
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
            double bigm = (leader->getGeneralUB() - leader->getGeneralLB())/2.0;
            // double bigm = (instance.getModel()->leader_ub - instance.getModel()->leader_lb)/2.0;
            
            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                leader->getGRBModel()->addConstr(w_s[k] <=  bigm*w[k]);
                leader->getGRBModel()->addConstr(w_s[k] >= -bigm*w[k]);
                leader->getGRBModel()->addConstr(w_s[k] <=  leader_obj_s - (Fw + Fs)/2.0 + bigm*(1.0 - w[k]));
                leader->getGRBModel()->addConstr(w_s[k] >=  leader_obj_s - (Fw + Fs)/2.0 - bigm*(1.0 - w[k]));

                // leader->getGRBModel()->addGenConstrIndicator(w[k], 0.0, w_s[k] == 0.0);
                // leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, w_s[k] == leader_obj_s - (Fw + Fs)/2.0);
            }

            GRBLinExpr expr_alpha_min = 0;
            GRBLinExpr expr_alpha_max = 0;

            for(int k = 0; k < input.getNbIntervalsGeneral(); ++k){
                double u = 0.0;
                if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenFragilePower)
                    u = (std::pow(instance.getGeneralAccep(pr,s), 1.0/input.getIntCoeffPowerFrgGeneral(k)) - 1.0);
                if(input.getTypeDepGeneral() == Input::TypesDepGeneral::GenStrongPower)
                    u = (std::pow(instance.getGeneralAccep(pr,s), 1.0/input.getIntCoeffPowerStrGeneral(k)) - 1.0);

                GRBLinExpr expr = (w_delta[k]*0.5*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) +
                                  (w_s[k]*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s));
                
                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) >= (1e-12)){
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_min[s]/numeric) >= expr);
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == instance.getStepSizeInterval());
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true)  leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == instance.getStepSizeInterval());

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addConstr(b_alpha_max[s][index_max] <= 1.0 - w[k]);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_min], 1.0, (alpha_min[s]/numeric) <= expr); 
                    
                    expr_alpha_min += (w_delta[k]*0.5*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) +
                                      (w_s[k]*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s));
                    expr_alpha_max += w[k]*instance.getStepSizeInterval();
                }

                if((instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k) <= (-1e-12))){
                    leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (alpha_max[s]/numeric) <= expr);
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == -instance.getStepSizeInterval());

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_max[s][index_max]/numeric) == expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrIndicator(w[k], 1.0, (bound_alpha_min[s][index_min]/numeric) == -instance.getStepSizeInterval());

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addConstr(b_alpha_min[s][index_min] <= 1.0 - w[k]);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_max], 1.0, (alpha_max[s]/numeric) >= expr);
                    
                    expr_alpha_min += w[k]*(-instance.getStepSizeInterval());
                    expr_alpha_max += (w_delta[k]*0.5*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s)*input.getIntCoeffSignalGeneral(k)) +
                                      (w_s[k]*u) / (instance.getGeneralCoeffAlphaObjLeader(pr,s));
                }  
            }

            // leader->getGRBModel()->addConstr((alpha_min[s]/numeric) >= expr_alpha_min);
            // leader->getGRBModel()->addConstr((alpha_max[s]/numeric) <= expr_alpha_max);

            // if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_min[s][index_min]/numeric) == expr_alpha_min);
            // if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_max[s][index_max]/numeric) == expr_alpha_max);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DepGeneralFollowerSolver::defineExplicitStep(){
    // Define constraints for alpha_min.
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        int nb_constrs_min = getNbConstrsMin(s);

        // Binary variable defining best lower bound on alpha_min.
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false) {
            b_alpha_min[s] = leader->getGRBModel()->addVars(nb_constrs_min,GRB_BINARY);

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

        // Continuous variable defining bounds.
        if((input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) ||
          (input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true)){
            bound_alpha_min[s] = leader->getGRBModel()->addVars(nb_constrs_min,GRB_CONTINUOUS);
            for(int i = 0; i < nb_constrs_min; ++i){
                bound_alpha_min[s][i].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
                bound_alpha_min[s][i].set(GRB_DoubleAttr_UB,0.0);
            }
        }

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
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);

                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr >= constr.rhs - 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr >= constr.rhs - (constr.rhs - constr.bound)*(1.0 - b_alpha_min[s][j])); 
                }
                if(constr.sense == '>'){
                    leader->getGRBModel()->addConstr(expr >= constr.rhs);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr <= constr.rhs + 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr <= constr.rhs + (constr.bound - constr.rhs)*(1.0 - b_alpha_min[s][j]));
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
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - instance.getModel()->follower_ub*(1.0 - b_alpha_min[s][j]));
                }else{
                    bound_expr += (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr(expr <= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][j], 1, expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - (instance.getModel()->follower_ub + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))*(1.0 - b_alpha_min[s][j]));
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
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_min[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_min[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_lb+j], 1, yi[i] + dir*(alpha_min[s]/numeric) <= var.lb + 1e-6);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_min[s][index_lb+j]));
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) >= var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_min[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_min[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_lb+j], 1, y[s-1][i] + dir*(alpha_min[s]/numeric) <= var.lb + 1e-6);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_min[s][index_lb+j]));
            }  
        }
        // Upper bound constraints defining alpha_min.
        for(int j = 0; j < instance.getGeneralNbBdubsAlpha(pr,s,0); ++j){
            int i = instance.getGeneralBdubsAlpha(pr,s,0,j);
            BilevelVariable var = instance.getModel()->follower_vars[i];
            
            double dir = instance.getGeneralDirection(pr,s,i);
            int index_ub = instance.getGeneralNbConstrsAlpha(pr,s,0) + instance.getGeneralNbBdlbsAlpha(pr,s,0);
            if(s == 0){
                leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) <= var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_min[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_min[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_ub+j], 1, yi[i] + dir*(alpha_min[s]/numeric) >= var.ub - 1e-6);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_min[s]/numeric) >= var.ub - (var.ub-var.lb)*(1.0 - b_alpha_min[s][index_ub+j]));
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) <= var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_min[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_min[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_min[s][index_ub+j], 1, y[s-1][i] + dir*(alpha_min[s]/numeric) >= var.ub - 1e-6);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_min[s]/numeric) >= var.ub - (var.ub - var.lb)*(1.0 - b_alpha_min[s][index_ub+j]));
            }
        }
        // Lower bound on alpha_min.
        int index = instance.getGeneralNbConstrsAlpha(pr,s,0) + 
                    instance.getGeneralNbBdlbsAlpha(pr,s,0) + 
                    instance.getGeneralNbBdubsAlpha(pr,s,0);
        leader->getGRBModel()->addConstr((alpha_min[s]/numeric) >= -instance.getStepSizeInterval());
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel) leader->getGRBModel()->addConstr((alpha_min[s]/numeric) <= -instance.getStepSizeInterval()*b_alpha_min[s][index]);
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_min[s][index]/numeric) == -instance.getStepSizeInterval());
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_min[s][index]/numeric) == -instance.getStepSizeInterval());

        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrMax(alpha_min[s], bound_alpha_min[s], nb_constrs_min, -instance.getStepSizeInterval()*numeric);
    }

    // Define constraints for alpha_max.
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        int nb_constrs_max = getNbConstrsMax(s);

        // Binary variable defining best upper bound on alpha_max.
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel){
            b_alpha_max[s] =  leader->getGRBModel()->addVars(nb_constrs_max,GRB_BINARY);

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

        // Continuous variable defining bounds on alpha_max.
        if((input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) ||
           (input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true)){
            bound_alpha_max[s] = leader->getGRBModel()->addVars(nb_constrs_max,GRB_CONTINUOUS);
            for(int i = 0; i < nb_constrs_max; ++i){
                bound_alpha_max[s][i].set(GRB_DoubleAttr_LB,0.0);
                bound_alpha_max[s][i].set(GRB_DoubleAttr_UB,GRB_INFINITY);
            }
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
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr >= constr.rhs - 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr >= constr.rhs - (constr.rhs - constr.bound)*(1.0 - b_alpha_max[s][j]));
                }
                if(constr.sense == '>'){
                    leader->getGRBModel()->addConstr(expr >= constr.rhs);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr <= constr.rhs + 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr <= constr.rhs + (constr.bound - constr.rhs)*(1.0 - b_alpha_max[s][j]));
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
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr >= (1.0 - input.getEpsNearOpt())*fs + input.getEpsNearOpt()*instance.getModel()->follower_ub - instance.getModel()->follower_ub*(1.0 - b_alpha_max[s][j]));
                }else{
                    bound_expr += (fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))/instance.getGeneralCoeffAlphaObjFollower(pr,s);
                    leader->getGRBModel()->addConstr(expr <= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][j]/numeric) == bound_expr);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][j], 1, expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - 1e-6);
                    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(expr >= fs + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) - (instance.getModel()->follower_ub + input.getEpsNearOpt()*(instance.getModel()->follower_ub - instance.getModel()->follower_lb))*(1.0 - b_alpha_max[s][j]));
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
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_max[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_max[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_lb+j], 1, yi[i] + dir*(alpha_max[s]/numeric) <= var.lb + 1e-6);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_lb+j]));                
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) >= var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_max[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_max[s][index_lb+j]/numeric) == var.lb);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_lb+j], 1, y[s-1][i] + dir*(alpha_max[s]/numeric) <= var.lb + 1e-6);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) <= var.lb + (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_lb+j]));
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
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_max[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(yi[i] + dir*(bound_alpha_max[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_ub+j], 1, yi[i] + dir*(alpha_max[s]/numeric) >= var.ub - 1e-6);    
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(yi[i] + dir*(alpha_max[s]/numeric) >= var.ub - (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_ub+j]));
            }else{
                leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) <= var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_max[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(bound_alpha_max[s][index_ub+j]/numeric) == var.ub);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == false) leader->getGRBModel()->addGenConstrIndicator(b_alpha_max[s][index_ub+j], 1, y[s-1][i] + dir*(alpha_max[s]/numeric) >= var.ub - 1e-6);
                if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == false && big_m == true) leader->getGRBModel()->addConstr(y[s-1][i] + dir*(alpha_max[s]/numeric) >= var.ub - (var.ub - var.lb)*(1.0 - b_alpha_max[s][index_ub+j]));
            }
        }
        // Upper bound on alpha_max.
        int index = instance.getGeneralNbConstrsAlpha(pr,s,1) + 
                    instance.getGeneralNbBdlbsAlpha(pr,s,1) + 
                    instance.getGeneralNbBdubsAlpha(pr,s,1);
        leader->getGRBModel()->addConstr((alpha_max[s]/numeric) <= instance.getStepSizeInterval());
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel) leader->getGRBModel()->addConstr((alpha_max[s]/numeric) >= instance.getStepSizeInterval()*b_alpha_max[s][index]);
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty) leader->getGRBModel()->addConstr((bound_alpha_max[s][index]/numeric) == instance.getStepSizeInterval());
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addConstr((bound_alpha_max[s][index]/numeric) == instance.getStepSizeInterval());

        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel && grb_minmax == true) leader->getGRBModel()->addGenConstrMin(alpha_max[s], bound_alpha_max[s], nb_constrs_max, instance.getStepSizeInterval());
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
    if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == true)){ eval = Fs_; var_eval = 0.0; f_eval = fs_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
    else if((std::abs(Fw_ - Fs_) <= 1e-3) && (input.isFollowerNearOpt() == false) && (std::abs(fw_eps_ - fs_eps_) <= 1e-3)){ eval = Fs_; var_eval = 0.0; f_eval = fs_eps_; f_var_eval = 0.0; ess = 1.0; rhat = 1.0; }
    else instance.evaluateDepGeneral(eval,var_eval,f_eval,f_var_eval,ess,rhat,leader->getX_(),yi_,Fs_,Fw_,fs_,input.getTypeDepGeneral(),input.getNbIntervalsGeneral(),input.getScalingGeneral());

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
    if(leader->getNbSol() > 0){
        if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty){
            for(int s = 0; s < instance.getNbScenarios(); ++s){
                if(y_[s]) delete[] y_[s];
                y_[s] = leader->getGRBModel()->get(GRB_DoubleAttr_X,y[s],instance.getModel()->nb_follower_vars);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DepGeneralFollowerSolver::solve(){
    if(input.getMethodDepGeneral() == Input::MethodDepGeneral::SingleLevel){
        leader->getGRBModel()->set(GRB_IntParam_Presolve, 1);
        AbstractFollowerSolver::solve();
    }else if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Relax){
        leader->getGRBModel()->set(GRB_IntParam_Presolve, 1);
        solve_relax();
        if(leader->getNbSol() > 0) computeStrongWeakInteriorSolutions();
    }else if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Penalty){
        solve_penalty();
        if(leader->getNbSol() > 0) computeStrongWeakInteriorSolutions();
    }else if(input.getMethodDepGeneral() == Input::MethodDepGeneral::LocalSearch){
        solve_local_search();
        if(leader->getNbSol() > 0) computeStrongWeakInteriorSolutions();
    }else if(input.getMethodDepGeneral() == Input::MethodDepGeneral::Callback){
        solve_callback();
        if(leader->getNbSol() > 0) computeStrongWeakInteriorSolutions();
    }
}

void DepGeneralFollowerSolver::solve_relax(){
    leader->getGRBModel()->set(GRB_IntParam_Presolve, 0);
    
    double start_time = time(NULL); 

    // Solving relaxed version of the problem
    leader->getGRBModel()->update();
    leader->getGRBModel()->optimize();

    leader->setStatus(statusFromGurobi(leader->getGRBModel()->get(GRB_IntAttr_Status)));
    if(leader->getStatus() == Status::Optimal || leader->getStatus() == Status::Time_Limit) upd_solution();
    else leader->verify_stat(); 

    // Getting feasible solution.
    if(leader->getNbSol() > 0){
        double x_obj_leader = 0.0;
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
            x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);

        double y_obj_leader = 0.0;
        instance.evaluateSAADepGeneral(y_obj_leader,pr,leader->getX_(),yi_,Fs_,Fw_,fs_,
               input.getTypeDepGeneral(),input.getNbIntervalsGeneral(),input.getScalingGeneral());
        
        leader->setUB((x_obj_leader+y_obj_leader));
    }

    leader->setTime(time(NULL) - start_time);
    if(input.getVerbose() >= 1) leader->print();

}

void DepGeneralFollowerSolver::solve_penalty(){
    double start_time = time(NULL); 
    leader->getGRBModel()->set(GRB_IntParam_OutputFlag, 0);

    double best_ub = std::numeric_limits<double>::infinity();
    double best_fs_ = std::numeric_limits<double>::infinity();;
    double * best_x_ = new double[instance.getModel()->nb_leader_vars];
    double * best_yi_ = new double[instance.getModel()->nb_follower_vars];
     
    double rho = rho_init;
    double curr_obj = 0.0;

    // Initialize problem.

    // Solve first subproblem. 
    leader->getGRBModel()->update();
    leader->getGRBModel()->optimize();

    Status sub_status = statusFromGurobi(leader->getGRBModel()->get(GRB_IntAttr_Status));
    if(sub_status == Status::Optimal || sub_status == Status::Time_Limit) {
        upd_solution();
        curr_obj = leader->getGRBModel()->get(GRB_DoubleAttr_ObjVal);
    }else leader->verify_stat();

    // Solve second subproblem.
    solve_sub_penalty(rho);

    // Upper bound.
    evaluateSAA(best_ub,best_fs_,best_x_,best_yi_);

    // Iterations for penalty.
    for(int it = 0; it < rho_max_it; ++it){
        // Iterations for partial minima.
        bool stop = false;
        while(stop == false){
            // Setting time limit.
            leader->getGRBModel()->set(GRB_DoubleParam_TimeLimit, std::max(0.0, input.getTimeLimit() - (time(NULL) - start_time)));

            // Solve first subproblem. 
            leader->getGRBModel()->update();
            leader->getGRBModel()->optimize();

            // Get first subproblem solution.
            double inf_norm = 0.0;
            double prev_obj = curr_obj;
            Status sub_status = statusFromGurobi(leader->getGRBModel()->get(GRB_IntAttr_Status));
            if(sub_status == Status::Optimal || sub_status == Status::Time_Limit) { 
                get_variables(inf_norm);
                curr_obj = leader->getGRBModel()->get(GRB_DoubleAttr_ObjVal);
            }else leader->verify_stat();

            if(sub_status == Status::Time_Limit) leader->setStatus(Status::Time_Limit);

            // Solve second subproblem. (dual of alpha problems).
            solve_sub_penalty(rho);

            // Compute upper bound.
            evaluateSAA(best_ub,best_fs_,best_x_,best_yi_);

            // Stop with partial minima solution or with time limit.
            if((inf_norm <= tol_part_min) || (std::abs((curr_obj - prev_obj)) <= tol_part_min_obj) || (leader->getStatus() == Status::Time_Limit)) stop = true; 

            if(input.getVerbose() >= 1) printf ("i %3d | %12.4f | %12.4f | %12.4f | %12.4f \n", it, inf_norm, std::abs((curr_obj - prev_obj)),  best_ub, time(NULL)-start_time);
        }
        if(leader->getStatus() == Status::Time_Limit) break;

        // Compute duality error.
        double duality_error = computeDualityError();
        if(input.getVerbose() >= 1) printf ("rho %12.4f | %12.4f | %12.4f | %12.4f \n", rho, duality_error, best_ub, time(NULL)-start_time);
        if(duality_error <= tol_bi_feas) {
            // Use status optimal meaining it converged within time limit and max nb iterations.
            leader->setStatus(Status::Optimal);
            break;
        }

        // Update penalty parameter.
        rho *= rho_mult;
    }
    leader->setTime(time(NULL) - start_time);

    // Updating best solution.
    leader->setUB(best_ub);
    leader->setX_(best_x_);

    fs_ = best_fs_;
    std::copy(best_yi_, best_yi_ + instance.getModel()->nb_follower_vars, yi_);

    if(best_x_) delete[] best_x_;
    if(best_yi_) delete[] best_yi_;

    if(leader->getStatus() == Status::Optimal){
        bool converged_sol_best = true;
        if(best_ub <= leader->getGRBModel()->get(GRB_DoubleAttr_ObjVal) - 1e-6) 
            converged_sol_best = false;
        else converged_sol_best = true;
    }

    // Verify if stopped by number of iterations.
    if(leader->getStatus() == Status::Not_Solved) leader->setStatus(Status::Iteration_Limit);

    if(input.getVerbose() >= 1) leader->print();
}

void DepGeneralFollowerSolver::solve_sub_penalty(double rho){
    // Solve second subproblem. (dual of alpha problems).
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        int index_bd_min, index_bd_max;
        double alpha_min_new, alpha_max_new;
        if(s == 0){
            alpha_min_new = instance.defineMinStep(index_bd_min, pr, s, leader->getX_(), yi_, fs_, Fs_, Fw_, input.getTypeDepGeneral(), input.getNbIntervalsGeneral(), input.getScalingGeneral(), true);
            alpha_max_new = instance.defineMaxStep(index_bd_max, pr, s, leader->getX_(), yi_, fs_, Fs_, Fw_, input.getTypeDepGeneral(), input.getNbIntervalsGeneral(), input.getScalingGeneral(), true);
        }else{
            alpha_min_new = instance.defineMinStep(index_bd_min, pr, s, leader->getX_(), y_[s-1], fs_, Fs_, Fw_, input.getTypeDepGeneral(), input.getNbIntervalsGeneral(), input.getScalingGeneral(), true);
            alpha_max_new = instance.defineMaxStep(index_bd_max, pr, s, leader->getX_(), y_[s-1], fs_, Fs_, Fw_, input.getTypeDepGeneral(), input.getNbIntervalsGeneral(), input.getScalingGeneral(), true);
        }

        // if(b_alpha_min_[s] != index_bd_min) {
        //     std::cout << s << " Min: " << b_alpha_min_[s] << " " << index_bd_min << std::endl;
        //     std::cout << s << " Alpha min: " << alpha_min_[s] << " " << alpha_min_new << std::endl; 
        // }
        // if(b_alpha_max_[s] != index_bd_max) {
        //     std::cout << s << " Max: " << b_alpha_max_[s] << " " << index_bd_max << std::endl;
        //     std::cout << s << " Alpha max: " << alpha_max_[s] << " " << alpha_max_new << std::endl; 
        // }

        b_alpha_min_[s] = index_bd_min;
        b_alpha_max_[s] = index_bd_max;

        // alpha_min_[s] = alpha_min_new;
        // alpha_max_[s] = alpha_max_new;

        for(int i = 0; i < getNbConstrsMin(s); ++i) 
            bound_alpha_min[s][i].set(GRB_DoubleAttr_Obj, 0.0);
        bound_alpha_min[s][index_bd_min].set(GRB_DoubleAttr_Obj, -rho);

        for(int i = 0; i < getNbConstrsMax(s); ++i) 
            bound_alpha_max[s][i].set(GRB_DoubleAttr_Obj, 0.0);
        bound_alpha_max[s][index_bd_max].set(GRB_DoubleAttr_Obj, rho);   
    }
}

void DepGeneralFollowerSolver::get_variables(double & inf_norm){
    // Get variables and inf. norm for leader variables.
    leader->get_variables(inf_norm);

    // Copy previous follower solution for comparison.
    double * prev_yi_ = new double[instance.getModel()->nb_follower_vars];
    std::copy(yi_, yi_ + instance.getModel()->nb_follower_vars, prev_yi_);
    
    std::vector<double*> prev_y_(instance.getNbScenarios(),NULL);
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        prev_y_[s] = new double[instance.getModel()->nb_follower_vars];
        std::copy(y_[s], y_[s] + instance.getModel()->nb_follower_vars, prev_y_[s]);
    }

    // Update current solution.
    if(yi_) { delete[] yi_; yi_ = NULL; }
    yi_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,yi,instance.getModel()->nb_follower_vars);
    
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        if(y_[s]) delete[] y_[s];
        y_[s] = leader->getGRBModel()->get(GRB_DoubleAttr_X,y[s],instance.getModel()->nb_follower_vars);
    }

    fs_ = fs.get(GRB_DoubleAttr_X);
    if(input.getTypeDepGeneral() != Input::TypesDepGeneral::Neutral){
        Fs_ = Fs.get(GRB_DoubleAttr_X);
        Fw_ = Fw.get(GRB_DoubleAttr_X);
    }  

    // Compute infinity norm of the difference between current and previous solutions.
    inf_norm = std::max(inf_norm, inf_norm_diff(prev_yi_, yi_, instance.getModel()->nb_follower_vars));
    for(int s = 0; s < instance.getNbScenarios(); ++s)
        inf_norm = std::max(inf_norm, inf_norm_diff(prev_y_[s], y_[s], instance.getModel()->nb_follower_vars));

    // Delete previous solution
    delete[] prev_yi_;
    for(int s = 0; s < instance.getNbScenarios(); ++s) delete[] prev_y_[s];
}

double DepGeneralFollowerSolver::computeDualityError(){
    if(alpha_min_) delete[] alpha_min_;
    if(alpha_max_) delete[] alpha_max_;
    alpha_min_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,alpha_min,instance.getNbScenarios());
    alpha_max_ = leader->getGRBModel()->get(GRB_DoubleAttr_X,alpha_max,instance.getNbScenarios());

    for(int s = 0; s < instance.getNbScenarios(); ++s){
        bound_alpha_min_[s] = bound_alpha_min[s][b_alpha_min_[s]].get(GRB_DoubleAttr_X);
        bound_alpha_max_[s] = bound_alpha_max[s][b_alpha_max_[s]].get(GRB_DoubleAttr_X);
    }

    double strong_duality_error = 0.0;
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        strong_duality_error = std::max(strong_duality_error, std::abs(alpha_min_[s] - bound_alpha_min_[s])/std::max(1e-12,std::abs(alpha_min_[s])));
        strong_duality_error = std::max(strong_duality_error, std::abs(bound_alpha_max_[s] - alpha_max_[s])/std::max(1e-12,std::abs(alpha_max_[s])));
    }
    return strong_duality_error;
}

void DepGeneralFollowerSolver::evaluateSAA(double & best_ub_, double & best_fs_, double * best_x_, double * best_yi_){
    double x_obj_leader = 0.0;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);

    double y_obj_leader = 0.0;
    instance.evaluateSAADepGeneral(y_obj_leader,pr,leader->getX_(),yi_,Fs_,Fw_,fs_,
        input.getTypeDepGeneral(),input.getNbIntervalsGeneral(),input.getScalingGeneral());
    
    if((x_obj_leader + y_obj_leader) <= best_ub_){
        best_ub_ = (x_obj_leader + y_obj_leader);

        best_fs_ = fs_;
        std::copy(leader->getX_(), leader->getX_() + instance.getModel()->nb_leader_vars, best_x_);
        std::copy(yi_, yi_ + instance.getModel()->nb_follower_vars, best_yi_);
    }
}

void DepGeneralFollowerSolver::solve_local_search(){
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
        BilevelVariable var = instance.getModel()->leader_vars[i];
        max_delta = std::min(max_delta,(int)((var.ub - var.lb)/2));
    }
    num_pairs = 5*instance.getModel()->nb_leader_vars;

    std::mt19937 rng(7919);
    std::uniform_int_distribution<int> sign_dist(0,1);
    std::uniform_int_distribution<int> var_dist(0,instance.getModel()->nb_leader_vars-1); 

    /////////////////////////////////////////////////////////////////////

    double start_time = time(NULL); 

    // Auxiliars for best solution.
    double best_ub_ = std::numeric_limits<double>::infinity();
    double best_fs_ = std::numeric_limits<double>::infinity();;  
    double * best_x_ = new double[instance.getModel()->nb_leader_vars];
    double * best_yi_ = new double[instance.getModel()->nb_follower_vars];

    // Auxiliar for current solution.
    double * test_x_ = new double[instance.getModel()->nb_leader_vars];

    // Initialize solution. 
    // Solution for optimistic follower with interior solution defined.
    leader->getGRBModel()->update(); 
    leader->getGRBModel()->optimize();

    Status sub_status = statusFromGurobi(leader->getGRBModel()->get(GRB_IntAttr_Status));
    if(sub_status == Status::Optimal || sub_status == Status::Time_Limit) upd_solution();
    else leader->verify_stat();
    evaluateSAA(best_ub_, best_fs_, best_x_, best_yi_);

    // Loop for local search.
    int it = 0;
    int delta = 1;
    while(delta <= max_delta){
        ++it;
        double prev_ub = best_ub_;
        
        if(input.getVerbose() >= 1) 
            printf ("delta: %3d | it: %3d | ub: %12.4f | time: %12.4f \n", delta, it, best_ub_, time(NULL)-start_time);

        // Including and removing a delta constant at each leader variable.
        std::copy(best_x_, best_x_ + instance.getModel()->nb_leader_vars, test_x_);
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
            BilevelVariable var = instance.getModel()->leader_vars[i];
            
            if((test_x_[i] + delta) <= var.ub){
                test_x_[i] += delta;
                if(verify_leader_feas(test_x_) == true){
                    bool bilevel_feasible = upd_follower_sol(test_x_);
                    if((bilevel_feasible == true) && (Fs_ <= best_ub_ - 1e-6)) {
                        evaluateSAA(best_ub_,best_fs_,best_x_,best_yi_);
                    }
                }
                test_x_[i] -= delta;
                if((time(NULL) - start_time) >= input.getTimeLimit()) break;
            }

            if((test_x_[i] - delta) >= var.lb){
                test_x_[i] -= delta;
                if(verify_leader_feas(test_x_) == true){
                    bool bilevel_feasible = upd_follower_sol(test_x_);
                    if((bilevel_feasible == true) && (Fs_ <= best_ub_ - 1e-6)) {
                        evaluateSAA(best_ub_,best_fs_,best_x_,best_yi_);
                    }
                }
                test_x_[i] += delta;
                if((time(NULL) - start_time) >= input.getTimeLimit()) break;
            }
        }

        // Sampling pairs to try to improve decision.
        for(int p = 0; p < num_pairs; ++p){
            int sign_i = (sign_dist(rng) == 0) ? -1 : 1;
            int sign_j = (sign_dist(rng) == 0) ? -1 : 1;

            int i = var_dist(rng);
            int j = var_dist(rng);
            while(i == j) j = var_dist(rng);

            const BilevelVariable & var_i = instance.getModel()->leader_vars[i];
            const BilevelVariable & var_j = instance.getModel()->leader_vars[j];

            if((test_x_[i] + sign_i*delta >= var_i.lb) && (test_x_[i] + sign_i*delta <= var_i.ub) &&
                (test_x_[j] + sign_j*delta >= var_j.lb) && (test_x_[j] + sign_j*delta <= var_j.ub)){

                test_x_[i] += sign_i * delta;
                test_x_[j] += sign_j * delta;

                if(verify_leader_feas(test_x_) == true){
                    bool bilevel_feasible = upd_follower_sol(test_x_);
                    if((bilevel_feasible == true) && (Fs_ <= best_ub_ - 1e-6)){
                        evaluateSAA(best_ub_,best_fs_,best_x_,best_yi_);
                    }
                }
                test_x_[i] -= sign_i * delta;
                test_x_[j] -= sign_j * delta;
            }
        }

        // Sampling pairs to try to improve decision.
        // for(int p = 0; p < num_pairs; ++p){
        //     int sign_i = (sign_dist(rng) == 0) ? -1 : 1;
        //     int sign_j = (sign_dist(rng) == 0) ? -1 : 1;  
        //     std::vector<int> feasible_i;
        //     std::vector<int> feasible_j;
        //     for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
        //         const BilevelVariable & var_i = instance.getModel()->leader_vars[i];
        //         if((test_x_[i] + sign_i*delta >= var_i.lb) && (test_x_[i] + sign_i*delta <= var_i.ub)) 
        //             feasible_i.push_back(i);
        //     }
        //     for(int j = 0; j < instance.getModel()->nb_leader_vars; ++j){
        //         const BilevelVariable & var_j = instance.getModel()->leader_vars[j];
        //         if((test_x_[j] + sign_j*delta >= var_j.lb) && (test_x_[j] + sign_j*delta <= var_j.ub)) 
        //             feasible_j.push_back(j);
        //     }
        //     if(feasible_i.empty() || feasible_j.empty()) continue;
        //     std::uniform_int_distribution<int> i_dist(0, feasible_i.size() - 1);
        //     std::uniform_int_distribution<int> j_dist(0, feasible_j.size() - 1);
        //     int i = feasible_i[i_dist(rng)];
        //     int j = feasible_j[j_dist(rng)];
        //     if (i == j) continue;
        //     test_x_[i] += sign_i * delta;
        //     test_x_[j] += sign_j * delta;
        //     if(verify_leader_feas(test_x_) == true){
        //         bool bilevel_feasible = upd_follower_sol(test_x_);
        //         if(bilevel_feasible){
        //             evaluateSAA(best_ub_,best_fs_,best_x_,best_yi_);
        //         }
        //     }
        //     test_x_[i] -= sign_i * delta;
        //     test_x_[j] -= sign_j * delta;
        // }

        if(best_ub_ >= prev_ub - 1e-6){
            delta += 1;
            if(delta > max_delta) leader->setStatus(Status::Optimal);
        }else{ 
            delta = 1; 
        }
        if(it >= nb_max_it) {
            leader->setStatus(Status::Iteration_Limit);
        }
        if((time(NULL) - start_time) >= input.getTimeLimit()){
            leader->setStatus(Status::Time_Limit);
        }
        if(leader->getStatus() == Status::Iteration_Limit || leader->getStatus() == Status::Time_Limit) break;
    }
    leader->setTime(time(NULL) - start_time);
    
    // Updating best solution.
    leader->setUB(best_ub_);
    leader->setX_(best_x_);

    fs_ = best_fs_;
    std::copy(best_yi_, best_yi_ + instance.getModel()->nb_follower_vars, yi_);

    if(best_x_) delete[] best_x_;
    if(best_yi_) delete[] best_yi_;
    if(test_x_) delete[] test_x_;

    if(input.getVerbose() >= 1) leader->print();
}

bool DepGeneralFollowerSolver::verify_leader_feas(double * test_x_){
    if(instance.getModel()->nb_leader_constrs > 0){
        for(int c = 0; c < instance.getModel()->nb_leader_constrs; ++c){
            BilevelConstraint constr = instance.getModel()->leader_constrs[c];
            
            double sum_x = 0.0;
            for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
                BilevelVariable var = instance.getModel()->leader_vars[i];
                sum_x += constr.coeffs[var.id]*test_x_[i];
            }

            if(constr.sense == '=') { if(sum_x >= constr.rhs + 1e-6 || sum_x <= constr.rhs - 1e-6) return false; }
            else if(constr.sense == '<'){ if(sum_x >= constr.rhs + 1e-6) return false; } 
            else if(constr.sense == '>'){ if(sum_x <= constr.rhs - 1e-6) return false; }
        }
    }
    return true;
}

bool DepGeneralFollowerSolver::upd_follower_sol(double * test_x_){
    // Update leader decision with
    leader->setX_(test_x_);

    // Solve follower problem.
    bool bilevel_feasible = solveFollowerProblem();

    // Solve interior, optimistic and pessimistic solutions.
    if(bilevel_feasible == true) {
        AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,true,true);
    }

    return bilevel_feasible;
}

void DepGeneralFollowerSolver::defineFollowerProblem(){
    model_follower = new GRBModel(*leader->getGRBEnv());

    // if(input.getVerbose() <= 1) model_follower->set(GRB_IntParam_OutputFlag, 0);
    // else model_follower->set(GRB_IntParam_OutputFlag, 1);
    model_follower->set(GRB_IntParam_OutputFlag, 0);
    model_follower->set(GRB_IntParam_Threads, input.getNbThreads());

    // Follower variables.
    y_follower = model_follower->addVars(instance.getModel()->nb_follower_vars,GRB_CONTINUOUS);
    for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
        BilevelVariable var = instance.getModel()->follower_vars[i];
        
        if(var.lb <= -instance.getModel()->inf) 
            y_follower[i].set(GRB_DoubleAttr_LB,-GRB_INFINITY);
        else y_follower[i].set(GRB_DoubleAttr_LB,var.lb);
        if(var.ub >=  instance.getModel()->inf) 
            y_follower[i].set(GRB_DoubleAttr_UB,GRB_INFINITY);
        else y_follower[i].set(GRB_DoubleAttr_UB,var.ub);

        y_follower[i].set(GRB_StringAttr_VarName,var.name);
        y_follower[i].set(GRB_CharAttr_VType,var.type);

        // Define follower objective.
        y_follower[i].set(GRB_DoubleAttr_Obj,var.obj_follower);
    }

    // Follower constraints.
    for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
        BilevelConstraint constr = instance.getModel()->follower_constrs[c];

        GRBLinExpr expr = 0;
        for(int i = 0; i < instance.getModel()->nb_follower_vars; ++i){
            BilevelVariable var = instance.getModel()->follower_vars[i];
            expr += constr.coeffs[var.id]*y_follower[i];
        }

        if(constr.sense == '=') 
            follower_constrs.push_back(model_follower->addConstr(expr == constr.rhs, constr.name));
        else if(constr.sense == '<'){
            follower_constrs.push_back(model_follower->addConstr(expr <= constr.rhs, constr.name));
        }else if(constr.sense == '>'){
            follower_constrs.push_back(model_follower->addConstr(expr >= constr.rhs, constr.name));
        }
    }
}

bool DepGeneralFollowerSolver::solveFollowerProblem(){
    // Update follower constraints according to current x_.
    for(int c = 0; c < instance.getModel()->nb_follower_constrs; ++c){
        const BilevelConstraint & constr = instance.getModel()->follower_constrs[c];

        double sum_x = 0.0;
        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i){
            BilevelVariable var = instance.getModel()->leader_vars[i];
            sum_x += constr.coeffs.at(var.id)*leader->getX_(i); 
        }
        follower_constrs[c].set(GRB_DoubleAttr_RHS, (constr.rhs - sum_x));
    }

    // Solve optimal follower and get solution.
    model_follower->update();
    model_follower->optimize();

    Status status_follower = statusFromGurobi(model_follower->get(GRB_IntAttr_Status));
    if(status_follower == Status::Optimal){   
        if(ys_) { delete[] ys_; ys_ = NULL; }
        ys_ = model_follower->get(GRB_DoubleAttr_X,y_follower,instance.getModel()->nb_follower_vars);
        fs_ = model_follower->get(GRB_DoubleAttr_ObjVal);
        return true;
    }
    else{ return false; } 
}

void DepGeneralFollowerSolver::defineCallback(){
    leader->defineBinaryVariables();

    theta = leader->getGRBModel()->addVar(leader->getGeneralLB(),GRB_INFINITY,1.0,GRB_CONTINUOUS);
    leader->getGRBModel()->addConstr(theta >= Fs);
    leader->getGRBModel()->addConstr(theta >= leader->getGeneralLB());

    leader->getGRBModel()->set(GRB_IntParam_LazyConstraints, 1);
}

void DepGeneralFollowerSolver::solve_callback(){
    double start_time = time(NULL);
    
    leader->getGRBModel()->setCallback(this);
    leader->getGRBModel()->update();
    leader->getGRBModel()->optimize();

    leader->setStatus(statusFromGurobi(leader->getGRBModel()->get(GRB_IntAttr_Status)));
    if(leader->getStatus() == Status::Optimal || leader->getStatus() == Status::Time_Limit) upd_solution();
    else leader->verify_stat(); 
    if(leader->getNbSol() > 0) {
        yi_ = new double[instance.getModel()->nb_follower_vars];
        AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,false,true); // yi_
    }

    leader->setTime(time(NULL) - start_time);
    if(input.getVerbose() >= 1) leader->print();
}

void DepGeneralFollowerSolver::callback(){
    try{
        if (where == GRB_CB_MIPSOL){ 
            // Getting current integer solution.
            double * x_ = getSolution(leader->getX(), instance.getModel()->nb_leader_vars);
            double * bx_ = getSolution(leader->getXBinary(), leader->getTotalNbXBinary());
            double theta_ = getSolution(theta);
            leader->setX_(x_);
            fs_ = getSolution(fs);
            if(ys_) { delete[] ys_; ys_ = NULL; }
            ys_ = getSolution(ys, instance.getModel()->nb_follower_vars);
            AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,true,true); // Fs_, Fw_, yi_ (Fs_ not guaranteed optimal in the current solution, should compute it here)

            // Compute correct expected value.
            double x_obj_leader = 0.0;
            for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
                x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*x_[i];

            double y_obj_leader = 0.0;
            instance.evaluateSAADepGeneral(y_obj_leader,pr,x_,yi_,Fs_,Fw_,fs_,
                input.getTypeDepGeneral(),input.getNbIntervalsGeneral(),input.getScalingGeneral());
            
            // Include cut to obtain correct expected value for current leader solution.
            if(y_obj_leader >= theta_ + 1e-6){
                int size = 0;
                GRBLinExpr expr = 0;
                for(int j = 0; j < leader->getTotalNbXBinary(); ++j){
                    if(bx_[j] >= 0.999) { expr += (y_obj_leader - leader->getGeneralLB())*leader->getXBinary(j); ++size; }
                    else if(bx_[j] <= 0.001) { expr -= (y_obj_leader - leader->getGeneralLB())*leader->getXBinary(j); }
                }
                expr -= (y_obj_leader - leader->getGeneralLB())*size;
                expr += y_obj_leader;
                addLazy(theta >= expr);
            }

            // New best solution found, pass it to the solver.
            if((x_obj_leader + y_obj_leader) <= getDoubleInfo(GRB_CB_MIPSOL_OBJBST) - 1e-6){
                setSolution(leader->getX(),x_,instance.getModel()->nb_leader_vars);
                setSolution(leader->getXBinary(),bx_,leader->getTotalNbXBinary());
                setSolution(theta,y_obj_leader);
                setSolution(Fs,Fs_);
                setSolution(fs,fs_);
                setSolution(ys,ys_,instance.getModel()->nb_follower_vars);
                if(input.isFollowerNearOpt() == true){
                    setSolution(ys_eps,ys_eps_,instance.getModel()->nb_follower_vars);
                    setSolution(fs_eps,fs_eps_);
                }
            }

            // Free memory.
            if(x_) delete[] x_;
            if(bx_) delete[] bx_;
        }
        if (where == GRB_CB_MIPNODE && cuts_mipnode == true){
            if(getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL){
                // Get current x_ and bx_.
                double * x_ = getNodeRel(leader->getX(), instance.getModel()->nb_leader_vars);
                double * bx_ = getNodeRel(leader->getXBinary(),leader->getTotalNbXBinary());
                
                // Verify if leader solution is integer.
                bool all_integer = true;
                for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i) {
                    if(std::abs(x_[i] - std::round(x_[i])) >= 1e-6) {
                        all_integer = false;
                        break;
                    }
                }
                if(all_integer == true){
                    for(int i = 0; i < leader->getTotalNbXBinary(); ++i){
                        if(std::abs(bx_[i] - std::round(bx_[i])) >= 1e-6) {
                            all_integer = false;
                            break;
                        }
                    }
                }

                if(all_integer == true){
                    double theta_ = getNodeRel(theta);
                    leader->setX_(x_);
                    bool optimal = solveFollowerProblem(); // fs_ (not well defined at mipnode)
                    if(optimal == true){
                        AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,true,true); // Fs_, Fw_, yi_ 

                        // Compute correct expected value.
                        double x_obj_leader = 0.0;
                        for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
                            x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*x_[i];

                        double y_obj_leader = 0.0;
                        instance.evaluateSAADepGeneral(y_obj_leader,pr,x_,yi_,Fs_,Fw_,fs_,
                            input.getTypeDepGeneral(),input.getNbIntervalsGeneral(),input.getScalingGeneral());
                        
                        // Include cut to obtain correct expected value for current leader solution.
                        if(y_obj_leader >= theta_ + 1e-6){
                            int size = 0;
                            GRBLinExpr expr = 0;
                            for(int j = 0; j < leader->getTotalNbXBinary(); ++j){
                                if(bx_[j] >= 0.999) { expr += (y_obj_leader - leader->getGeneralLB())*leader->getXBinary(j); ++size; }
                                else if(bx_[j] <= 0.001) { expr -= (y_obj_leader - leader->getGeneralLB())*leader->getXBinary(j); }
                            }
                            expr -= (y_obj_leader - leader->getGeneralLB())*size;
                            expr += y_obj_leader;
                            addLazy(theta >= expr);
                        }

                        if((x_obj_leader + y_obj_leader) <= getDoubleInfo(GRB_CB_MIPNODE_OBJBST) - 1e-6){
                            setSolution(leader->getX(),x_,instance.getModel()->nb_leader_vars);
                            setSolution(leader->getXBinary(),bx_,leader->getTotalNbXBinary());
                            setSolution(theta,y_obj_leader);
                            setSolution(Fs,Fs_);
                            setSolution(fs,fs_);
                            setSolution(ys,ys_,instance.getModel()->nb_follower_vars);
                            if(input.isFollowerNearOpt() == true){
                                setSolution(ys_eps,ys_eps_,instance.getModel()->nb_follower_vars);
                                setSolution(fs_eps,fs_eps_);
                            }
                            // setSolution(Fw,Fw_);
                            // setSolution(yw,yw_,instance.getModel()->nb_follower_vars);
                        }
                    }
                }
                
                // Free memory.
                delete[] x_;
                delete[] bx_;
            }
        }
    }
    catch (GRBException e){
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (const std::runtime_error& e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Standard exception: " << e.what() << std::endl;
    }
    catch (...){
        std::cout << "Error during callback" << std::endl;
    }
}
