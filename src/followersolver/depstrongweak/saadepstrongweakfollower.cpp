#include "saadepstrongweakfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

void SAADepStrWkFollowerSolver::defineLeaderObj(){
    z = leader->getGRBModel()->addVars(instance.getStrongWeakNbScenarios(pr),GRB_BINARY);

    bool sequence = false;
    if(sequence == true){
        for(int s = 0; s < instance.getStrongWeakNbScenarios(pr); ++s){
            double a_s = instance.getStrongWeakScenario(pr,s);
            double eps = input.getEpsBigM();
            leader->getGRBModel()->addConstr(fs <= a_s + ((instance.getModel()->follower_ub + 0.001)-a_s)*(1.0 - z[s]));
            // leader->getGRBModel()->addConstr(fs >= (a_s + eps) - (a_s + eps - (instance.getModel()->follower_lb - 0.1))*(z[s]));
        }

        // for(int s = 0; s < instance.getStrongWeakNbScenarios(pr); ++s){
        //     double a_s = instance.getStrongWeakScenario(pr,s);
        //     double eps = input.getEpsBigM();
        //     leader->getGRBModel()->addGenConstrIndicator(z[s], 1, fs <= a_s);
        //     leader->getGRBModel()->addGenConstrIndicator(z[s], 0, fs >= a_s + eps);
        // }

        for(int s = 0; s < instance.getStrongWeakNbScenarios(pr)-1; ++s){
            leader->getGRBModel()->addConstr(z[s] >= z[s+1]);
        }

    }else{
        GRBLinExpr expr = 0;
        for(int s = 0; s < instance.getStrongWeakNbScenarios(pr); ++s) expr += z[s];
        leader->getGRBModel()->addConstr(expr <= 1.0);

        for(int s = 0; s < instance.getStrongWeakNbScenarios(pr); ++s){
            double a_s = instance.getStrongWeakScenario(pr,s);
            leader->getGRBModel()->addConstr(fs <= a_s + ((instance.getModel()->follower_ub + 0.001)-a_s)*(1.0 - z[s]));
            // leader->getGRBModel()->addGenConstrIndicator(z[s], 1, fs <= a_s);
        }
    }

    GRBVar *delta = leader->getGRBModel()->addVars(instance.getStrongWeakNbScenarios(pr),GRB_CONTINUOUS);
    double bigm = leader->getGeneralUB() - leader->getGeneralLB();
    // double bigm = instance.getModel()->leader_ub - instance.getModel()->leader_lb;
    
    for(int s = 0; s < instance.getStrongWeakNbScenarios(pr); ++s){
        leader->getGRBModel()->addConstr(delta[s] <= (Fw - Fs) + 0.0001);
        leader->getGRBModel()->addConstr(delta[s] >= (Fw - Fs) - bigm*(1.0 - z[s]));
        leader->getGRBModel()->addConstr(delta[s] <= bigm*z[s]);

        // leader->getGRBModel()->addGenConstrIndicator(z[s], 1, delta[s] == (Fw - Fs));
        // leader->getGRBModel()->addGenConstrIndicator(z[s], 0, delta[s] == 0.0);
    }

    Fw.set(GRB_DoubleAttr_Obj,1.0);
    for(int s = 0; s < instance.getStrongWeakNbScenarios(pr); ++s) {
        if(sequence == true) delta[s].set(GRB_DoubleAttr_Obj,-(1.0/(double)instance.getNbScenarios()));
        else delta[s].set(GRB_DoubleAttr_Obj,-((instance.getStrongWeakWeight(pr,s))/(double)instance.getNbScenarios()));
    }
}

void SAADepStrWkFollowerSolver::evaluate(double & mean, double & variance, double & f_mean, double & f_variance){
    double x_obj_leader = 0.0;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);

    double beta_ = instance.getStrongWeakProbab(fs_);
    
    mean = x_obj_leader + (beta_*Fs_ + (1.0-beta_)*Fw_);
    variance = 0.0;
    
    f_mean = fs_;
    if(input.isFollowerNearOpt() == true) 
        f_mean = (beta_*fs_eps_ + (1.0-beta_)*fw_eps_);
    f_variance = 0.0;
}

void SAADepStrWkFollowerSolver::computeStrongWeakInteriorSolutions(){
    if(z[0].get(GRB_DoubleAttr_X) <= 0.001) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,false,true);
    else if(z[instance.getStrongWeakNbScenarios(pr)-1].get(GRB_DoubleAttr_X) >= 0.999) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,true,true);
    else AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,false,true);
}