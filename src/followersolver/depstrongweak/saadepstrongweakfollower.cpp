#include "saadepstrongweakfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

void SAADepStrWkFollowerSolver::defineLeaderObj(){
    z = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_BINARY);

    for(int s = 0; s < instance.getNbScenarios(); ++s){
        double a_s = instance.getStrongWeakScenario(pr,s);
        double eps = input.getEpsBigM();
        leader->getGRBModel()->addConstr(fs <= a_s + (instance.getModel()->follower_ub-a_s)*(1.0 - z[s]));
        leader->getGRBModel()->addConstr(fs >= (a_s + eps) - (a_s + eps - instance.getModel()->follower_lb)*(z[s]));
    }

    for(int s = 0; s < instance.getNbScenarios(); ++s){
        double a_s = instance.getStrongWeakScenario(pr,s);
        double eps = input.getEpsBigM();
        leader->getGRBModel()->addGenConstrIndicator(z[s], 1, fs <= a_s);
        leader->getGRBModel()->addGenConstrIndicator(z[s], 0, fs >= a_s + eps);
    }

    for(int s = 0; s < instance.getNbScenarios()-1; ++s){
        leader->getGRBModel()->addConstr(z[s] >= z[s+1]);
    }

    GRBVar *delta = leader->getGRBModel()->addVars(instance.getNbScenarios(),GRB_CONTINUOUS);
    double bigm = instance.getModel()->leader_ub - instance.getModel()->leader_lb;
    for(int s = 0; s < instance.getNbScenarios(); ++s){
        leader->getGRBModel()->addConstr(delta[s] <= (Fw - Fs));
        leader->getGRBModel()->addConstr(delta[s] >= (Fw - Fs) - bigm*(1.0 - z[s]));
        leader->getGRBModel()->addConstr(delta[s] <= bigm*z[s]);

        leader->getGRBModel()->addGenConstrIndicator(z[s], 1, delta[s] == (Fw - Fs));
        leader->getGRBModel()->addGenConstrIndicator(z[s], 0, delta[s] == 0.0);
    }

    Fw.set(GRB_DoubleAttr_Obj,1.0);
    for(int s = 0; s < instance.getNbScenarios(); ++s) {
        delta[s].set(GRB_DoubleAttr_Obj,-(1.0/(double)instance.getNbScenarios()));
    }

}

double SAADepStrWkFollowerSolver::evaluate(){
    double beta_ = instance.getStrongWeakProbab(fs_);
    eval_avg = (beta_*Fs_ + (1.0-beta_)*Fw_);

    double eval = eval_avg;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        eval += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);
    return eval;
}

void SAADepStrWkFollowerSolver::computeStrongWeakInteriorSolutions(){
    if(z[0].get(GRB_DoubleAttr_X) <= 0.001) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,false,true);
    else if(z[instance.getNbScenarios()-1].get(GRB_DoubleAttr_X) >= 0.999) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,true,true);
    else AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,false,true);
}