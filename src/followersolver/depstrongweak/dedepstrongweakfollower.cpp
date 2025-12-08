#include "dedepstrongweakfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

void DEDepStrWkFollowerSolver::defineLeaderObj(){
    beta = leader->getGRBModel()->addVar(0.0,1.0,0.0,GRB_CONTINUOUS,"beta");
    obj  = leader->getGRBModel()->addVar(-GRB_INFINITY,GRB_INFINITY,1.0,GRB_CONTINUOUS,"obj");
    leader->getGRBModel()->addQConstr(beta*Fs + Fw - beta*Fw == obj);
    
    if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional){
        leader->getGRBModel()->addConstr(beta == (instance.getModel()->follower_ub - fs)/(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
    }else{
        // Computing interval values for PWL.
        int nb_intervals = input.getNbPWLIntervalsStrongWeak();
        double interval = (instance.getModel()->follower_ub - instance.getModel()->follower_lb)/(double)nb_intervals;
        
        double * follower_value = new double[nb_intervals+1];
        double * cooperation_level = new double[nb_intervals+1];

        follower_value[0] = instance.getModel()->follower_lb;
        cooperation_level[0] = instance.getStrongWeakProbab(follower_value[0]);
        for(int i = 1; i <= nb_intervals; ++i){
            follower_value[i] = follower_value[i-1] + interval;
            cooperation_level[i] = instance.getStrongWeakProbab(follower_value[i]);
        }

        // Add piecewise linear constraint beta = beta(follower_obj_value)
        leader->getGRBModel()->addGenConstrPWL(fs, beta, (nb_intervals + 1), follower_value, cooperation_level, "PWL_beta");

        delete[] follower_value;
        delete[] cooperation_level;
    }
}

void DEDepStrWkFollowerSolver::evaluate(double & mean, double & variance){
    double x_obj_leader = 0.0;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);

    double beta_ = instance.getStrongWeakProbab(fs_);
    
    mean = x_obj_leader + (beta_*Fs_ + (1.0-beta_)*Fw_);
    variance = 0.0;
}

void DEDepStrWkFollowerSolver::computeStrongWeakInteriorSolutions(){
    double beta_ = instance.getStrongWeakProbab(fs_);
    if(beta_ <= 0.0001) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,false,true);
    else if(beta_ >= 0.9999) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,true,true);
    else AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,false,true);
}