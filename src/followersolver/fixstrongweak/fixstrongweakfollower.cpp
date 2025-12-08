#include "fixstrongweakfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

void FixStrongWeakFollowerSolver::defineLeaderObj(){
    Fs.set(GRB_DoubleAttr_Obj,input.getFixCoopLevel());
    Fw.set(GRB_DoubleAttr_Obj,(1.0-input.getFixCoopLevel()));
}

void FixStrongWeakFollowerSolver::evaluate(double & mean, double & variance){  
    double x_obj_leader = 0.0;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        x_obj_leader += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);
    
    mean = x_obj_leader + (Fs_*input.getFixCoopLevel() + Fw_*(1.0-input.getFixCoopLevel()));
    variance = 0.0;
}

void FixStrongWeakFollowerSolver::computeStrongWeakInteriorSolutions(){
    if(input.getFixCoopLevel() <= 0.0001) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,false,true);
    else if(input.getFixCoopLevel() >= 0.9999) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,true,true);
    else AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,false,true);
}