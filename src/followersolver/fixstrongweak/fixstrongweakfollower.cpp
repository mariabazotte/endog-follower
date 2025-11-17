#include "fixstrongweakfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

void FixStrongWeakFollowerSolver::defineLeaderObj(){
    Fs.set(GRB_DoubleAttr_Obj,input.getFixCoopLevel());
    Fw.set(GRB_DoubleAttr_Obj,(1.0-input.getFixCoopLevel()));
}

double FixStrongWeakFollowerSolver::evaluate(){ 
    eval_avg = (Fs_*input.getFixCoopLevel() + Fw_*(1.0-input.getFixCoopLevel()));
    
    double eval = eval_avg;
    for(int i = 0; i < instance.getModel()->nb_leader_vars; ++i)
        eval += instance.getModel()->leader_vars[i].obj_leader*leader->getX_(i);
    return eval; 
}

void FixStrongWeakFollowerSolver::computeStrongWeakSolutions(){
    if(input.getFixCoopLevel() <= 0.0001) AbstractFollowerSolver::computeStrongWeakSolutions(true,false);
    if(input.getFixCoopLevel() >= 0.9999) AbstractFollowerSolver::computeStrongWeakSolutions(false,true);
}