#include "dedepstrongweakfollower.hpp"
#include "../../leadersolver/leadersolver.hpp"

JULIA_DEFINE_FAST_TLS

void DEDepStrWkFollowerSolver::defineLeaderObj(){
    beta = leader->getGRBModel()->addVar(-0.01,1.01,0.0,GRB_CONTINUOUS,"beta");
    obj  = leader->getGRBModel()->addVar(-GRB_INFINITY,GRB_INFINITY,1.0,GRB_CONTINUOUS,"obj");
    leader->getGRBModel()->addQConstr(beta*Fs + Fw - beta*Fw == obj,"bilinear_objective");
    
    if(input.getTypeDepStrongWeak() == Input::TypesDepStrongWeak::Proportional){
        leader->getGRBModel()->addConstr(beta == (instance.getModel()->follower_ub - fs)/(instance.getModel()->follower_ub - instance.getModel()->follower_lb));
    }else{
        computePointsPWL();
    }
}

void DEDepStrWkFollowerSolver::fixedPointsPWL() {
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

void DEDepStrWkFollowerSolver::computePointsPWL() {
    // Initializing julia.
    jl_init();              
    jl_eval_string("using PiecewiseLinApprox"); 

    // Loading modules.
    // jl_module_t* PWL_module = jl_get_module(jl_main_module, "PiecewiseLinApprox");
    jl_module_t* PWL_module = (jl_module_t*)jl_eval_string("PiecewiseLinApprox"); 
    if(!PWL_module) throw std::runtime_error("Could not get Julia modules.");

    // Loading functions.
    jl_function_t* absolute_func  = jl_get_function(PWL_module, "Absolute");
    // jl_function_t* relative_func  = jl_get_function(PWL_module, "Relative");
    jl_function_t* linearize_func = jl_get_function(PWL_module, "Linearize");
    if(!absolute_func || !linearize_func) throw std::runtime_error("Could not load julia functions.");

    // Calling relative/absolute error functions.
    jl_value_t* absolute_obj = jl_call1(absolute_func, jl_box_float64(input.absErrorPWLStrongWeak()));
    // jl_value_t* relative_obj = jl_call1(relative_func, jl_box_float64(input.relErrorPWLStrongWeak()));
    
    // Args for linearize function.
    jl_value_t** args;
    JL_GC_PUSHARGS(args, 4);
    args[0] = jl_eval_string((instance.getStrongWeakStringFunction()).c_str());
    args[1] = jl_box_float64(0.0);
    args[2] = jl_box_float64(1.0);
    args[3] = absolute_obj;
    // args[3] = relative_obj;
    
    // Calling linearize function.
    double time_ = time(NULL);
    jl_value_t* pwl_array = jl_call(linearize_func, args, 4);
    time_pwl = time(NULL) - time_;

    size_t nb_intervals = jl_array_len(pwl_array);
    jl_value_t** pwl_data = jl_array_data(pwl_array, jl_value_t*);

    double * follower_value = new double[nb_intervals+1];
    double * cooperation_level = new double[nb_intervals+1];

    follower_value[0] = instance.getModel()->follower_lb;
    cooperation_level[0] = instance.getStrongWeakProbab(follower_value[0]);
    for (size_t i = 0; i < nb_intervals; ++i) {
        jl_value_t* piece = pwl_data[i];
        jl_value_t* xMax_val = jl_get_field(piece, "xMax");
        double xMax = jl_unbox_float64(xMax_val);

        follower_value[i+1] = xMax*(instance.getModel()->follower_ub - instance.getModel()->follower_lb) + instance.getModel()->follower_lb;
        cooperation_level[i+1] = instance.getStrongWeakProbab(follower_value[i+1]);
    }

    // Add piecewise linear constraint beta = beta(follower_obj_value)
    leader->getGRBModel()->addGenConstrPWL(fs, beta, (nb_intervals + 1), follower_value, cooperation_level, "PWL_beta");

    delete[] follower_value;
    delete[] cooperation_level;

    jl_atexit_hook(0);
}

void DEDepStrWkFollowerSolver::evaluate(double & mean, double & variance, double & f_mean, double & f_variance){
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

void DEDepStrWkFollowerSolver::computeStrongWeakInteriorSolutions(){
    double beta_ = instance.getStrongWeakProbab(fs_);
    if(beta_ <= 0.0001) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(true,false,true);
    else if(beta_ >= 0.9999) AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,true,true);
    else AbstractFollowerSolver::computeStrongWeakInteriorSolutions(false,false,true);
}

void DEDepStrWkFollowerSolver::solve(){
    leader->getGRBModel()->set(GRB_DoubleParam_TimeLimit, std::max(0.0,input.getTimeLimit() - time_pwl));
    AbstractFollowerSolver::solve();
    leader->setTime(leader->getTime()+time_pwl);
}