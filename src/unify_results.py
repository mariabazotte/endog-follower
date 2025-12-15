

############################################################################################################
#                                              Instances
############################################################################################################

instances = []
seeds     = [i in range(1,11)]

############################################################################################################
#                                              Behaviors
############################################################################################################

# Follower near/optimality
near_opt     = [True]
eps_near_opt = [0.01, 0.05]
mul_near_opt = [True]

near_opt_cfg = []

# Follower behaviors
behaviors  = [0, 1, 2] 
spec_types = {0: [0], 1: [0,1,2,3], 2: [0,1,2]} 

name_behaviors  = ["fixstrongweak","depstrongweak","depgeneral"]
name_spec_types = { 0: ["-"], 1: ["proportional","threshold","strong","fragile"], 2: ["neutral","proportional","fragile"]}

# Parameters for each behavior type.
fix_strwk_coop = [1.0, 0.7, 0.5, 0.3, 0.0]
dep_strwk_thr  = [0.5, 2.0, 5.0, 10.0]
dep_strwk_str  = [0.5, 2.0, 5.0, 10.0]
dep_strwk_frg  = [0.5, 2.0, 5.0, 10.0]
dep_gen_nbint  = [4, 5]
dep_gen_scal   = [0.5, 2.0, 5.0, 10.0]
dep_gen_expl   = [True]

behavior_params = {
    0: {0: fix_strwk_coop},
    1: {0: [0.0], 1: dep_strwk_thr, 2: dep_strwk_str, 3: dep_strwk_frg},
    2: {0: [[e, 0, 0.0] for e in dep_gen_expl], 
        1: [[e, i, s] for e in dep_gen_expl for i in dep_gen_nbint for s in dep_gen_scal], 
        2: [[e, i, s] for e in dep_gen_expl for i in dep_gen_nbint for s in dep_gen_scal]}
}

# Complete configuration of behaviors.
behaviors_config = [(bv,tp,pr) for bv in behaviors for tp in spec_types[bv] for pr in behavior_params[bv][tp]]

def get_name_behavior_config(behavior_config):
    bv = behavior_config[0]
    tp = behavior_config[1]
    pr = behavior_config[2]
    return name_behaviors[bv] + "_" + name_spec_types[bv][tp] + "_" + str(params)

############################################################################################################
#                                              Scenarios
############################################################################################################

# SAA problem defined for behaviors 1 (Decision-dependent strong-weak) and 2 (Decision-dependent general).
# We solve behavior 0 (Fixed strong-weak) only with deterministic equivalent formulation.
nb_saa_problems  = {1: 5, 2: 5}
nb_saa_scenarios = {1: 500, 2: 500}
nb_saa_thinning  = {1: 1, 2: 5}

# Evaluation with sampling only defined for behavior 2 (Decision-dependent general).
# For behaviors 1 and 2 we can enumerate the endogenous scenarios and compute the 
# endogenous probability, while for behavior 3 we need to use MCMC (sampling from convex set).
nb_eval_problems  = 5
nb_eval_scenarios = 10000
nb_eval_thinning  = 50
coord_har = True
metrp_has = True








time = 


def gettstudent():
    if(len(seeds)==1):
        return 0
    if(len(seeds)==2):
        return 12.71
    if(len(seeds)==3):
        return 4.30
    if(len(seeds)==4):
        return 3.18
    if(len(seeds)==5):
        return 2.78
    if(len(seeds)==6):
        return 2.57
    if(len(seeds)==7):
        return 2.45
    if(len(seeds)==8):
        return 2.36
    if(len(seeds)==7):
        return 2.31
    if(len(seeds)==10):
        return 2.26
    return 0






folders_behavior = ["depgeneral/","depstrongweak/","fixstrongweak/"]
folvers_solution = ["dep/","saa/"]




def get_file_name():
    pass



############################################################################################################
#                                              Dictionaries
############################################################################################################

dict_sol_all = dict()

def init_dict_sol_all(bv,tp,params,instance,seed):
    dict_sol_all[bv]["instance"].append(instance) 
    dict_sol_all[bv]["seed"].append(seed)
    dict_sol_all[bv]["type"].append(name_spec_types[bv][tp])
    dict_sol_all[bv]["params"].append(params)

def add_dict_sol_all_dep(bv,tp,df):
    sol = "_dep"
    dict_sol_all[bv]["status"+sol].append(int(df["STATUS"][0]))
    dict_sol_all[bv]["lb"+sol].append(float(df["LB"][0]))
    dict_sol_all[bv]["ub"+sol].append(float(df["UB"][0]))
    dict_sol_all[bv]["gap"+sol].append(float(df["GAP"][0]))
    dict_sol_all[bv]["obj_follower"+sol].append(float(df["OBJ_FOLLOWER"][0]))
    
    if near_opt:
        dict_sol_all[bv]["obj_near_opt_follower"+sol].append(float(df["OBJ_OPT_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_pes_follower"+sol].append(float(df["OBJ_PES_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_mean_follower"+sol].append(float(df["OBJ_MEAN_FOLLOWER"][0]))
    
    dict_sol_all[bv]["nb_bnb"+sol].append(int(df["NB_BNB"][0]))
    dict_sol_all[bv]["time"+sol].append(float(df["TIME"][0]))

def add_dict_sol_all_saa(bv,tp,df):
    sol = "_saa"
    dict_sol_all[bv]["status"+sol].append(int(df["STATUS"][0]))
    
    dict_sol_all[bv]["lb"+sol].append(float(df["LB"][0]))
    dict_sol_all[bv]["std_lb"+sol].append(math.sqrt(float(df["VAR_LB"][0])))
    dict_sol_all[bv]["perc_std_lb"+sol].append(math.sqrt(float(df["VAR_LB"][0]))/math.abs(float(df["LB"][0])))
    dict_sol_all[bv]["stat_lb"+sol].append(float(df["STAT_LB"][0]))

    dict_sol_all[bv]["ub"+sol].append(float(df["UB"][0]))
    dict_sol_all[bv]["std_ub"+sol].append(math.sqrt(float(df["VAR_UB"][0])))
    dict_sol_all[bv]["perc_std_ub"+sol].append(math.sqrt(float(df["VAR_UB"][0]))/math.abs(float(df["UB"][0])))
    dict_sol_all[bv]["stat_ub"+sol].append(float(df["STAT_UB"][0]))

    dict_sol_all[bv]["gap"+sol].append(float(df["GAP"][0]))
    dict_sol_all[bv]["abs_gap"+sol].append(float(df["GAP"][0])*float(df["UB"][0]))
    dict_sol_all[bv]["std_gap"+sol].append(math.sqrt(float(df["VAR_GAP"][0])))
    dict_sol_all[bv]["perc_std_gap"+sol].append(math.sqrt(float(df["VAR_GAP"][0]))/(float(df["GAP"][0])*float(df["UB"][0])))
    dict_sol_all[bv]["stat_gap"+sol].append(float(df["STAT_GAP"][0]))

    dict_sol_all[bv]["obj_follower"+sol].append(float(df["OBJ_FOLLOWER"][0]))

    if near_opt:
        dict_sol_all[bv]["obj_near_opt_follower"+sol].append(float(df["OBJ_OPT_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_pes_follower"+sol].append(float(df["OBJ_PES_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_mean_follower"+sol].append(float(df["OBJ_MEAN_FOLLOWER"][0]))

    dict_sol_all[bv]["nb_bnb"+sol].append(int(df["NB_BNB"][0]))
    dict_sol_all[bv]["var_nb_bnb"+sol].append(float(df["VAR_NB_BNB"][0]))
    dict_sol_all[bv]["nb_pr_not_opt"+sol].append(int(df["NB_NOT_OPT"][0]))

    dict_sol_all[bv]["time"+sol].append(int(df["TIME"][0]))
    dict_sol_all[bv]["time_lb"+sol].append(int(df["TIME_LB"][0]))
    dict_sol_all[bv]["time_eval"+sol].append(int(df["TIME_EVAL"][0]))
    
def add_dict_sol_all_comp(bv):
    if bv != 1:
        raise ValueError("This function computes information available only for behavior 1 (Decision-dependent strong-weak).") 

    dict_sol_all[bv]["gap_lb"].append(100*(dict_sol_all[bv]["lb_dep"][-1] - dict_sol_all[bv]["lb_saa"][-1])/dict_sol_all[bv]["lb_dep"][-1])
    dict_sol_all[bv]["gap_ub"].append(100*(dict_sol_all[bv]["ub_dep"][-1] - dict_sol_all[bv]["ub_saa"][-1])/dict_sol_all[bv]["ub_dep"][-1])

dict_sol_mean = dict()

def init_dict_sol_mean(bv,tp,params,instance):
    dict_sol_mean[bv]["instance"].append(instance) 
    dict_sol_mean[bv]["type"].append(name_spec_types[bv][tp])
    dict_sol_mean[bv]["params"].append(params)

def add_dict_sol_mean(bv):
    if bv != 2: 
        for opt in ["gap_dep","nb_bnb_dep","time_dep"]:
            dict_sol_mean[bv][opt].append(np.mean(dict_sol_all[bv][opt][-len(seeds):]))
            dict_sol_mean[bv]["int"+opt].append(str(round(dict_sol_mean[bv][opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(dict_sol_all[bv][opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))
    if bv != 0:
        for opt in ["gap_saa","nb_bnb_saa","time_saa","time_lb_saa","time_eval_saa","nb_pr_not_opt","perc_std_lb_saa","perc_std_ub_saa","perc_std_gap_saa"]:
            dict_sol_mean[bv][opt].append(np.mean(dict_sol_all[bv][opt][-len(seeds):]))
            dict_sol_mean[bv]["int"+opt].append(str(round(dict_sol_mean[bv][opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(dict_sol_all[bv][opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))

        for opt in ["perc_std_lb_saa","perc_std_ub_saa","perc_std_gap_saa"]:
            dict_sol_mean[bv]["min"+opt].append(np.min(dict_sol_all[bv][opt][-len(seeds):]))
    if bv == 1:
        for opt in ["gap_lb","gap_ub"]:
            dict_sol_mean[bv][opt].append(np.mean(dict_sol_all[bv][opt][-len(seeds):]))
            dict_sol_mean[bv]["int"+opt].append(str(round(dict_sol_mean[bv][opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(dict_sol_all[bv][opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))

dict_sol_values_all = dict()

def add_dict_sol_values_all(instance,seed):
    dict_sol_values_all["instance"].append(instance)
    dict_sol_values_all["seed"].append(seed)

    nb_cfgs = len(behaviors_config)

    for cf in range(nb_cfgs):
        bv = behavior_config[cf][0]
        name = get_name_behavior_config(behavior_config[cf])
        if bv == 0:
            dict_sol_values_all["leader_obj_"+name].append(dict_sol_all[bv]["ub_dep"][-nb_cfgs]/dict_sol_all[bv]["ub_dep"][-nb_cfgs + cf])
        else:
            dict_sol_values_all["leader_obj_"+name].append(dict_sol_all[bv]["ub_dep"][-nb_cfgs]/dict_sol_all[bv]["ub_saa"][-nb_cfgs + cf])

    for cf in range(nb_cfgs):
        bv = behavior_config[cf][0]
        name = get_name_behavior_config(behavior_config[cf])
        if bv == 0:
            dict_sol_values_all["follower_obj_"+name].append(dict_sol_all[bv]["obj_follower_dep"][-nb_cfgs]/dict_sol_all[bv]["obj_follower_dep"][-nb_cfgs + cf])
        else:
            dict_sol_values_all["follower_obj_"+name].append(dict_sol_all[bv]["obj_follower_dep"][-nb_cfgs]/dict_sol_all[bv]["obj_follower_saa"][-nb_cfgs + cf])

    for cf in range(nb_cfgs):
        bv = behavior_config[cf][0]
        name = get_name_behavior_config(behavior_config[cf])
        if bv == 0:
            dict_sol_values_all["follower_obj_near_mean_"+name].append(dict_sol_all[bv]["obj_near_mean_follower_dep"][-nb_cfgs]/dict_sol_all[bv]["obj_near_mean_follower_dep"][-nb_cfgs + cf])
        else:
            dict_sol_values_all["follower_obj_near_mean_"+name].append(dict_sol_all[bv]["obj_near_mean_follower_dep"][-nb_cfgs]/dict_sol_all[bv]["obj_near_mean_follower_saa"][-nb_cfgs + cf])

dict_sol_values_mean = dict()

def add_dict_sol_values_mean(instance):
    dict_sol_values_mean["instance"].append(instance)

    nb_cfgs = len(behaviors_config)

    for cf in range(nb_cfgs):
        bv = behavior_config[cf][0]
        opt = "leader_obj_" + get_name_behavior_config(behavior_config[cf])

        dict_sol_values_mean[opt].append(np.mean(dict_sol_values_all[opt][-len(seeds):]))
        dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(dict_sol_values_all[opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))

    for cf in range(nb_cfgs):
        bv = behavior_config[cf][0]
        opt = "follower_obj_" + get_name_behavior_config(behavior_config[cf])

        dict_sol_values_mean[opt].append(np.mean(dict_sol_values_all[opt][-len(seeds):]))
        dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(dict_sol_values_all[opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))

    for cf in range(nb_cfgs):
        bv = behavior_config[cf][0]
        opt = "follower_obj_near_mean_" + get_name_behavior_config(behavior_config[cf])

        dict_sol_values_mean[opt].append(np.mean(dict_sol_values_all[opt][-len(seeds):]))
        dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(dict_sol_values_all[opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))



dict_comp_all = dict()

def add_dict_comp_all(bv_cfg,instance,seed,df):

    for idx, row in df.iterrows():
        bv = int(row["FOLLOWER_BEHAVIOR"])
        tp = 0 if row["COOPERATION_TYPE"] == "-" else int(row["COOPERATION_TYPE"])
        pr = [float(p) for p in row["PARAMS"].split("/")]



        column = name_behaviors[bv] + "_" + name_spec_types[tp] + "_" + str(params)

        dict_comp_all("leader_obj_mean_"+column).append(float(row["LEADER_OBJ_MEAN"]))
        dict_comp_all("leader_obj_std_"+column).append(float(row["LEADER_OBJ_VARIANCE"]))
        dict_comp_all("leader_obj_std_"+column).append(float(row["LEADER_VARS_ESS"]))
        dict_comp_all("leader_obj_std_"+column).append(float(row[LEADER_VARS_RHAT
        dict_comp_all("leader_obj_std_"+column).append(float(row[FOLLOWER_OBJ_OPT

        dict_comp_all("leader_obj_std_"+column).append(float(row[FOLLOWER_OBJ_NEAROPT_MEAN
        dict_comp_all("leader_obj_std_"+column).append(float(row[FOLLOWER_OBJ_NEAROPT_VARIANCE



dict_comp_mean = dict()











dict_comp_all = dict()
dict_comp_mean = dict()

def get_table_comp_overall_mean():
def get_table_()















average 

part 1) computational time + estimators value  -> like paper 1

compare leader and follower optimal values according to behavior -> optmistic is 100% and then the remaining are compared accordingly
generate tables
part 2)

compare optimal value -> you want to generate images

part 3)
computational time
