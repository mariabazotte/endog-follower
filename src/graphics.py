
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

############################################################################################################
#                                              Behaviors
############################################################################################################

# Follower near/optimality
near_opt     = 1
eps_near_opt = 0.01
mul_near_opt = 0

# Alignment between leader and follower objectives information.
plot_align = True
align_types = 3
align_ranges = [[-1.0,-0.1],[-0.1,0.1],[0.1,1.0]]

# Follower behaviors
behaviors  = [0, 1, 2] 
spec_types = {0: [0], 1: [0,1,2,3,4,5], 2: [0,1,2,3,4]} 

name_behaviors  = ["fixstrongweak","depstrongweak","depgeneral"]
name_spec_types = { 0: [""], 
                    1: ["proportional","threshold","strong","fragile","strong_power","fragile_power"], 
                    2: ["neutral","proportional","fragile","fragile_power","strong_power"]}

# Parameters for each behavior type.
fix_strwk_coop  = [1.0, 0.7, 0.5, 0.3, 0.0]
dep_strwk_thr = [2.0, 5.0, 10.0]
dep_strwk_str = [2.0, 5.0, 10.0]
dep_strwk_frg = [2.0, 5.0, 10.0]
dep_strwk_str_pow = [2.0, 5.0, 10.0]
dep_strwk_frg_pow = [2.0, 5.0, 10.0]
dep_gen_nbint  = [-1.0]
dep_gen_scal_frg = [10.0, 20.0, 40.0, 80.0]
dep_gen_scal_frg_pow = [10.0, 20.0]
dep_gen_scal_str_pow = [10.0, 20.0]

behavior_params = {
    0: {0: fix_strwk_coop},
    1: {0: [0.0], 1: dep_strwk_thr, 2: dep_strwk_str, 3: dep_strwk_frg, 4: dep_strwk_str_pow, 5: dep_strwk_frg_pow},
    2: {0: [[0.0, 0.0]], 1: [[i, 1.0] for i in dep_gen_nbint], 2: [[i, s] for i in dep_gen_nbint for s in dep_gen_scal_frg],3: [[i, s] for i in dep_gen_nbint for s in dep_gen_scal_frg_pow], 4: [[i, s] for i in dep_gen_nbint for s in dep_gen_scal_str_pow]}
}

# Complete configuration of behaviors.
behaviors_config = [(bv,tp,pr) for bv in behaviors for tp in spec_types[bv] for pr in behavior_params[bv][tp]]

def get_name_behavior_config(bv_config):
    bv = bv_config[0]
    tp = bv_config[1]
    pr = bv_config[2]
    return name_behaviors[bv] + "_" + name_spec_types[bv][tp] + "_" + str(pr)

name_dep_strwk = ['Prp','Thr','Str','Frg','Str-p','Frg-p']
name_dep_gen   = ['Ntr','Prp','Frg','Frg','Str']
def get_name_behavior_config_graphic(bv_config):
    bv = bv_config[0]
    tp = bv_config[1]
    pr = bv_config[2]

    if bv == 0:
        val = float(pr)
        pr_format = f"{int(val)}" if val.is_integer() else f"{val}"
        return rf'$\beta^{{{pr_format}}}$'
    if bv == 1:
        val = float(pr)
        pr_format = f"{int(val)}" if val.is_integer() else f"{val}"
        if tp == 0:
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}(\cdot)$'
        if tp == 1:
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}^{{{pr_format}}}(\cdot)$' #, $\delta = {pr}$'
        if tp == 2 or tp == 3: 
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}^{{{pr_format}}}(\cdot)$' #, $a = {pr}$'
        if tp == 4 or tp == 5: 
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}^{{{pr_format}}}(\cdot)$' #, $p = {pr}$'
    if bv == 2:
        val = float(pr[1])
        pr_format = f"{int(val)}" if val.is_integer() else f"{val}"
        if tp == 0 or tp == 1:
            return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}}}(\cdot)$'
        if tp == 2 or tp == 4:
            return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}}}^{{{pr_format}}}(\cdot)$' #, $K = {pr[0]}$, $\gamma = {pr[1]}$'
        if tp == 3:
            return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}_{{+}}}}^{{{pr_format}}}(\cdot)$' # $K = {pr[0]}$, $\gamma = {pr[1]}$'

label_behaviors = [get_name_behavior_config_graphic(bv_cfg) for bv_cfg in behaviors_config]

def get_name_behavior_config_general_graphic(bv,tp):
    if bv == 0:
        return rf'$\beta^{{{tp}}}$'
    if bv == 1:
        return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}(\cdot)$'
    if bv == 2:
        if tp == 3:
            return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}_{{+}}}}(\cdot)$'
        else:
            return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}}}(\cdot)$'

############################################################################################################
#                                              Directory
############################################################################################################

# Define file and save paths.
name_folder = "../results/"
# name_folder = "../results-thesis/"

file_path = name_folder + "compiled/results_nearopt" + str(near_opt) 
if near_opt:
    file_path += "_eps" + str(eps_near_opt) + "_mult" + str(mul_near_opt)
file_path += '.xlsx'

name_common = "nearopt" + str(near_opt)
if near_opt:
    name_common += "_eps" + str(eps_near_opt)

save_path = name_folder + "compiled/graphics/"
sol_values_path = save_path + 'sol_values/'
cmp_values_path = save_path + 'cmp_values/'
perf_path = save_path + 'perf_swd/'
perf_path_bv = save_path + 'perf_all_bv/'

############################################################################################################
#                                              Graphics
############################################################################################################

colors = {'leader_obj_': '#0072B2','follower_obj_': '#009E73','follower_obj_near_mean_': '#E69F00'}
styles = {'leader_obj_':'--','follower_obj_':'-.','follower_obj_near_mean_':':'}
hatch_styles = {'leader_obj_':'/','follower_obj_':'x','follower_obj_near_mean_':'o'}

def plot_sol_metric(y_column, y_column_name, align):
    """
    Plots the average and quantiles of a behavior comparison gap metric (y_column) 
    by follower behavior, and saves the figure.
    
    Parameters:
    - y_column: Name of the column on the excel file to use for the y-axis. ('leader_obj_','follower_obj_','follower_obj_near_mean_')
    - y_column_name: Name to display on the y-axis of the figure.
    - align: [-1,0,1,2,3] -1 for all instances 
    """

    # Set global font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        'text.usetex': True,
        'font.size': 15,
        'axes.titlesize': 14,
        'axes.labelsize': 15,
        'xtick.labelsize': 15,
        'ytick.labelsize': 15,
        'legend.fontsize': 15,
        'axes.linewidth': 1,
        'lines.markersize': 3
    })

    if y_column not in ['leader_obj_','follower_obj_','follower_obj_near_mean_']:
        raise ValueError(f"Invalid y_column: {y_column}")

    # Columns correspnding to metric y_column
    columns = [y_column+get_name_behavior_config(bv_cfg) for bv_cfg in behaviors_config]

    # Reading sheet.
    df = pd.read_excel(file_path, sheet_name="sol_values_all")

    if align >= 0 and align <= align_types:
        df = df[(df["align"] >= align_ranges[align][0]) & (df["align"] <= align_ranges[align][1] - 1e-12)]

    df = df[columns]

    # Ignore the whole line (instance) if at least one configuration has "-".
    # df = df[~(df == "-").any(axis=1)]
    # df = df.apply(pd.to_numeric)

    # Ignore only the specific instance configuration pair with "-".
    df = df.replace("-", np.nan).infer_objects(copy=False)
    df = df.apply(pd.to_numeric)

    df = df.clip(lower=-100,upper=100)

    # print(df.mean())       # line plot values
    # print(df.describe())   # boxplot min, 25%, 50%, 75%, max

    ###############################################################

    y = df.mean()
    y_std = df.std()
    x = range(0,20*len(y),20)

    n = df.count()
    y_se = y_std / np.sqrt(n)
    ci = 1.96 * y_se  # 95% confidence interval

    plt.figure(figsize=(11, 4.8))

    plt.plot(x, y, label="Mean", color=colors[y_column], linestyle=styles[y_column], marker='o', linewidth=1)
    plt.fill_between(x, y - ci, y + ci, color=colors[y_column], alpha=0.25)
    # plt.errorbar(x, y, yerr=ci, fmt='o', color=colors[y_column], label="Mean", capsize=3)

    # plt.xlabel(r'Follower behavior')
    # plt.ylabel(f'{y_column_name}')

    plt.xticks(x, label_behaviors, rotation=65)
    plt.margins(x=0.01)

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(sol_values_path + name_common + "_align" + str(align) + f"_{y_column}" + "line_mean.png", dpi=300)
    plt.close()

    ###############################################################

    plt.figure(figsize=(12, 4.8))

    plt.plot(x, df.median(), label="Median", color=colors[y_column], linestyle=styles[y_column], marker='o', linewidth=1)
    plt.fill_between(x, df.quantile(0.25), df.quantile(0.75), color=colors[y_column], alpha=0.25)

    plt.xlabel(r'Follower behavior')
    # plt.ylabel(f'{y_column_name}')
    
    plt.xticks(x, label_behaviors, rotation=65)
    plt.margins(x=0.01)

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # plt.savefig(sol_values_path + name_common + "_align" + str(align) + f"_{y_column}" + "line_median.png", dpi=300)
    plt.close()

    ############################################################################

    plt.figure(figsize=(12, 4.8))

    sns.boxplot(data=df, color=colors[y_column], showfliers=False)

    plt.xticks(ticks=range(len(columns)), labels=label_behaviors, rotation=65)
    plt.grid(True)
    plt.tight_layout()

    # for i, tick in enumerate(plt.gca().get_xticklabels()):
    #     if i < 5:
    #         tick.set_color('red')
    #     elif i < 21:
    #         tick.set_color('blue')
    #     else:
    #         tick.set_color('green')

    plt.savefig(sol_values_path + name_common + "_align" + str(align) + f"_{y_column}" + "boxplot.png", dpi=300)
    plt.close()

    ############################################################################

behavior_params_acc = {
    1: {0: 0.0, 1: 10.0, 2: 10.0, 3: 10.0, 4: 10.0, 5: 10.0},
    2: {0: [0.0, 0.0], 1: [-1.0, 1.0], 2: [-1.0, 80.0],3: [-1.0, 20.0], 4: [-1.0, 20.0]}
}
spec_types_plot_acc = {0: fix_strwk_coop, 1: spec_types[1], 2: spec_types[2]}
behavior_spec_colors_acc = {0: {0.0:"#CC79A7", 0.3:"#F0E442", 0.5:"#999999", 0.7:"#56B4E9", 1.0:"#D55E00"},
                        1: {0:"#882255", 1:"#44AA99", 2:"#DDCC77", 3:"#AA4499", 4:"#117733", 5:"#332288"},
                        2: {0:"#661100", 1:"#6699CC", 2:"#AA4466", 3:"#4477AA", 4:"#228833"}}
behavior_spec_styles_acc = {0: {0.0:"-", 0.3:"--", 0.5:"-.", 0.7:":", 1.0:(0, (5, 2))}, 
                        1: {0:"-", 1:"--", 2:"-.", 3:":", 4:(0, (5, 2)), 5:(0, (3, 1, 1, 1))},
                        2: {0:"-", 1:"--", 2:"-.", 3:":", 4:(0, (5, 2))}}

# Accumulative version: selected subtypes for each behavior.
def plot_sol_metric_acc_per_bv(y_column, y_column_name, bv, align):
    """
    Plots the accumulated gap metric by percentage of instance for the different subtypes of
    follower behavior {bv} for instances with respective alignment {align}, and saves the figures. 
    
    Parameters:
    - y_column: Name of the column on the excel file to use for the y-axis. ('leader_obj_','follower_obj_','follower_obj_near_mean_')
    - y_column_name: Name to display on the y-axis of the figure.
    - bv in [0,1,2]: Behavior general type.
    - align: [-1,0,1,2,3] -1 for all instances.
    """ 

    # Set font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        'text.usetex': True,
        'font.size': 7.5,
        'axes.titlesize': 8,
        'axes.labelsize': 7.5,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'axes.linewidth': 0.5,
        'lines.linewidth': 1.0,
        'lines.markersize': 2.0,
    }) 

    if y_column not in ['leader_obj_','follower_obj_','follower_obj_near_mean_']:
        raise ValueError(f"Invalid y_column: {y_column}")

    # Reading sheet.
    df = pd.read_excel(file_path, sheet_name="sol_values_all")

    if align >= 0 and align <= 3:
        df = df[(df["align"] >= align_ranges[align][0]) & (df["align"] <= align_ranges[align][1] - 1e-12)]

    ###############################################################

    plt.figure(figsize=(3.4, 2.4))

    # Types
    for tp in spec_types_plot_acc[bv]:
        # Names.
        column = None
        column_label = None
        
        if bv == 0:
            column = y_column+get_name_behavior_config((bv,0,tp))
            column_label = get_name_behavior_config_graphic((bv,0,tp))

        if bv == 1 or bv == 2:
            column = y_column+get_name_behavior_config((bv,tp,behavior_params_acc[bv][tp]))
            column_label = get_name_behavior_config_graphic((bv,tp,behavior_params_acc[bv][tp]))

        # Get dataframe specific column.
        df_copy = df.copy()
        df_copy[column] = pd.to_numeric(df_copy[column].replace("-", np.nan), errors="coerce")
        df_copy["align"] = pd.to_numeric(df_copy["align"].replace("-", np.nan), errors="coerce")

        df_copy = df_copy.dropna(subset=[column, "align"])
        df_copy[column] = df_copy[column].clip(lower=-100, upper=100)

        # if align == -1:
        #     y = df_copy.sort_values(by="align")[column]
        # else:
        y = df_copy.sort_values(by=column)[column]

        x = [i/y.shape[0]*100 for i in range(1,y.shape[0]+1)]

        plt.plot(x, y, label=column_label, color=behavior_spec_colors_acc[bv][tp], linestyle=behavior_spec_styles_acc[bv][tp], linewidth=1)

    plt.xlabel(r'Instances (\%)')
    plt.ylabel(f'{y_column_name}')
    
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(sol_values_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + "_bv" + str(bv) +  "_align" + str(align) + "_accumulated_" + f"{y_column}.png", dpi=300)
    plt.close()

    ############################################################################

colors_cmp = {'gap_leader_': '#0072B2','gap_follower_': '#009E73','gap_follower_near_': '#E69F00'}
styles_cmp = {'gap_leader_':'--','gap_follower_':'-.','gap_follower_near_':':'}
hatch_styles_cmp = {'gap_leader_':'/','gap_follower_':'x','gap_follower_near_':'o'}
def plot_cmp_metric(bv_cfg, y_column, y_column_name, align):
    """
    Plots the average and quantiles of a metric (y_column) by follower behavior, and saves the figure.
    
    Parameters:
    - y_column: Name of the column to use for the y-axis. ('gap_leader_','gap_follower_','gap_follower_near_')
    - y_column_name: Name to display on the y-axis 
    - align: [-1,0,1,2,3] -1 for all instances
    """

    # Set global font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        'text.usetex': True,
        'font.size': 15,
        'axes.titlesize': 14,
        'axes.labelsize': 15,
        'xtick.labelsize': 15,
        'ytick.labelsize': 15,
        'legend.fontsize': 15,
        'axes.linewidth': 1,
        'lines.markersize': 3
    })

    # Columns corresponding to metric y_column
    columns = [y_column+get_name_behavior_config(column_bv_cfg) for column_bv_cfg in behaviors_config]
    rows = get_name_behavior_config(bv_cfg)

    # Reading sheet.
    df = pd.read_excel(file_path, sheet_name="comp_values_all")

    if align >= 0 and align <= 3:
        df = df[(df["align"] >= align_ranges[align][0]) & (df["align"] <= align_ranges[align][1] - 1e-12)]

    df = df[df["config"] == rows]
    df = df[columns]

    # Ignore the whole line (instance) if at least one configuration has "-".
    # df = df[~(df == "-").any(axis=1)]
    # df = df.apply(pd.to_numeric)

    # Ignore only the specific instance configuration pair with "-".
    df = df.replace("-", np.nan).infer_objects(copy=False)
    df = df.apply(pd.to_numeric)

    # df = df.clip(lower=0,upper=100)
    df = df.clip(upper=100)

    ###############################################################

    y = df.mean()
    y_std = df.std()
    x = range(len(y))

    n = df.count()
    y_se = y_std / np.sqrt(n)
    ci = 1.96 * y_se  # 95% confidence interval

    plt.figure(figsize=(12, 4.8))
    # plt.plot(x, y, label="Mean", color=colors_cmp[y_column], linestyle=styles_cmp[y_column], marker='o', linewidth=1)
    # plt.fill_between(x, y - ci, y + ci, color=colors_cmp[y_column], alpha=0.25)

    for xi, yi, ci_i in zip(x, y, ci):
        # draw CI as box
        plt.fill_between(
            [xi - 0.2, xi + 0.2],
            yi - ci_i,
            yi + ci_i,
            color=colors_cmp[y_column],
            alpha=0.3
        )

    # mean as point
    plt.plot(x, y, 'o', color=colors_cmp[y_column],markersize=6)

    # plt.xlabel(r'Follower behavior')
    # plt.ylabel(f'{y_column_name}')

    plt.xticks(x, label_behaviors, rotation=65)
    plt.margins(x=0.01)
    plt.yticks([0, 1, 2, 3, 4, 5, 6, 7, 8 , 9, 10, 11, 12])

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    for i, tick in enumerate(plt.gca().get_xticklabels()):
        if i < 5:
            tick.set_color('red')
        elif i < 21:
            tick.set_color('blue')
        else:
            tick.set_color('green')

    plt.savefig(cmp_values_path + name_common + "_align" + str(align) + "_bv" + str(bv_cfg) + f"_{y_column}" + "plot_mean.png", dpi=300)
    plt.close()

    ###############################################################

    plt.figure(figsize=(12, 4.8))

    plt.plot(x, df.median(), label="Median", color=colors_cmp[y_column], linestyle=styles_cmp[y_column], marker='o', linewidth=1)
    plt.fill_between(x, df.quantile(0.25), df.quantile(0.75), color=colors_cmp[y_column], alpha=0.25)

    # plt.xlabel(r'Follower behavior')
    # plt.ylabel(f'{y_column_name}')
    
    plt.xticks(x, label_behaviors, rotation=65)
    plt.margins(x=0.01)

    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # plt.savefig(cmp_values_path + name_common + "_align" + str(align) + "_bv" + str(bv_cfg)  + f"_{y_column}" + "line_median_.png", dpi=300)
    plt.close()

    ############################################################################

def plot_cmp_metric_two(bv_cfg_y, bv_cfg_x, column, column_name, align):
    """
    Plots the average of a metric (column) by follower behavior (comparison of two fixed behaviors), and saves the figure.
    
    Parameters:
    - column: Name of the column to use for the y-axis. ('gap_leader_','gap_follower_','gap_follower_near_')
    - column_name: Name to display on the y-axis 
    - align: [-1,0,1,2,3] -1 for all instances
    """

    # Set global font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        'text.usetex': True,
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'axes.linewidth': 1,
        'lines.markersize': 3
    })

    # Columns corresponding to metric column 
    columns = [column + get_name_behavior_config(column_bv_cfg) for column_bv_cfg in behaviors_config]

    # Rows for y and x axis.
    rows_y = get_name_behavior_config(bv_cfg_y)
    rows_x = get_name_behavior_config(bv_cfg_x)

    # Reading sheet.
    df = pd.read_excel(file_path, sheet_name="comp_values_all")

    if align >= 0 and align <= 3:
        df = df[(df["align"] >= align_ranges[align][0]) & (df["align"] <= align_ranges[align][1] - 1e-12)]

    df_y = df[df["config"] == rows_y]
    df_y = df_y[columns]

    df_x = df[df["config"] == rows_x]
    df_x = df_x[columns]

    # Ignore only the specific instance configuration pair with "-".
    df_y = df_y.replace("-", np.nan).infer_objects(copy=False)
    df_y = df_y.apply(pd.to_numeric)
    df_y = df_y.clip(upper=100)

    df_x = df_x.replace("-", np.nan).infer_objects(copy=False)
    df_x = df_x.apply(pd.to_numeric)
    df_x = df_x.clip(upper=100)

    ###############################################################

    y = df_y.mean()
    x = df_x.mean()

    plt.figure(figsize=(5, 5))

    # plt.plot(x, y, 'o', color=colors_cmp[column])
    for xi, yi, label in zip(x, y, label_behaviors):
        plt.plot(xi, yi, 'o', color=colors_cmp[column])
        plt.text(xi, yi, label, fontsize=8)

    lims = [
        min(min(x), min(y)),
        max(max(x), max(y))
    ]

    plt.plot(lims, lims, 'k--', linewidth=1, color="black")

    # plt.xlim(lims)
    # plt.ylim(lims)

    plt.xlabel(rf'{get_name_behavior_config_graphic(bv_cfg_x)}')
    plt.ylabel(rf'{get_name_behavior_config_graphic(bv_cfg_y)}')

    plt.legend()
    plt.grid(True)
    plt.tight_layout(pad=1.5)

    plt.savefig(cmp_values_path + name_common + "comp_align" + str(align) + "_bvx" + str(bv_cfg_x) + "_bvy" + str(bv_cfg_y) + f"_{column}" + "plot_mean.png", dpi=300)
    plt.close()

solver_columns = ['_saa','_dep']
gap_columns = ['gap_1','gap']
time_columns = ['time_lb','time']
solver_names = [r'$\texttt{SAA}$',r'$\texttt{DE}$']
solver_colors = ['#E69F00',  '#377EB8']
solver_styles = ['-', '--']

# Graphic between gap de and gap saa problem 1.
def plot_perf(tp):
    """
    Plots the accumulated value of a metric (y_column) by follower behavior, and saves the figure.
    
    opt. gap, time, and ub gap by percentage of instance for the different methods 
    (deterministic equivalent and transformed saa) for the different types {tp} of the decision-dependent 
    strong weak behavior, and saves the figures. 
    
    Parameters:
    - tp in [-1,0,1,2,3,4,5]: -1 i for all types
    """  

    # Set global font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        'text.usetex': True,
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 10,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 12,
        'axes.linewidth': 1,
        'lines.markersize': 3
    })  

    df = pd.read_excel(file_path, sheet_name="sol_all_1")
    if tp != -1:
        df = df[df["type"] == name_spec_types[1][tp]]
    x = [i/df.shape[0]*100 for i in range(1,df.shape[0]+1)]
    
    ################################
    
    # Optimality gap.
    plt.figure(figsize=(3.4, 2.4))

    for sol in [0,1]:
        y_column = gap_columns[sol] + solver_columns[sol]
        df[y_column] = df[y_column].clip(upper=100) # max gap 100% (avoid problems as we have positive and negative values)
        y = df[y_column].sort_values()

        plt.plot(x, y, label=solver_names[sol], color=solver_colors[sol], linestyle=solver_styles[sol], linewidth=1.5)
    
    plt.xlabel(r'Instance-Behavior-Parameter Configurations (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{opt}$ (\%)')
    
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_opt_gap_tp" + str(tp) + ".png", dpi=300)
    plt.close() 
    
    ##############################
        
    # Time
    plt.figure(figsize=(3.4, 2.4))
    for sol in [0,1]:
        y_column = time_columns[sol] + solver_columns[sol]
        y = df[y_column].sort_values()

        plt.plot(x, y, label=solver_names[sol], color=solver_colors[sol], linestyle=solver_styles[sol], linewidth=1.5)
        
    plt.xlabel(r'Instance-Behavior-Parameter Configurations (\%)')
    # plt.ylabel(r'Time (seconds)')
    
    plt.legend()
    plt.yticks([0, 1800, 3600, 5400, 7200, 9000, 10800])

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_time_tp" + str(tp) + ".png", dpi=300)
    plt.close()
    
    ############################ 
    
    # Oportunity gap
    # plt.figure(figsize=(4, 2.5))
    plt.figure(figsize=(3.2, 2.4))
    y_column = 'gap_ub'
    df[y_column] = df[y_column].clip(lower=-100,upper=100) # max gap 100% (avoid problems as we have positive and negative values)
    y = df[y_column].sort_values()

    plt.plot(x, y, label = solver_names[sol], color=solver_colors[sol], linestyle=solver_styles[sol], linewidth=1.5)
            
    plt.xlabel(r'Instance-Behavior-Parameter Configurations (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{ub^*}(\%)$')

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    # plt.savefig(perf_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + "_ub_gap_" + str(tp) + ".png", dpi=300)
    plt.close() 

behavior_names = [r'$\beta$',r'$\beta(\cdot)$',r'$\pi(\cdot)$']
behavior_styles = ['-', '--', '-.']
behavior_colors = ["#882255", "#44AA99", "#DDCC77"]

# Graphic time and gaps per instance for different behaviors: 3 lines for behavior general types.
def plot_perf_all_bv():
    """
    Plots the accumulated (mean per instance and behavior general type) opt. gap and time by percentage of instance for the three
    different general behavior types (fixed strong weak, decision-dependent strong weak,
    generalized decision-dependent), and saves the figures. 

    """ 

    # Set font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        'text.usetex': True,
        'font.size': 7.5,
        'axes.titlesize': 8,
        'axes.labelsize': 7.5,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'axes.linewidth': 0.5,
        'lines.linewidth': 1.0,
        'lines.markersize': 2.0,
    })

    df = list()
    df.append(pd.read_excel(file_path, sheet_name="sol_all_0"))
    df.append(pd.read_excel(file_path, sheet_name="sol_all_1"))
    df.append(pd.read_excel(file_path, sheet_name="sol_all_2"))

    x = [i/df[0][["instance", "seed"]].drop_duplicates().shape[0]*100 for i in range(1,df[0][["instance", "seed"]].drop_duplicates().shape[0]+1)]

    ################################
    
    # Optimality gap.
    plt.figure(figsize=(3.4, 2.4))

    for bv in behaviors:
        y = None
        if bv == 0:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["gap_dep"].mean())["gap_dep"].sort_values()
            y = y.clip(upper=100)
        if bv == 1 or bv == 2:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["gap_1_saa"].mean())["gap_1_saa"].sort_values()
            y = y.clip(upper=100)

        plt.plot(x, y, label=behavior_names[bv], color=behavior_colors[bv], linestyle=behavior_styles[bv], linewidth=1.5)
    
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{opt}$ (\%)')

    plt.yticks([0, 10, 20, 30])
    
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_opt_gap_bv_all.png", dpi=300)
    plt.close() 

    ##############################
        
    # Time
    plt.figure(figsize=(3.4, 2.4))
    for bv in behaviors:
        y = None
        if bv == 0:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["time_dep"].mean())["time_dep"].sort_values()
        if bv == 1:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["time_lb_saa"].mean())["time_lb_saa"].sort_values()
        if bv == 2:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["time_lb_saa"].mean())["time_lb_saa"].sort_values()
    
        plt.plot(x, y, label=behavior_names[bv], color=behavior_colors[bv], linestyle=behavior_styles[bv], linewidth=1.5)
        
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'Time (seconds)')

    plt.legend()
    plt.yticks([0, 1800, 3600, 5400, 7200, 9000, 10800])

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_time_bv_all.png", dpi=300)
    plt.close()
    
    ############################ 

spec_types_plot = {0: fix_strwk_coop, 1: spec_types[1], 2: spec_types[2]}
behavior_spec_names = {bv: {tp: get_name_behavior_config_general_graphic(bv,tp) for tp in spec_types_plot[bv]} for bv in behaviors}
behavior_spec_colors = {0: {0.0:"#CC79A7", 0.3:"#F0E442", 0.5:"#999999", 0.7:"#56B4E9", 1.0:"#D55E00"},
                        1: {0:"#882255", 1:"#44AA99", 2:"#DDCC77", 3:"#AA4499", 4:"#117733", 5:"#332288"},
                        2: {0:"#661100", 1:"#6699CC", 2:"#AA4466", 3:"#4477AA", 4:"#228833"}}
behavior_spec_styles = {0: {0.0:"-", 0.3:"--", 0.5:"-.", 0.7:":", 1.0:(0, (5, 2))}, 
                        1: {0:"-", 1:"--", 2:"-.", 3:":", 4:(0, (5, 2)), 5:(0, (3, 1, 1, 1))},
                        2: {0:"-", 1:"--", 2:"-.", 3:":", 4:(0, (5, 2))}}

# Graphic time and gaps per instance for different behaviors: Multiple lines for subypes inside general type of behavior.
def plot_perf_spec_bv(bv,sol=""):
    """
    Plots the accumulated (mean per instance and behavior specific subtypes) opt. gap and time by percentage of instance for 
    the different subtypes of each general behavior type, and saves the figures. 

    """ 

    # Set font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        'text.usetex': True,
        'font.size': 7.5,
        'axes.titlesize': 8,
        'axes.labelsize': 7.5,
        'xtick.labelsize': 6.5,
        'ytick.labelsize': 6.5,
        'legend.fontsize': 6.5,
        'axes.linewidth': 0.5,
        'lines.linewidth': 1.0,
        'lines.markersize': 2.0,
    })

    df = pd.read_excel(file_path, sheet_name="sol_all_"+ str(bv))

    x = [i/df[["instance", "seed"]].drop_duplicates().shape[0]*100 for i in range(1,df[["instance", "seed"]].drop_duplicates().shape[0]+1)]

    ################################

    # Optimality gap.
    plt.figure(figsize=(2.2, 1.7))

    for tp in spec_types_plot[bv]:
        y = None
        if bv == 0:
            df_tp = df[df["params"] == tp]
            y = (df_tp.groupby(["instance","seed"], as_index=False)["gap_dep"].mean())["gap_dep"].sort_values()
        if bv == 2:
            df_tp = df[df["type"] == name_spec_types[bv][tp]]
            y = (df_tp.groupby(["instance", "seed"], as_index=False)["gap_1_saa"].mean())["gap_1_saa"].sort_values()
        if bv == 1:
            df_tp = df[df["type"] == name_spec_types[bv][tp]]
            if sol == "_dep":
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["gap_dep"].mean())["gap_dep"].sort_values().clip(upper=100) # max gap 100%
            if sol == "_saa":
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["gap_1_saa"].mean())["gap_1_saa"].sort_values()

        plt.plot(x, y, label=behavior_spec_names[bv][tp], color=behavior_spec_colors[bv][tp], linestyle=behavior_spec_styles[bv][tp], linewidth=1)
    
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{opt}$ (\%)')

    if (bv == 0) or (bv == 2) or (bv == 1 and sol == "_saa"):
        if bv == 0:
            plt.yticks([0, 20, 40, 60, 80, 100])
            plt.ylim(-2, 102)
        else: 
            plt.yticks([0, 20, 40, 60, 80, 100], [])
            plt.ylim(-2, 102)
    if (bv == 1 and sol == "_dep"):
        plt.yticks([0, 20, 40, 60, 80, 100])
    
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_opt_gap_bv" + str(bv) + sol + ".png", dpi=300)
    plt.close() 

    ##############################
    
    # Time
    plt.figure(figsize=(2.2, 1.7))
    for tp in spec_types_plot[bv]:
        y = None
        if bv == 0:
            df_tp = df[df["params"] == tp]
            y = (df_tp.groupby(["instance", "seed"], as_index=False)["time_dep"].mean())["time_dep"].sort_values() / 100
        if bv == 1:
            df_tp = df[df["type"] == name_spec_types[bv][tp]]
            if sol == "_dep":
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["time_dep"].mean())["time_dep"].sort_values() / 100
            if sol == "_saa":
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["time_lb_saa"].mean())["time_lb_saa"].sort_values() / 100
        if bv == 2:
            df_tp = df[df["type"] == name_spec_types[bv][tp]]
            y = (df_tp.groupby(["instance", "seed"], as_index=False)["time_lb_saa"].mean())["time_lb_saa"].sort_values() / 100
    
        plt.plot(x, y, label=behavior_spec_names[bv][tp], color=behavior_spec_colors[bv][tp], linestyle=behavior_spec_styles[bv][tp], linewidth=1)
        
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'Time (100 seconds)')
    # plt.yscale('log') 

    plt.legend()
    if (bv == 0) or (bv == 2) or (bv == 1 and sol == "_saa"):
        if bv == 0:
            plt.yticks([0, 18, 36, 54, 72, 90, 108])
        else:
            plt.yticks([0, 18, 36, 54, 72, 90, 108], [])
    if (bv == 1 and sol == "_dep"):
        plt.yticks([0, 18, 36, 54, 72, 90, 108])

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_time_bv" + str(bv) + sol + ".png", dpi=300)
    plt.close()

############################################################################################################
#                                              Main
############################################################################################################

def main():
    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    os.makedirs(os.path.dirname(sol_values_path), exist_ok=True)
    os.makedirs(os.path.dirname(cmp_values_path), exist_ok=True)
    os.makedirs(os.path.dirname(perf_path), exist_ok=True)
    os.makedirs(os.path.dirname(perf_path_bv), exist_ok=True)

    # Plotting mean and average of optimal value comparison gap metrics.
    plot_sol_metric('leader_obj_',r'$\overline{\texttt{gap}}^{\mu}$ (\%)',-1)
    plot_sol_metric('follower_obj_',r'$\overline{\texttt{gap}}_{\varphi}^{\mu}$ (\%)',-1)
    plot_sol_metric('follower_obj_near_mean_',r'$\overline{\texttt{gap}}_{f}^{\mu}$ (\%)',-1)
    if plot_align == True:
        for align in range(0,align_types):
            plot_sol_metric('follower_obj_',r'$\overline{\texttt{gap}}_{\varphi}^{\mu}$ (\%)',align)
            plot_sol_metric('follower_obj_near_mean_',r'$\overline{\texttt{gap}}_{f}^{\mu}$ (\%)',align)

    # Plotting optimal solution comparison based on support-endogenous (fixed strong-weak) behaviors and neutral approach only and align = -1 and leader objective.
    for pr in behavior_params[0][0]:
        plot_cmp_metric((0,0,pr),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',-1)
    plot_cmp_metric((1,0,0.0),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',-1)
    plot_cmp_metric((2,0,[0.0, 0.0]),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',-1)

    # for bv in behaviors_config:
    #     plot_cmp_metric(bv,'gap_leader_',r'$\overline{\texttt{gap}}_{F}^{b \leftarrow b}$ (\%)',-1)

    # if plot_align == True:
    #     for align in range(0,align_types):
    #         for pr in behavior_params[0][0]:
    #             plot_cmp_metric((0,0,pr),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',align)
    #         plot_cmp_metric((1,0,0.0),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',align)
    #         plot_cmp_metric((2,0,[0.0, 0.0]),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',align)
    #         for bv in behaviors_config:
    #             plot_cmp_metric(bv,'gap_leader_',r'$\overline{\texttt{gap}}_{F}^{b \leftarrow b}$ (\%)',align)
    
    # plot_cmp_metric_two((0,0,0.5),(1,0,0.0),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',-1)
    # plot_cmp_metric_two((0,0,0.5),(2,0,[0.0, 0.0]),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',-1)
    # plot_cmp_metric_two((0,0,0.5),(1,2,10.0),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',-1)
    # plot_cmp_metric_two((0,0,0.5),(1,3,10.0),'gap_leader_',r'$\overline{\texttt{gap}}^{\mu \leftarrow \tilde{\mu}}$ (\%)',-1)

    # Plotting performance comparison DE and SAA for decision-dependent strong weak case. (all subtypes together)
    plot_perf(-1) 
    
    # Plotting performance comparison each type and subtype of behavior.
    plot_perf_all_bv()
    for bv in behaviors:
        if bv == 1:
            plot_perf_spec_bv(bv,"_saa")
            plot_perf_spec_bv(bv,"_dep")
        else:
            plot_perf_spec_bv(bv)

    plot_extra = False
    if plot_extra == True:
        # Plotting accumulated version of optimal value comparison gap metrics.
        for bv in behaviors:
            plot_sol_metric_acc_per_bv('follower_obj_',r'$\overline{\texttt{gap}}_{\varphi}^{\mu}$ (\%)',bv,-1)
            plot_sol_metric_acc_per_bv('follower_obj_near_mean_',r'$\overline{\texttt{gap}}_{f}^{\mu}$ (\%)',bv,-1)
            for align in range(0,align_types):
                plot_sol_metric_acc_per_bv('follower_obj_',r'$\overline{\texttt{gap}}_{f}^{b \mid opt}$ (\%)',bv,align)
                plot_sol_metric_acc_per_bv('follower_obj_near_mean_',r'$\overline{\texttt{gap}}_{f_\epsilon}^{b \mid opt}$ (\%)',bv,align)
    
        # Plotting comparison for all behaviors, align values and metrics.  
        for bv in behaviors_config:
            for align in range(-1,align_types):
                plot_cmp_metric(bv,'gap_leader_',r'$\overline{\texttt{gap}}_{F}^{b \leftarrow b}$ (\%)',align)
                plot_cmp_metric(bv,'gap_follower_',r'$\overline{\texttt{gap}}_{f}^{b \leftarrow b}$ (\%)',align)
                plot_cmp_metric(bv,'gap_follower_near_',r'$\overline{\texttt{gap}}_{f_\epsilon}^{b \leftarrow b}$ (\%)',align)
        
        # Plotting performance comparison DE and SAA for decision-dependent strong weak case. (separated subtypes)
        for tp in spec_types[1]:
            plot_perf(tp)

if __name__ == "__main__":
    main()