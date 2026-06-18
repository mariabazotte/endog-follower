#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=20:00:00
#SBATCH --array=1-850
#SBATCH --output=../slurm/strwk/fix/%A_%a.out

i=1
verbose=1
nbthreads=1
timelimit=10800

nbvalidateproblems=4
nbvalidatescenarios=10000
nbvalidatethinning=50 

# "denegre-miblp_20_15_50_0110_10_" "denegre-miblp_20_20_50_0110_5_" "denegre-miblp_20_20_50_0110_10_" "denegre-miblp_20_20_50_0110_15_" "denegre-milp_30_30_30_2310_5_" "denegre-milp_30_30_30_2310_10_" 
# "denegre-miblp_20_20_50_0110_5_e1_" "denegre-milp_30_30_30_2310_5_e1_" "denegre-milp_30_30_30_2310_5_e2_" "denegre-milp_30_30_30_2310_10_e1_" "denegre-milp_30_30_30_2310_10_e2_"
# "xuwang-bmilplib_10_" "xuwang-bmilplib_60_" "xuwang-bmilplib_110_" "xuwang-bmilplib_160_" "xuwang-bmilplib_210_" "xuwang-bmilplib_260_"

for instance in "denegre-miblp_20_15_50_0110_10_" "denegre-miblp_20_20_50_0110_5_" "denegre-miblp_20_20_50_0110_10_" "denegre-miblp_20_20_50_0110_15_" "denegre-milp_30_30_30_2310_5_" "denegre-milp_30_30_30_2310_10_" "denegre-miblp_20_20_50_0110_5_e1_" "denegre-milp_30_30_30_2310_5_e1_" "denegre-milp_30_30_30_2310_5_e2_" "denegre-milp_30_30_30_2310_10_e1_" "denegre-milp_30_30_30_2310_10_e2_" "xuwang-bmilplib_10_" "xuwang-bmilplib_60_" "xuwang-bmilplib_110_" "xuwang-bmilplib_160_" "xuwang-bmilplib_210_" "xuwang-bmilplib_260_"
do
for epsnearopt in 0.01
do
for seed in 1 2 3 4 5 6 7 8 9 10
do  
    # Fixed Strong-Weak.
    for fixcoop in 0.0 0.3 0.5 0.7 1.0
    do
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 0 -fix_strwk $fixcoop -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
        fi
        (( i = $i +1 ))
    done
done
done
done