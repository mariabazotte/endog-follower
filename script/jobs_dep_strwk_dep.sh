#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=20:00:00
#SBATCH --array=1-3230
#SBATCH --output=../slurm/strwk/dep/%A_%a.out

i=1
verbose=1
nbthreads=1
timelimit=10800

nbvalidateproblems=4
nbvalidatescenarios=10000
nbvalidatethinning=50

# "denegre-miblp_20_15_50_0110_10_" "denegre-miblp_20_20_50_0110_5_" 
# "denegre-miblp_20_20_50_0110_10_" "denegre-miblp_20_20_50_0110_15_" 
# "denegre-milp_30_30_30_2310_5_" "denegre-milp_30_30_30_2310_10_"
# "denegre-miblp_20_20_50_0110_5_e1_" "denegre-milp_30_30_30_2310_5_e1_" 
# "denegre-milp_30_30_30_2310_5_e2_" "denegre-milp_30_30_30_2310_10_e1_" 
# "denegre-milp_30_30_30_2310_10_e2_"
# "xuwang-bmilplib_10_" "xuwang-bmilplib_60_" 
# "xuwang-bmilplib_110_" "xuwang-bmilplib_160_" 
# "xuwang-bmilplib_210_" "xuwang-bmilplib_260_"

for instance in "denegre-miblp_20_15_50_0110_10_" "denegre-miblp_20_20_50_0110_5_" "denegre-miblp_20_20_50_0110_10_" "denegre-miblp_20_20_50_0110_15_" "denegre-milp_30_30_30_2310_5_" "denegre-milp_30_30_30_2310_10_" "denegre-miblp_20_20_50_0110_5_e1_" "denegre-milp_30_30_30_2310_5_e1_" "denegre-milp_30_30_30_2310_5_e2_" "denegre-milp_30_30_30_2310_10_e1_" "denegre-milp_30_30_30_2310_10_e2_" "xuwang-bmilplib_10_" "xuwang-bmilplib_60_" "xuwang-bmilplib_110_" "xuwang-bmilplib_160_" "xuwang-bmilplib_210_" "xuwang-bmilplib_260_"
do
for epsnearopt in 0.01
do
for seed in 1 2 3 4 5 6 7 8 9 10
do  
    # Decision-dependent Strong-Weak.
    for approach in 1
    do
        # Proportional
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 1 -dep_strwk 0 -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
        fi
        (( i = $i +1 ))

        # Threshold
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 1 -dep_strwk 1 -thr_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Strong/Flexible
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 1 -dep_strwk 2 -str_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Fragile/Inflexible
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 1 -dep_strwk 3 -frag_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Strong power
        for param in 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 1 -dep_strwk 4 -str_pow_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Fragile power
        for param in 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 1 -dep_strwk 5 -frag_pow_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done
    done
done
done
done