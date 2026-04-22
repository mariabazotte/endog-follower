#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=5:00:00
#SBATCH --array=1-920
#SBATCH --output=../slurm/test_%A_%a.out
#SBATCH --account=def-vidalthi

i=1
verbose=1
nbthreads=1
timelimit=7200

nbproblems=5
nbscenarios=500 
nbthinning=2

nbvalidateproblems=5
nbvalidatescenarios=10000
nbvalidatethinning=50

instance="denegre-miblp_20_15_50_0110_10_"

for epsnearopt in 0.01
do
for seed in 1 2 3 4 5 6 7 8 9 10
do  
    # Fixed Strong-Weak.
    for fixcoop in 0.0 0.3 0.5 0.7 1.0
    do
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 0 -fix_strwk $fixcoop -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
        fi
        (( i = $i +1 ))
    done

    # Decision-dependent Strong-Weak.
    for approach in 0 1
    do
        # Proportional
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 1 -dep_strwk 0 -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
        fi
        (( i = $i +1 ))

        # Threshold
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 1 -dep_strwk 1 -thr_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Strong/Flexible
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 1 -dep_strwk 2 -str_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Fragile/Inflexible
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 1 -dep_strwk 3 -frag_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Strong power
        for param in 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 1 -dep_strwk 4 -str_pow_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        # Fragile power
        for param in 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 1 -dep_strwk 5 -frag_pow_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done
    done

    # Decision-dependent general.

    # Neutral 
    if [ $SLURM_ARRAY_TASK_ID -eq $i ]
    then
        ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 2 -dep_gen 0 -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
    fi
    (( i = $i +1 ))

    for nbint in 4 8 12
    do 
        # Proportional
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 2 -dep_gen 1 -nb_intv $nbint -scal_param 1.0 -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning  -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
        fi
        (( i = $i +1 ))

        for scaling in 0.5 2.0 5.0 10.0 20.0 40.0 80.0
        do  
            # Fragile
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 2 -dep_gen 2 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        for scaling in 1.0 2.0 5.0 10.0
        do  
            # Fragile power
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 2 -dep_gen 3 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done

        for scaling in 1.0 2.0 5.0 10.0 
        do  
            # Strong power
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}.txt" -follower 2 -dep_gen 4 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning
            fi
            (( i = $i +1 ))
        done
    done
done
done