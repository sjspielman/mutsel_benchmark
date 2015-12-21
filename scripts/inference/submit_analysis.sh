#! /bin/bash
# Submit MutSel model inferences across datasets


function submit_swmutsel {
    TOPDIR=$1
    DATA=$2
    TREE=$3

    for M in nopenal mvn1 mvn10 mvn100 d0.01 d0.1 d1.0; do
        sed -i "s/-N JOB/-N sw_${DATA}_${M}/" swmutsel/swmutsel.qsub
        qsub swmutsel/swmutsel.qsub $TOPDIR $DATA $TREE $M
        sed -i "s/-N sw_${DATA}_${M}/-N JOB/" swmutsel/swmutsel.qsub
    done
}

function submit_pbmutsel {
    TOPDIR=$1
    DATA=$2
    TREE=$3

    sed -i "s/-N JOB/-N pb_${DATA}/" phylobayes/phylobayes.qsub
    qsub phylobayes/phylobayes.qsub $TOPDIR $DATA $TREE
    sed -i "s/-N pb_${DATA}/-N JOB/" phylobayes/phylobayes.qsub
}

 
# Inferences for simulated data
TOPDIR=$HOME/mutsel_benchmark/data/simulation
TREE="n512_bl0.5.tre"
for NAME in 1B4T_A 1G58_B 1GV3_A 1IBS_A 1R6M_A 1RII_A 1V9S_B 1W7W_B 2BCG_Y 2CFE_A 2FLI_A; do 
    for DEL in strong weak; do
        
        DATA=${NAME}_del${DEL}
        
        submit_swmutsel $TOPDIR $DATA $TREE
        submit_pbmutsel $TOPDIR $DATA $TREE

    done
done



# Inferences for empirical data
TOPDIR=$HOME/mutsel_benchmark/data/empirical
for NAME in PF00106 PF00149 PF00188 PF00300 PF00512 PF00753 PF01261 PF01551 PF01636 PF03144 PF00141 PF00158 PF00226 PF00482 PF00520 PF01061 PF01546 PF01584 PF02775 PF04542 PF00168 PF00486 PF00535 PF01590 PF03466 PF00271 PF00501 PF00571 PF00593 PF00126 PF01266 PF01336 PF01926 PF02518 PF04055 PF07715; do 
        
    DATA=$NAME
    TREE=${NAME}.tre
        
    submit_swmutsel $TOPDIR $DATA $TREE
    submit_pbmutsel $TOPDIR $DATA $TREE

    done
done
