#! /bin/bash
# Submit MutSel model inferences across datasets


function submit_swmutsel {
    DATA=$1
    TREE=$2

    for M in nopenal mvn1 mvn10 mvn100 d0.01 d0.1 d1.0; do
        sed -i "s/-N JOB/-N sw_${DATA}_${M}/" swmutsel/swmutsel.qsub
        qsub swmutsel/swmutsel.qsub $TOPDIR $DATA $TREE $M
        sed -i "s/-N sw_${DATA}_${M}/-N JOB/" swmutsel/swmutsel.qsub
    done
}

function submit_pbmutsel {
    DATA=$1
    TREE=$2

    sed -i "s/-N JOB/-N pb_${DATA}/" phylobayes/phylobayes.qsub
    qsub phylobayes/phylobayes.qsub $DATA $TREE
    sed -i "s/-N pb_${DATA}/-N JOB/" phylobayes/phylobayes.qsub
}

 
TREE="n512_bl0.5.tre"
for NAME in 1B4T_A 1G58_B 1GV3_A 1IBS_A 1R6M_A 1RII_A 1V9S_B 1W7W_B 2BCG_Y 2CFE_A 2FLI_A; do 
    for DEL in strong weak weakweak; do
        
        DATA=${NAME}_del${DEL}
        
        submit_swmutsel $DATA $TREE
        submit_pbmutsel $DATA $TREE

    done
done
