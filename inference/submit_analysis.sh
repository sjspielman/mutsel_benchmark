#! /bin/bash
# Submit MutSel model inferences across datasets


function submit_swmutsel {
    DATA=$1
    TREE=$2
    
    for M in mvn1 mvn10 mvn100 d1.0 d0.1 d0.01 nopenal; do
        sed -i "s/-N JOB/-N sw_${DATA}_${M}/" swmutsel.qsub
        qsub swmutsel.qsub $TOPDIR $DATA $TREE $M
        sed -i "s/-N sw_${DATA}_${M}/-N JOB/" swmutsel.qsub
    done
}

function submit_pbmutsel {
    DATA=$1
    TREE=$2

    sed -i "s/-N JOB/-N pb_${DATA}/" phylobayes.qsub
    qsub phylobayes.qsub $DATA $TREE
    sed -i "s/-N pb_${DATA}/-N JOB/" phylobayes.qsub
}

 
for NAME in 1B4T_A 1G58_B 1GV3_A 1IBS_A 1R6M_A 1RII_A 1V9S_B 1W7W_B 2BCG_Y 2CFE_A 2FLI_A; do 
    for BL in 0.01 0.5; do
        for DEL in strong weak; do
        
            DATA=${NAME}_del${DEL}_bl${BL}
            TREE=n512_bl${BL}.tre
            submit_pbmutsel $DATA $TREE
            submit_swmutsel $DATA $TREE
        done
    done
done


for NAME in HA NP; do 
    for BL in 0.01 0.5; do
        
         DATA=${NAME}_bl${BL}
         TREE=n512_bl${BL}.tre
         submit_pbmutsel $DATA $TREE
         submit_swmutsel $DATA $TREE
    done
done

