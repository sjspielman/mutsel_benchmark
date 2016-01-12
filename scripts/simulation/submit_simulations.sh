# Submit simulations to phylocluster:

for DATA in 1B4T_A 1G58_B 1GV3_A 1IBS_A 1R6M_A 1RII_A 1V9S_B 1W7W_B 2BCG_Y 2CFE_A 2FLI_A; do 
    for DEL in strong weak weakweak; do
        NAME=${DATA}_del${DEL}
        qsub simulate_alignments.qsub $NAME
    done
done

