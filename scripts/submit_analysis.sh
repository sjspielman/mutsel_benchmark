#! /bin/bash

# Submit a full analysis for a given empirical dataset. It is assumed that the alignment is in two formats, named <dataset name>.fasta and <dataset name>.phy, and the tree is named <dataset name>.tre.

if [ $# != 2 ]; then
    echo "Usage: sh submit_analysis.sh <dataset name> <type>, where <dataset name> is the data name and <type> is either sim (for simulated) or emp (for empirical)."
    exit 1
fi
DATA=$1
TYPE=$2


if [[ $TYPE == "sim" ]]; then 
    TOPDIR=$HOME/mutsel_bench/data/simulation
    ALN=${DATA}_simulated
    TREE="ntaxa1024_bl0.5.tre"
else
    TOPDIR=$HOME/mutsel_bench/data/empirical
    ALN=$DATA
    TREE=$DATA.tre
    
fi 

# FEL 1rate
#sed -i "s/-N JOB/-N fel1_${DATA}/" hyphy/fel1.qsub
#qsub hyphy/fel1.qsub $TOPDIR $ALN $TREE
#sed -i "s/-N fel1_${DATA}/-N JOB/" hyphy/fel1.qsub

# FEL 2rate, sjs
#sed -i "s/-N JOB/-N fel2_${DATA}/" hyphy/fel2.qsub
#qsub hyphy/fel2.qsub $TOPDIR $ALN $TREE
#sed -i "s/-N fel2_${DATA}/-N JOB/" hyphy/fel2.qsub


# FEL 2rate, default
#sed -i "s/-N JOB/-N fel2d_${DATA}/" hyphy/fel2_default.qsub
#qsub hyphy/fel2_default.qsub $TOPDIR $ALN $TREE
#sed -i "s/-N fel2d_${DATA}/-N JOB/" hyphy/fel2_default.qsub


# swmutsel
#sed -i "s/-N JOB/-N sw_${DATA}/" swmutsel/swmutsel.qsub
#qsub swmutsel/swmutsel.qsub $TOPDIR $ALN $TREE
#sed -i "s/-N sw_${DATA}/-N JOB/" swmutsel/swmutsel.qsub

 
# phylobayes
sed -i "s/-N JOB/-N pb_${DATA}/" phylobayes/phylobayes.qsub
qsub phylobayes/phylobayes.qsub $TOPDIR $ALN $TREE
sed -i "s/-N pb_${DATA}/-N JOB/" phylobayes/phylobayes.qsub



