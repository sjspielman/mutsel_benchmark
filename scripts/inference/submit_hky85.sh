#! /bin/bash
# Fit HKY85 model to empirical datasets, for "empirical" dN/dS prediction.

if [ $# != 2 ]; then
    echo "Usage: sh submit_hky85.sh <dataset name>
fi
DATA=$1
TYPE=$2


TOPDIR=$HOME/mutsel_benchmark/data/empirical
ALN=$DATA
TREE=$DATA.tre

sed -i "s/-N JOB/-N hky85_${DATA}/" hyphy_hky85/hky85.qsub
qsub hyphy_hky85/hky85.qsub $TOPDIR $ALN $TREE
sed -i "s/-N hky85_${DATA}/-N JOB/" hyphy_hky85/hky85.qsub


