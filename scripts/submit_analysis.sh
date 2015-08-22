#! /bin/bash

# Submit a full analysis for a given datafile. It is assumed that the alignment is named <dataset name>.phy and the tree is named <dataset name>.tre

if [ $# != 1 ]; then
    echo "Usage: sh submit_analysis.sh <dataset name>"
    exit 1
fi
DATA=$1

# FEL1,2
sed -i "s/-N JOB/-N ${DATA}_fel1/" fel1.qsub
qsub fel1.qsub $DATA
sed -i "s/-N ${DATA}_fel1/-N JOB/" fel1.qsub

sed -i "s/-N JOB/-N ${DATA}_fel2/" fel2.qsub
qsub fel2.qsub $DATA
sed -i "s/-N ${DATA}_fel2/-N JOB/" fel2.qsub

# swmutsel, with estimated mutation rates
sed -i "s/-N JOB/-N ${DATA}_estmu_sw/" swmutsel.qsub
qsub swmutsel.qsub $DATA 1
sed -i "s/-N ${DATA}_estmu_sw/-N JOB/" swmutsel.qsub

# phylobayes
sed -i "s/-N JOB/-N ${DATA}_pb/" phylobayes.qsub
qsub phylobayes.qsub $DATA
sed -i "s/-N ${DATA}_pb/-N JOB/" phylobayes.qsub






# swmutsel, with symmetric mutation rates
#sed -i "s/-N JOB/-N ${DATA}_treeonly_sw/" swmutsel_qsub.sh
#qsub swmutsel_qsub.sh $DATA 0
#sed -i "s/-N ${DATA}_treeonly_sw/-N JOB/" swmutsel_qsub.sh
# 

