#! /bin/bash

# Submit a full analysis for a given datafile. It is assumed that the alignment is named <dataset name>.phy and the tree is named <dataset name>.tre

if [ $# != 1 ]; then
    echo "Usage: sh submit_analysis.sh <dataset name>"
    exit 1
fi
DATA=$1

# # SLAC
sed -i "s/-N JOB/-N ${DATA}_omega/" slac_qsub.sh
qsub slac_qsub.sh $DATA
sed -i "s/-N ${DATA}_omega/-N JOB/" slac_qsub.sh

# swmutsel, with asymmetric mutation rates
sed -i "s/-N JOB/-N ${DATA}_estmu_sw/" swmutsel_qsub.sh
qsub swmutsel_qsub.sh $DATA 1
sed -i "s/-N ${DATA}_estmu_sw/-N JOB/" swmutsel_qsub.sh

# swmutsel, with symmetric mutation rates
sed -i "s/-N JOB/-N ${DATA}_treeonly_sw/" swmutsel_qsub.sh
qsub swmutsel_qsub.sh $DATA 0
sed -i "s/-N ${DATA}_treeonly_sw/-N JOB/" swmutsel_qsub.sh
# 
# # phylobayes
sed -i "s/-N JOB/-N ${DATA}_pb/" phylobayes_qsub.sh
qsub phylobayes_qsub.sh $DATA
sed -i "s/-N ${DATA}_pb/-N JOB/" phylobayes_qsub.sh



