#! /bin/bash

# Submit a full analysis for a given datafile. It is assumed that the alignment is in two formats, named <dataset name>.fasta and <dataset name>.phy, and the tree is named <dataset name>.tre.

if [ $# != 1 ]; then
    echo "Usage: sh submit_analysis.sh <dataset name>"
    exit 1
fi
DATA=$1

# FEL 1rate
#sed -i "s/-N JOB/-N ${DATA}_fel1/" fel1.qsub
#qsub fel1.qsub $DATA
#sed -i "s/-N ${DATA}_fel1/-N JOB/" fel1.qsub

# FEL 2rate, sjs
#sed -i "s/-N JOB/-N ${DATA}_fel2/" fel2.qsub
#qsub fel2.qsub $DATA
#sed -i "s/-N ${DATA}_fel2/-N JOB/" fel2.qsub


# FEL 2rate, default
sed -i "s/-N JOB/-N ${DATA}_fel2_default/" fel2_default.qsub
qsub fel2_default.qsub $DATA
sed -i "s/-N ${DATA}_fel2_default/-N JOB/" fel2_default.qsub


# swmutsel
#sed -i "s/-N JOB/-N ${DATA}_sw/" swmutsel.qsub
#qsub swmutsel.qsub $DATA
#sed -i "s/-N ${DATA}_sw/-N JOB/" swmutsel.qsub

 
# phylobayes
#sed -i "s/-N JOB/-N ${DATA}_pb/" phylobayes.qsub
#qsub phylobayes.qsub $DATA
#sed -i "s/-N ${DATA}_pb/-N JOB/" phylobayes.qsub



