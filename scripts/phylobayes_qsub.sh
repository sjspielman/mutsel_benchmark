#$ -N JOB
#$ -e e_$JOB_NAME
#$ -o o_$JOB_NAME
#$ -S /bin/bash
#$ -q wilke
#$ -m beas
#$ -pe serial 16
source ~/.bashrc

# define files, variables used in analysis
CPU=16
DATA=$1
ALN=$DATA.phy
TREE=$DATA.tre

PBEXEC=$HOME/bin/pb_mpi1.5a/data/pb_mpi
READPBEXEC=$HOME/bin/pb_mpi1.5a/data/readpb_mpi

DATADIR=$HOME/mutsel_bench/data
SCRIPTDIR=$HOME/mutsel_bench/scripts

# create and enter working directory
WDIR=/state/partition1/sjs3495/$JOB_NAME-$JOB_ID
mkdir -p $WDIR
cd $WDIR

# copy files, executables here
cp $DATADIR/$ALN .
cp $DATADIR/$TREE .
cp $SCRIPTDIR/run_phylobayes.py .
cp $SCRIPTDIR/dnds_mutsel_functions.py .
cp $PBEXEC .
cp $READPBEXEC .

# run analysis
module load python
export PYTHONPATH=$HOME/bin/lib/python2.7/site-packages # to get biopython behaving w/ phylip-relaxed and to get dendropy
python run_phylobayes.py $ALN $TREE $CPU ${DATA}_phylobayes

FINALDIR=$HOME/mutsel_bench/results/
mkdir -p $FINALDIR
cp *dnds.txt $FINALDIR
cp *.aap $FINALDIR/
cp *.trace $FINALDIR

cd ..
cp -r $WDIR /share/WilkeLab/work/sjs3495/
rm -r $WDIR



