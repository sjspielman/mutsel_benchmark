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
ESTMU=$2
ALN=$DATA.phy
TREE=$DATA.tre

SWEXEC=$HOME/bin/swmutsel-alpha_jdk1.6.jar

DATADIR=$HOME/MutSel_benchmark/data
SCRIPTDIR=$HOME/MutSel_benchmark/scripts

# create and enter working directory
WDIR=/state/partition1/sjs3495/$JOB_NAME-$JOB_ID
mkdir -p $WDIR
cd $WDIR

# copy files, executables here
cp $DATADIR/$ALN .
cp $DATADIR/$TREE .
cp $SCRIPTDIR/* .
cp $SWEXEC swmutsel.jar
cp $HOME/bin/bin/HYPHYMP .


# run analysis
module load java
module load python
export PYTHONPATH=$HOME/bin/lib/python2.7/site-packages # to get biopython behaving w/ phylip-relaxed
python run_swmutsel.py $ALN $TREE $CPU ${DATA}_swmutsel $ESTMU

FINALDIR=$HOME/MutSel_benchmark/results/
mkdir -p $FINALDIR
cp *dnds.txt $FINALDIR
cp *MLE.txt $FINALDIR
cp *fitness.txt $FINALDIR
cp hyout.txt $FINALDIR/${DATA}_swmutsel_hyout.txt

cd
cp -r $WDIR /share/WilkeLab/work/sjs3495/mutsel_bench_data
rm -r $WDIR
