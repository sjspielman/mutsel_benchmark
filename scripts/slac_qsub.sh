#$ -N JOB
#$ -e e_$JOB_NAME
#$ -o o_$JOB_NAME
#$ -S /bin/bash
#$ -q wilke
#$ -m beas
#$ -pe serial 15
source ~/.bashrc

# define files, variables used in analysis
CPU=15
DATA=$1
TREE=$2
ALN=$DATA.fasta
OUTFILE=${DATA}_slac.txt

DATADIR=$HOME/sitewise_dnds_mutsel_exp/data
RESULTDIR=$HOME/sitewise_dnds_mutsel_exp/results
SCRIPTDIR=$HOME/sitewise_dnds_mutsel_exp/scripts/hyphy

# create and enter working directory
WDIR=/state/partition1/sjs3495/$JOB_NAME-$JOB_ID
mkdir -p $WDIR
cd $WDIR

# copy files, hyphy executable here
cp $DATADIR/alignments/$ALN temp.fasta
cp $DATADIR/trees/$TREE temp.tre
cp $SCRIPTDIR/autoSLAC.bf .
cp $SCRIPTDIR/parse_slac.py .
cp $SCRIPTDIR/replace_placeholder.py .
cp $HOME/bin/bin/HYPHYMP .

#FEL
python replace_placeholder.py autoSLAC.bf $WDIR
./HYPHYMP CPU=$CPU autoSLAC.bf

# Parse SLAC analysis
module load python
python parse_slac.py $OUTFILE
sed -i "s/ //g" $OUTFILE # remove spaces, so that R can read it in properly

mkdir -p $RESULTDIR
cp $OUTFILE $RESULTDIR


# Check if a previous run exists and if so remove it. Then save and cleanup.
#cd $RDIR
#ls | grep $JOB_NAME | xargs rm -r
#cp -r $WDIR .
rm -r $WDIR
