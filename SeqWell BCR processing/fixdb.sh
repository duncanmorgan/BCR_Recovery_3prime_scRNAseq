#!/bin/sh
# request Bourne shell as shell for job
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m be
#$ -pe whole_nodes 1
#SBATCH --exclude=c[10]

###################  
module load python3
module load igblast
module load usearch
module load muscle
module load fastxtoolkit
module load presto
module load changeo

cd $5
echo $5

MakeDb.py igblast -i file_igblast.fmt7 -s file.fa -r ~/data/changeodb/germlines/imgt/$4/vdj/ --log mdb.log --extended --failed --partial

#ParseLog.py -l BC1.log -f ID SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
#ParseLog.py -l BC2.log -f ID SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
#ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE

python3 ~/data/211108Lov/Rh22/compile_logs.py
