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

awk 'NR%4==2' $1 | cut -c -20 > BC.txt
python3 fixumi.py
awk 'NR%4==0' $1 | cut -c -20 > BC_Q.txt

awk 'NR%4==2' $3 > R2.txt
awk 'NR%4==0' $3 > R2_Q.txt

paste -d '' BC_fixed.txt R2.txt > BCR2.txt
paste -d '' BC_Q.txt R2_Q.txt > BCR2_Q.txt

awk 'NR%4==2' $2 > R1.txt
awk 'NR%4==0' $2 > R1_Q.txt

paste -d '' BC_fixed.txt R1.txt > BCR1.txt
paste -d '' BC_Q.txt R1_Q.txt > BCR1_Q.txt


awk 'NR%4==1' $3 > R2_header.txt
awk 'NR%4==3' $3 > R2_Qheader.txt
 
awk 'NR%4==1' $2 > R1_header.txt
awk 'NR%4==3' $2 > R1_Qheader.txt
 
paste -d '\n' R2_header.txt BCR2.txt R2_Qheader.txt BCR2_Q.txt > BC_R2.fastq
paste -d '\n' R1_header.txt BCR1.txt R1_Qheader.txt BCR1_Q.txt > BC_R1.fastq

fastx_trimmer -l 250 -i BC_R1.fastq -o BC_R1_trim.fastq -f 1 -Q33  
fastx_trimmer -l 250 -i BC_R2.fastq -o BC_R2_trim.fastq -f 1 -Q33  


FilterSeq.py quality -s BC_R1_trim.fastq -q 25 --outname Filtered_R1 --log FS1.log --failed
FilterSeq.py quality -s BC_R2_trim.fastq -q 25 --outname Filtered_R2 --log FS2.log --failed

MaskPrimers.py score -s  Filtered_R1_quality-pass.fastq --outname R1PRIMER --failed --mode mask --barcode --start 20 -p ../IG_primer.fasta --log MP1.log 
MaskPrimers.py score -s  Filtered_R2_quality-pass.fastq --outname R2PRIMER --failed --mode mask --start 20 -p ../vh_primers.fasta --log MP2.log --barcode --maxerror .4

ParseHeaders.py rename -s R2PRIMER_primers-pass.fastq -f PRIMER -k VPRIMER
ParseHeaders.py rename -s R1PRIMER_primers-pass.fastq -f PRIMER -k ISOTYPE

PairSeq.py -1 R1PRIMER_primers-pass_reheader.fastq -2 R2PRIMER_primers-pass_reheader.fastq --2f VPRIMER --coord illumina --1f ISOTYPE

BuildConsensus.py -s R1PRIMER_primers-pass_reheader_pair-pass.fastq --bf BARCODE --pf ISOTYPE -n 5 --maxerror 0.50 --maxgap 0.5 --outname r1P --log BC1.log
BuildConsensus.py -s R2PRIMER_primers-pass_reheader_pair-pass.fastq --bf BARCODE -n 5 --pf VPRIMER  --maxerror 0.50 --maxgap 0.5 --outname r2P --log BC2.log

PairSeq.py -1 r2P_consensus-pass.fastq -2 r1P_consensus-pass.fastq --coord presto --failed

AssemblePairs.py align -1 r2P_consensus-pass_pair-pass.fastq -2 r1P_consensus-pass_pair-pass.fastq --log AP.log  --coord presto --rc tail --outname pairs

paste - - - - < pairs_assemble-pass.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > file.fa
AssignGenes.py igblast -s file.fa --organism $4 --loci ig --format blast -b ~/data/changeodb
MakeDb.py igblast -i file_igblast.fmt7 -s file.fa -r ~/share/germlines/imgt/$4/vdj/ --log mdb.log --extended --failed --partial

ParseLog.py -l BC1.log -f ID SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
ParseLog.py -l BC2.log -f ID SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE

python3 compile_logs.py
