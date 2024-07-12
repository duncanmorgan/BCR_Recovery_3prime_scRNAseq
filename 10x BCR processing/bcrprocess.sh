#!/bin/sh
# request Bourne shell as shell for job
#$ -S /bin/sh
#$ -cwd
#$ -n 16
#$ -V
#$ -m be
#$ -pe whole_nodes 1
###################  

module load python3
module load igblast
module load usearch
module load muscle
module load fastxtoolkit
module load presto
module load changeo

cd $5

# print out which folder is being processed 
pwd 

# parse the fastq file and extrapolate sequences (line2) and their corresponding quality scores (line4)
# do it on read 1 (cell barcode + UMI=28 bp for 10x Genomics, reads from 3' end)
awk 'NR%4==2' $1 | cut -c -28 > BC.txt
awk 'NR%4==0' $1 | cut -c -28 > BC_Q.txt

# do the same for Read 3 (Nextera primer, start at V region towards D, J, C)
awk 'NR%4==2' $3 > R2.txt
awk 'NR%4==0' $3 > R2_Q.txt

# Hamming distance correction on UMIs
python3 ../jz_fixumi.py
paste -d '' BC_fixed.txt R2.txt > BCR2.txt
paste -d '' BC_Q.txt R2_Q.txt > BCR2_Q.txt

# do the same for Read 2 (start at constant region towards J,D,V)
awk 'NR%4==2' $2 > R1.txt
awk 'NR%4==0' $2 > R1_Q.txt

# Reconstruct the fastq file for each sample 
paste -d '' BC_fixed.txt R1.txt > BCR1.txt
paste -d '' BC_Q.txt R1_Q.txt > BCR1_Q.txt

awk 'NR%4==1' $3 > R2_header.txt
awk 'NR%4==3' $3 > R2_Qheader.txt
 
awk 'NR%4==1' $2 > R1_header.txt
awk 'NR%4==3' $2 > R1_Qheader.txt
 
paste -d '\n' R2_header.txt BCR2.txt R2_Qheader.txt BCR2_Q.txt > BC_R2.fastq
paste -d '\n' R1_header.txt BCR1.txt R1_Qheader.txt BCR1_Q.txt > BC_R1.fastq

# Trim R1 and R2 down to 250 bp (this number depends on Phred scores, often drop below 30 at 250bp) 
# output is trimmed R1 and R2
fastx_trimmer -l 250 -i BC_R1.fastq -o BC_R1_trim.fastq -f 1 -Q33  
fastx_trimmer -l 250 -i BC_R2.fastq -o BC_R2_trim.fastq -f 1 -Q33  

# Filter R1 and R2 based on the quality scores; -q specifies the quality score cut off
FilterSeq.py quality -s BC_R1_trim.fastq -q 25 --outname Filtered_R1 --log FS1.log --failed
FilterSeq.py quality -s BC_R2_trim.fastq -q 25 --outname Filtered_R2 --log FS2.log --failed

# Match and mask the given primers 
MaskPrimers.py score -s  Filtered_R1_quality-pass.fastq --outname R1PRIMER --failed --mode mask --barcode --start 28 -p ../vh_IG_primer.fasta --log MP1.log --maxerror .4
MaskPrimers.py score -s  Filtered_R2_quality-pass.fastq --outname R2PRIMER --failed --mode mask --start 28 -p ../vh_primers.fasta --log MP2.log --barcode --maxerror .4

ParseHeaders.py rename -s R2PRIMER_primers-pass.fastq -f PRIMER -k VPRIMER
ParseHeaders.py rename -s R1PRIMER_primers-pass.fastq -f PRIMER -k ISOTYPE

# Sorts and matches sequence records with matching coordinates across files. If a molecule is in R2 but not in R1, count as fail.  
PairSeq.py -1 R1PRIMER_primers-pass_reheader.fastq -2 R2PRIMER_primers-pass_reheader.fastq --2f VPRIMER --coord illumina --1f ISOTYPE --failed

# Cluster the sequences 
ClusterSets.py set -s R1PRIMER_primers-pass_reheader_pair-pass.fastq --outname R1 -k CLUSTER -f BARCODE --log CS1.log
ClusterSets.py set -s R2PRIMER_primers-pass_reheader_pair-pass.fastq --outname R2 -k CLUSTER -f BARCODE --log CS2.log

# parse the header to combine cluster and barcode 
ParseHeaders.py merge -s R1_cluster-pass.fastq -f BARCODE CLUSTER -k BARCODE_CLUSTER --act cat -o R1-cluster-pass-reheader.fastq 
ParseHeaders.py merge -s R2_cluster-pass.fastq -f BARCODE CLUSTER -k BARCODE_CLUSTER --act cat -o R2-cluster-pass-reheader.fastq

# build consensus among the BARCODE_CLUSTER copy ISOTYPE and VPRIMER over 
BuildConsensus.py -s R1-cluster-pass-reheader.fastq --bf BARCODE_CLUSTER --pf ISOTYPE --cf ISOTYPE VPRIMER --act set set -n 1 --maxerror 0.50 --maxgap 0.5 --outname R1-cluster --log BC1.log
BuildConsensus.py -s R2-cluster-pass-reheader.fastq --bf BARCODE_CLUSTER --pf VPRIMER --cf ISOTYPE VPRIMER --act set set -n 1 --maxerror 0.50 --maxgap 0.5 --outname R2-cluster --log BC2.log

# calculate and fileter clusters in one step
python3 ../filter_cluster.py R1-cluster_consensus-pass.fastq R1-filtered-cluster-consensus.fastq
python3 ../filter_cluster.py R2-cluster_consensus-pass.fastq R2-filtered-cluster-consensus.fastq

PairSeq.py -1 R1-filtered-cluster-consensus.fastq -2 R2-filtered-cluster-consensus.fastq --coord presto --failed

AssemblePairs.py align -1 R1-filtered-cluster-consensus_pair-pass.fastq -2 R2-filtered-cluster-consensus_pair-pass.fastq --log AP.log  --coord presto --rc tail --outname pairs

paste - - - - < pairs_assemble-pass.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > file.fa

AssignGenes.py igblast -s file.fa --organism $4 --loci ig --format blast -b ../igblast 

MakeDb.py igblast -i file_igblast.fmt7 -s file.fa --log mdb.log --extended --failed --partial -r ../germlines/imgt/$4/vdj/ 

ParseLog.py -l BC1.log -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
ParseLog.py -l BC2.log -f ID BARCODE SEQCOUNT PRIMER PRCOUNT PRCONS PRFREQ CONSCOUNT
ParseLog.py -l AP.log -f ID LENGTH OVERLAP ERROR PVALUE

# lastly compile the logs to generate "final_calls.csv"
python3 ../compile_logs.py
