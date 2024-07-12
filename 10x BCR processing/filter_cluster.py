# Created by: Jason Zhang
# Last modified: March 30 2023
##################################################
# usage:
# python3 filter_cluster R1-cluster_consensus-pass.fastq R1-filtered-cluster-consensus.fastq
#           argv[0]     |         argv[1]               |      argv[2]
##################################################

from Bio import SeqIO
import pandas as pd
import sys
import re

ifile_name=sys.argv[1]
ofile_name=sys.argv[2]

# read in the fastq file 
records = list(SeqIO.parse(ifile_name,"fastq"))
print("\nCalculating CLUSTER_FRACTION for", len(records),"records in file:",ifile_name,'\n')

# parse out the info from id  
df=pd.DataFrame(index=[records[i].id for i in range(0, len(records))])
df['BARCODE'] = [re.split(r'\d', x, 1)[0] for x in df.index]
df['CONSCOUNT'] = [int(x.split('CONSCOUNT=')[1].split('|')[0]) for x in df.index]


# get the total count for each barcode  
barcode_total_count = pd.DataFrame(df.groupby(['BARCODE'])['CONSCOUNT'].sum(), index=None).reset_index()
barcode_total_count.rename(columns={"CONSCOUNT": "total_count"}, inplace=True)

# save the index (aka record.id as a column), otherwise the index will be removed after merging 
df.reset_index(inplace=True)

# merge the barcode count onto the df 
df = df.merge(barcode_total_count, how='left',on='BARCODE' )

# put the index back 
df.set_index('index', inplace=True)

# calculate the fraction of each cluster 
df['CLUSTER_FRACTION'] = df['CONSCOUNT']/df['total_count']

# only keep the cluster with the greatest cluster fraction for each molecule 
df = df.loc[df.groupby(["BARCODE"])["CLUSTER_FRACTION"].idxmax()]

# Open the input and output file handles
with open(ifile_name, "r") as input_handle, open(ofile_name, "w") as output_handle:
    # keep track of number of records kept 
    n_output = 0

    for i,record in enumerate(SeqIO.parse(input_handle, "fastq")):
    
        # parce out BARCODE info 
        info = record.id

        if info in df.index:
            n_output += 1 
            # first remove the number after barcode
            line= record.id.split('|')
            line[0] = ''.join([i for i in line[0] if not i.isdigit()]) 
            record.id='|'.join(line)
            # update the record description with cluster_fraction info 
            record.description = record.id
            SeqIO.write(record, output_handle, "fastq")
            

print('Kept',n_output,"out of",i, "records to file:",ofile_name,'\n')
