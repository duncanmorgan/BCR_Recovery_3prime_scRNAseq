import pandas as pd
import os
import re

df = pd.read_csv('file_igblast_db-pass.tab', index_col = 0, sep = '\t')
df.rename(index={'SEQUENCE_ID':'BARCODE'}, inplace=True)

bc1_log = pd.read_csv('BC1_table.tab', index_col = 0, sep = '\t')
# remove the cluster number in the BARCODE column 
bc1_log['BARCODE'] = [re.split(r'\d', x, 1)[0] for x in bc1_log.BARCODE]
# Calculate the sum 
sumCount = pd.DataFrame(bc1_log.groupby(['BARCODE'])['CONSCOUNT'].sum())
sumCount.reset_index(inplace=True)
sumCount.rename(columns={'CONSCOUNT':'total_count'}, inplace=True)
# populate bc1_log
bc1_log = bc1_log.merge(sumCount,how='left', on='BARCODE')
# calculate fraction 
bc1_log['FRACTION'] = bc1_log['CONSCOUNT']/bc1_log['total_count']
# only keep the cluster with the greatest cluster fraction for each molecule 
bc1_log = bc1_log.loc[bc1_log.groupby(["BARCODE"])["FRACTION"].idxmax()]

# set index to BARCODE without cluster number 
bc1_log.set_index('BARCODE', inplace=True)
df['ISOTYPE'] = bc1_log.PRCONS[df.index]
df['ISOTYPE_FREQ'] = bc1_log.PRFREQ[df.index]
df['R1CONSCOUNT'] = bc1_log.CONSCOUNT[df.index]
df['R1_CLUSTER_FRACTION'] = bc1_log.FRACTION[df.index]

bc2_log = pd.read_csv('BC2_table.tab', index_col = 0, sep = '\t')
bc2_log['BARCODE'] = [re.split(r'\d', x, 1)[0] for x in bc2_log.BARCODE]
sumCount = pd.DataFrame(bc2_log.groupby(['BARCODE'])['CONSCOUNT'].sum())
sumCount.reset_index(inplace=True)
sumCount.rename(columns={'CONSCOUNT':'total_count'}, inplace=True)
bc2_log = bc2_log.merge(sumCount,how='left', on='BARCODE')
bc2_log['FRACTION'] = bc2_log['CONSCOUNT']/bc2_log['total_count']
bc2_log = bc2_log.loc[bc2_log.groupby(["BARCODE"])["FRACTION"].idxmax()]

bc2_log.set_index('BARCODE', inplace=True)
df['VPRIMER'] = bc2_log.PRCONS[df.index]
df['VPRIMER_FREQ'] = bc2_log.PRFREQ[df.index]
df['R2CONSCOUNT'] = bc2_log.CONSCOUNT[df.index]
df['R2_CLUSTER_FRACTION'] = bc2_log.FRACTION[df.index]

ap_log = pd.read_csv('AP_table.tab', index_col = 0, sep = '\t')
df['ERROR'] = ap_log.ERROR[df.index]
df['LENGTH'] = ap_log.LENGTH[df.index]
df['OVERLAP'] = ap_log.OVERLAP[df.index]
df['PVAL'] = ap_log.PVALUE[df.index]

outfile = os.getcwd().split('/')[-1]+'_final_calls.csv'
df.to_csv(outfile, sep = '\t')
print('compile_log.py Done!')
