import pandas as pd

df = pd.read_csv('file_igblast_db-pass.tab', index_col = 0, sep = '\t')

bc1_log = pd.read_csv('BC1_table.tab', index_col = 0, sep = '\t')
df['ISOTYPE'] = bc1_log.PRCONS[df.index]
df['ISOTYPE_FREQ'] = bc1_log.PRFREQ[df.index]
df['R1CONSCOUNT'] = bc1_log.CONSCOUNT[df.index]
df['R1SEQCOUNT'] = bc1_log.SEQCOUNT[df.index]


bc2_log = pd.read_csv('BC2_table.tab', index_col = 0, sep = '\t')
df['VPRIMER'] = bc2_log.PRCONS[df.index]
df['VPRIMER_FREQ'] = bc2_log.PRFREQ[df.index]
df['R2CONSCOUNT'] = bc2_log.CONSCOUNT[df.index]
df['R2SEQCOUNT'] = bc2_log.SEQCOUNT[df.index]

ap_log = pd.read_csv('AP_table.tab', index_col = 0, sep = '\t')
df['ERROR'] = ap_log.ERROR[df.index]
df['LENGTH'] = ap_log.LENGTH[df.index]
df['OVERLAP'] = ap_log.OVERLAP[df.index]
df['PVAL'] = ap_log.PVALUE[df.index]

df.to_csv('final_calls.csv', sep = '\t')