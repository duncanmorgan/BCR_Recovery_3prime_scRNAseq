import pandas as pd
import numpy as np
import glob

files = glob.glob('*/*final_call*')

for curr in files:
	print(curr)
	sample_name = curr.split('/')[0]
	curr_df = pd.read_csv(curr, sep = '\t')
	curr_df['Sample'] = sample_name
	curr_df['ID'] = [curr_df.Sample[x] + curr_df.SEQUENCE_ID[x] for x in range(0, curr_df.shape[0])]
	curr_df.index = curr_df.ID
	

	
	if curr == files[0]:
		df = curr_df
	else:
		df = pd.concat([df, curr_df])


print(df.NP1_LENGTH.head())
df.to_csv('combined_results.tab', sep = '\t', index = False)
