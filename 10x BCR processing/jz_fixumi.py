import pandas as pd
import umi_tools
import numpy as np
from umi_tools import UMIClusterer
from collections import Counter
import time
start_time = time.time()
print(f'Starting time: {start_time/60}')

umis = np.loadtxt('BC.txt', dtype=bytes)
read_time = time.time()
print(f'Read in BC.txt used {(read_time-start_time)/60} minutes')

clusterer = UMIClusterer(cluster_method='directional')

umi_counts = Counter(umis)
print(len(umi_counts))

clustered_umis = clusterer(umi_counts, threshold=1)
cluster_time = time.time()
print(f'Clustering used {(cluster_time-read_time)/60} minutes')

# Create a dictionary using set operations
correction_dict = {each: correct for correct, *others in clustered_umis for each in others}
correction_dict.update({correct: correct for correct, *_ in clustered_umis})

# Use list comprehension to create corrected_umis
corrected_umis = [correction_dict[umi] for umi in umis]

corrected_counts = Counter(corrected_umis)
print(len(corrected_counts))

corrected_umis = [x.decode('utf-8') for x in corrected_umis]
fix_time = time.time()
print(f'Fixing UMI used {(fix_time-cluster_time)/60} minutes')

np.savetxt('BC_fixed.txt', corrected_umis, fmt='%s')

end_time = time.time()
elapsed_time = (end_time - start_time)
print(f"Elapsed time: {elapsed_time/60} minutes")
