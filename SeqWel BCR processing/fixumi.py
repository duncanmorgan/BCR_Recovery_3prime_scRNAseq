import pandas as pd
import umi_tools
import numpy as np


umis = np.loadtxt('BC.txt', dtype = bytes)

from umi_tools import UMIClusterer
clusterer = UMIClusterer(cluster_method = 'directional')

from collections import Counter
umi_counts = Counter(umis)
print(len(umi_counts))

clustered_umis = clusterer(umi_counts, threshold = 1)
corrected_umis = umis.copy()
print('clustering complete, fixing UMIs')


corrected_umis = umis.copy()

dict = {}
for x in clustered_umis:
	correct = x[0]
	if (len(x) > 1):
		for each in x[1:]:
			dict[each] = correct
	dict[correct] = correct
	

for i in range(0, len(corrected_umis)):

	corrected_umis[i] = dict[umis[i]]
	
corrected_counts = Counter(corrected_umis)
print(len(corrected_counts))

corrected_umis = [x.decode('utf-8') for x in corrected_umis]

np.savetxt('BC_fixed.txt', corrected_umis, fmt = '%s')