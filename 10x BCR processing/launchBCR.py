import os
import glob
import shutil

files = glob.glob('*.fastq')
samples = list(set([x.split('_')[0] for x in files]))
samples.sort()
species = 'human' 

for curr in samples:
	# make directory
	os.makedirs(curr, exist_ok = True)
	curr_files = glob.glob(curr + '*.fastq')
	curr_files.sort()
	print(curr)

	# copy fastq files
	for each in curr_files:
		os.system('cp ' + each + ' ' + curr + '/')

	# copy script
	os.system('cp bcrprocess.sh ' + curr)

	# launch job
	os.system('sbatch bcrprocess.sh ' + curr_files[0] + ' ' + curr_files[1] + ' ' + curr_files[2] + ' ' + species  + ' '  + os.getcwd() + '/' + curr)
	
