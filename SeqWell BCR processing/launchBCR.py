import os
import glob
import shutil

files = glob.glob('*.fastq')
samples = list(set([x.split('_')[0] for x in files]))

for curr in samples:
	# make directory
	os.makedirs(curr, exist_ok = True)
	curr_files = glob.glob(curr + '*.fastq')
	curr_files.sort()
	# copy fastq files
	for each in curr_files:
		os.system('cp ' + each + ' ' + curr + '/')
		
	# copy fixumi
	os.system('cp ~/data/fixumi.py ' + curr)
	
	# copy script
	os.system('cp bcrprocess.sh ' + curr)
	
	# launch job
	os.system('sbatch bcrprocess.sh ' + curr_files[0] + ' ' + curr_files[1] + ' ' + curr_files[2] + ' human ' + os.getcwd() + '/' + curr)
	