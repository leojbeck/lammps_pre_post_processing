# Folder level Post Processing script
# Leo Beck
# Created: May 8th, 2024

# This script takes a folder name / directory, then searches every subdirectory & file within
# It then writes a csv for each output file, then writes the last line of each to a master csv
# It also writes the back half average to a separate master csv (WORK IN PROGRESS)

# To run, use the following in the command line:
# python last_avg_pproc_fold.py <folder> <(optional) number>

# <(optional) number>: 0 for printing master file to cwd, 1 for everything to cwd

import csv
import sys
import os
import numpy as np


# Check if the filename is provided as a command line argument
if len(sys.argv) < 2:
	print("Usage: python script_name.py <folder>")
	sys.exit(1)

# Get the LAMMPS output filename from command line argument
#lammps_output_file = sys.argv[1]
folder = sys.argv[1]
pfold = os.path.abspath(folder)
curfold = os.path.basename(os.path.normpath(pfold))
print(curfold)

# Change master csv print location to cwd
#if len(sys.argv) > 3:
    #if ((sys.argv[3] == '0') or (sys.argv[3] == '1')):
        #pfold = os.path.abspath(os.getcwd())



# Initialize variables
data = []
ll_data = []
avg_data = []
avg_data_tot = []
start_index = None
end_index = 0
es = 0

file_names = []
sys_names = []


# Walks through all the subdirectories, makes a list of file names for later
# This is the big part that makes this work in the parent directory - very cool
for subdir, dirs, files in os.walk(folder):
	for file in files:
		# Check if it is an output file
		if file.endswith('.out'):
			# Put subdir + file to get the filepath 
			file_names.append(os.path.join(subdir, file))


# Iterate through the files in the list of names
for lammps_output_file in file_names:
	# Open the LAMMPS output file
	with open(lammps_output_file, 'r') as file:
		lines = file.readlines()
		endline = len(lines)

	# Gather list of system names from the end of the first line in the output
	sys_names.append(open(lammps_output_file).readline().split()[-1])

	# Not in the regular pproc - Used for multiple iterations of step / loop runs
	num_runs = len([i for j,i in enumerate(lines) if i.strip().startswith('Step')])
	print('\n')
	print(lammps_output_file)
	print('Number of runs in file: ' + str(num_runs))
	
	# Iterate over output runs - normally 1 time, you'll know if you're making more
	# Idea to add this from Jordan's sims
	start_index = min([j for j,i in enumerate(lines) if i.strip().startswith('Step')])
	end_index = start_index + 1
	for k in range(0, num_runs):
		# Find the start and end indices for the data section
		# Starts at the previous end
		
		for i, line in enumerate(lines[end_index+1:endline], end_index+1):
			end_check = [line.strip().startswith('Loop'), line.strip().startswith('srun'), line.strip().startswith('slurmstepd:'), line.strip().startswith('WARNING:'), line.strip().startswith('ERROR ')]
			# Check line start (only the line start has capitalized 'Step'
			if line.strip().startswith('Step'):
				start_index = i
				print('Start index ' + str(k) + ': ' + str(i))
			elif  sum(end_check) > 0:
				end_index = i
				print('End index ' + str(k) + ': ' + str(i))
				break
			
			if i == endline:
				end_index = i + 1

		# Extract data for each time step	
		if ((end_index - start_index) > 3):
			for j in range(start_index, end_index - 1):
				line_data = lines[j].strip().split()
				data.append(line_data)
				# Check if it is the last line
				if j == end_index - 2:
					# save the last line`
					last_data = line_data
					# Append file name to the line
					last_data.append(lammps_output_file)
					last_data.append("L")
					last_data.append(sys_names[-1])
					ll_data.append(last_data)
					#print("L")


		# Write data to CSV file
		# Change individual csv print location to cwd
		#if ((sys.argv[2] == 1) or (sys.argv[2] == 2)):
			#pifold = os.path.abspath(os.getcwd())
		output_filename = lammps_output_file.replace('.out', '_'+str(k)+'.csv')
		# Checks if you want to print everything to cwd, then changes filepath if you do
		if sys.argv[2] == '1':
			output_filename = pfold + '/' + os.path.basename(output_filename)
		with open(output_filename, 'w') as csv_file:
			csv_writer = csv.writer(csv_file)
			csv_writer.writerows(data)

		# Get the back half averages, append to big data
		# Actually using last quarter rn
		if len(data) > 0:
			# Get midpoint length (to average back half)
			mid = len(data) // 2
			print("Mid: ", mid)
			# Transform list of lists into floats
			data = data[mid:-1]
			data_tonum = []
			for lst in data:
				data_tonum.append([float(i) for i in lst])
			# Get average of each column of floats
			avg_data = np.mean(data_tonum, axis = 0)
			# Turn back into list of strings
			avg_data = list(avg_data)
			avg_data = [str(i) for i in avg_data]
			# append file info to end
			avg_data.append(lammps_output_file)
			avg_data.append("A")
			avg_data.append(sys_names[-1])
			avg_data_tot.append(avg_data)
			#print(avg_data)

		# Clear data so it can write individually
		data = []

		print('Data extracted and saved to ' + output_filename)

# If statement to check if the folder starts with a period (e.g. '.', '../fold1/')
#ll_filename = folder + folder.replace('/', '') + 'mast_output.csv'
ll_filename = pfold + '/' + curfold + 'mast_output.csv'

with open(ll_filename, 'w') as csv_file:
	csv_writer = csv.writer(csv_file)
	csv_writer.writerows(ll_data)
	csv_writer.writerows(avg_data_tot)
	

print('Master output sent to ' + ll_filename)
