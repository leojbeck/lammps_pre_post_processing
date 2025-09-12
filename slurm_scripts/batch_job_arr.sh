#!/bin/bash
#SBATCH --job-name=testrun
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=ucb553_asc1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
##SBATCH --exclusive
#SBATCH --output=testout_%A_%a.out
#SBATCH --mail-type=END,FAIL
##SBATCH --constraint=ib
#SBATCH --array=0-2

# LAMMPS input script name
inp_name=(*.inp)

ml purge
ml intel/2022.1.2
ml impi/2021.5.0
ml lammps/2Aug23

# Grab all subdirectories and store them in an array
dirs=($(ls -d */ | sort))         # Ensure consistent ordering

# Get subdirectory name for this array task
subdir=${dirs[$SLURM_ARRAY_TASK_ID]}
subdir=${subdir%/}                # Remove trailing slash
basename="${subdir}"
# Extract the base filename
datafile="${subdir}.data"

echo "Running LAMMPS in $subdir with data file $datafile"

# Change into the subdirectory
cd "$subdir"

# Call LAMMPS and pass the structure file as a variable
lmp -var datafile "$datafile" -var basename "$basename" -in ../$inp_name
