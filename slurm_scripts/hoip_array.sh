#!/bin/bash
## alternative to sbatch --export=NONE: '#!/bin/bash -i'
## the -i fixes a module load issue by making this script 'interactive'

## Leo Beck - hois_array.sh
## February 11, 2025

## This bash script takes a folder of car and mdf files, then preps everything and 
## starts the simulations. Useful for many structures with the same simulation parameters
## To run, type './bash_file.sh <folder>/'

## Specific steps:
##	1) initialize variables, then create a list of names. This is done by checking the name of
##		every file in the <folder> (dir=$1), removing the extension, then creating a list of 
##		the unique values
##	2) Iterate over each filename
##		a) create .data using the .car and .mdf (the ./msi2lmp.exe line)
##		b) Create folder for that filename, then move the .car, .mdf, and .data in
##		c) Move the files into the filename folder
##	3) Start the job using slurm array job script

## Exits script if any statement returns non-true return (prevents error snowballing)
set -e 

## Used to create more descriptive names (e.g. filename_aniso_5M)
## (Not utilized at the moment)
postroot=""

## Initialize variables (dir=$1 grabs <folder>)
JOB_NAMES_LIST=()
JOB_NAMES=()
dir=$1
temp=""
job_script_name=(*job*.script)
echo "Job script: $job_script_name"
##lammps_inp_name=(*.inp)
frc_name=(*.frc)
echo "Forcefield file: $frc_name"
frc_name=${frc_name%.*}

# makes sure dir has / at the end
dir=${dir%/}
dir+="/"

## Make sure msi2lmp.exe is executable
chmod 777 msi2lmp.exe || echo "msi2lmp.exe not found."


## First for loop - creates list of unique filenames (without extensions)
for f in "$dir"*; do
  ## Check if file ends with car or mdf
  if [[ $f == *".car" ]] || [[ $f == *".mdf" ]]
  then  
    ## Get rid of leading filepath (if it exists)
    temp=${f##*/}
    ## Get rid of file extension
    temp=${temp%.*}
    echo ${temp}
    ##temp=$(basename "$f" | sed 's\/.[^.]*$//')
    ## Add filename to list of names
    JOB_NAMES_LIST+=($temp)
  else
    echo "Skipping $f"
  fi
done

## The part that creates the list of unique names
##JOB_NAMES=$(printf "%s\n" "${JOB_NAMES_LIST[@]}" | sort -u)
mapfile -t JOB_NAMES < <(printf "%s\n" "${JOB_NAMES_LIST[@]}" | sort -u)

for str in ${JOB_NAMES[@]}; do
  name=$str$postroot

  ## do the msi thing
  ## Currently has the frc in the same folder as bash_file.sh (hence $PWD).
  ## Change to your liking
  ./msi2lmp.exe $dir$str -c 1 -p 0 -frc $PWD/cvff_iff_v1_5_G_poly_solv_ions_per_v55 -i

  ## make and move the files
  mkdir $dir$str
  mv $dir$str.* $dir$str/.
done

## Edit the job array line to match number of subdirs
## Find max index of job array (0 to n-1)
NUM_JOBS=${#JOB_NAMES[@]}
MAX_IDX=$((NUM_JOBS - 1))

## Swap line starting with #SBATCH --array=
sed -i "s/^#SBATCH --array=.*/#SBATCH --array=0-${MAX_IDX}/" $job_script_name

## Start the job
sbatch --export=NONE $job_script_name
## --export=NONE fixes the module load issue
## Basically this avoids exporting the login node environment to the compute nodes since login has none of the modules
