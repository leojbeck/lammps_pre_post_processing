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
frc_file=(*.frc)

## Check that job script and frc are grabbed correctly
if [[ "$job_script_name" == "*job*.script" ]] ; then
    echo "Batch job script file not found."
    exit 1
elif (( ${#job_script_name[@]} > 1 )); then
    echo "Multiple job scripts found."
    exit 1
else
    echo "Batch job script file: $job_script_name"
fi

if [[ "$frc_file" == "*.frc" ]] ; then
    echo "Forcefield file not found."
    exit 1
elif (( ${#frc_file[@]} > 1 )); then
    echo "Error: Multiple forcefield files found: ${frc_file[*]}"
    exit 1
else 
    echo "Forcefield file: $frc_file"
fi
frc_file=${frc_file%.*}

# makes sure dir has / at the end
dir=${dir%/}
dir+="/"

## Make sure msi2lmp.exe is executable
chmod 777 msi2lmp.exe || (echo "msi2lmp.exe not found."; exit 2)


## First for loop - creates list of unique filenames (without extensions)
for f in "$dir"*; do
  ## Check if file ends with car or mdf
  if [[ $f == *".car" ]] || [[ $f == *".mdf" ]]
  then  
    ## Get rid of leading filepath (if it exists)
    temp=${f##*/}
    # echo ${temp}
    ## Get rid of file extension
    temp=${temp%.*}
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
  echo "Prepping ${str}"
  # Check to make sure both car and mdf are in folder.
  [ ! -f $dir$str.car ] && echo "${str}.car file not found!"
  [ ! -f $dir$str.mdf ] && echo "${str}.mdf file not found!"
  
  ## Currently has the frc in the same folder as bash_file.sh (hence $PWD).
  ./msi2lmp.exe $dir$str -c 1 -p 0 -frc $PWD/$frc_file -i || (echo "msi2lmp failed for $str."; exit 1)

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
