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
job_script_name="run_lrht_array.script"
lammps_inp_name="org_hoip.inp"

# makes sure dir has / at the end
dir=${dir%/}
dir+="/"

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
JOB_NAMES=$(printf "%s\n" "${JOB_NAMES_LIST[@]}" | sort -u)

## Start of Step 2
count=0
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

## Start the job
sbatch --export=NONE $dir$job_script_name
## --export=NONE fixes the module load issue
## Basically this avoids exporting the login node environment to the compute nodes
## since login has none of the modules
