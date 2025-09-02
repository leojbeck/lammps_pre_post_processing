This folder, slurm_scripts, is used for prepping and submitting batch jobs on a slurm system. 

hoip_array.sh
This script takes a directory with .car and .mdf files, then creates subdirectories for each pair, converts them to .data files using msi2lmp.exe (in the working directory), then gets ready to submit a batch job array.

slurm_job_array.sh
This submit a job array via slurm for every subdirectory in the working directory. 
