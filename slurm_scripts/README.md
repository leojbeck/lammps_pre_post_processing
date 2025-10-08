This folder, slurm_scripts, is used for prepping and submitting batch jobs on a slurm system. 

hoip_array.sh
This script takes a directory with .car and .mdf files, then creates subdirectories for each pair, converts them to .data files using msi2lmp.exe (in the working directory), then gets ready to submit a batch job array.

slurm_job_array.sh
This submits a job array via slurm for every subdirectory in the working directory. 

slurm_auto_rerun.sh
This submits a job that is expected to take longer than the wall time. When it approaches wall time, it cancels the running job and restarts the job. Needs to be paired with an input file that checks if the restart files exist (see lammps_auto_restart.inp)

lammps_auto_restart.inp
lammps input script that checks if the restart files exist. If so, it starts from the restart files. Otherwise, it will start from the data file. 
