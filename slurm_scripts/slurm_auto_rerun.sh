#!/bin/bash
#SBATCH --job-name=auto_restart_job     # Job name
#SBATCH --partition=amilan              # Partition (atesting, amilan, atesting_a100, ....)
#SBATCH --qos=normal                    # Quality of service (normal, testing, long, condo) - only set long or condo for Blanca
##SBATCH --account=<acc_name>            # Account name (Ascent allocation ucb553_asc1, normal)
#SBATCH --output=job_output_%j.txt      # Output file (%j includes job_id)
#SBATCH --time=00:02:00                 # Wall time (HH:MM:SS)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # CPUs per task
##SBATCH --mem=16GB                      # Memory per node (Don't usually need to clarify)
##SBATCH --exclusive
##SBATCH --mail-type=END,FAIL
##SBATCH --constraint=ib
#SBATCH --requeue                       # Allow job requeuing 

# Define input and output filename, etc
inp_name=(*.inp)
data_files=(*.data)
output_file=$(awk -F'=' '/#SBATCH --output=/ {print $2}' < "${BASH_SOURCE[0]}" | awk '{print $1}')
#output_file=$(grep -oP 'output=\K[^ ]+' "${BASH_SOURCE[0]}")
echo ${output_file}
#output_file="job_output.txt"

foss=20 # factor of safety (seconds) - how many seconds before the walltime you want to timeout

# Load necessary modules
ml purge
ml intel/2022.1.2
ml impi/2021.5.0
ml lammps/2Aug23
##export OMP_NUM_THREADS=1 # with MPI, you just want one thread per task

# Define rerun time: 20 seconds before walltime
walltime_seconds=$(scontrol show job $SLURM_JOB_ID | grep -oP 'TimeLimit=\K[0-9:]+' | awk -F: '{print ($1 * 3600) + ($2 * 60) + $3}')
rt=$((walltime_seconds - $foss)) # Timeout set to 20 seconds before the walltime

# Timeout execution with polling
timeout ${rt}s mpirun -np $SLURM_NTASKS lmp_mpi -var datafile "$datafile" -var basename "$basename" -in $inp_name >> ${output_file} 2>&1
job_pid=$!  # Capture the process ID
wait $job_pid # Wait for an exit code from <job_pid>
exit_code=$?  # Capture the exit code of `wait`

# Check if the exit code is the one produced from the timeout
if [[ $exit_code == 124 ]]; then
    echo "Job timed out, requeuing..." >> ${output_file}
    #scontrol requeue $SLURM_JOB_ID
    scancel $job_pid
    sbatch $0 >> $output_file 2>&1
    exit 0  # Ensure the current job exits cleanly
fi

echo "Job completed successfully with exit code ${exit_code}." >> ${output_file}
