# lammps_pre_post_processing
Pre and post processing scripts for running lammps simulations through Slurm. Scripts are divided into subdirectories. Most are focused specifically on Hybrid Organic Inorganic Perovskites (HOIPs).

structure_editors/ contains scripts that directly create or edit input files. 

  -  triple-halid_perovskite_builder.py creates lammps data files of CsPb(I_{1-x-y}Br_{y}Cl_{x})_3 perovskites.
  -  halogen_positioning_HOIP.py edits a master .car .mdf pair to create all the mono and di halogenated placement combinations on the organic rings. 

slurm_scripts/ contains scripts that are used to prep / submit jobs while on the HPC.

  -  hoip_array.sh takes a folder of .car .mdf pairs and preps and runs lammps simulations for all of them. It converts each to a .data file using msi2lmp.exe, creates a designated folder, and submits the job via slurm.
  -  slurm_job_array.sh is a job array batch script for slurm to run lammps simulations. 

analysis_scripts/ contains post-processing files.

  -  last_avg_pproc_fold.py takes a directory input and searches every subdirectory to create a csv from each lammps .out file, then creates a master csv that contains the last line and back half average from each individual .out it parses.
  -  perov_struct_analysis.py is a structural analysis script for HOIPs. Takes a lammps trajectory file (.dump) as input. Finds (among other things) lattice parameters, layer heights, organic tilts, and octahedral tilts. Outputs plots from averages and time series. 
