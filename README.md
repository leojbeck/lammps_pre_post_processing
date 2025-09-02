# lammps_pre_post_processing
Pre and post processing scripts for running LAMMPS simulations through slurm. Scripts are divided into subdirectories. Most are focused specifically on Hybrid Organic Inorganic Perovskites (HOIPs).

## structure_editors 
 *contains scripts that directly create or edit MD input files.*

  -  **triple-halid_perovskite_builder.py** creates lammps data files of CsPb(I_{1-x-y}Br_{y}Cl_{x})_3 cubic perovskites. It is intended as a standalone script, so the percentages for Br and Cl, cell_parameters, element charges, and masses are written in variables. The charges on the swapped atoms affect the charges on the neighboring atoms, split by the amount of nearest neighbors. 
  -  **halogen_positioning_HOIP.py** edits a master .car .mdf pair to create all the mono and di halogenated placement combinations on the organic rings. The placements and halogen charges can be easily edited.

## slurm_scripts 
*contains scripts that are used to prep / submit jobs while on the HPC.*

  -  **hoip_array.sh** takes a folder of .car .mdf pairs and preps and runs lammps simulations for all of them. It converts each to a .data file using msi2lmp.exe, creates a designated folder, and submits the job via slurm.
  -  **slurm_job_array.sh** is a job array batch script for slurm to run lammps simulations. 

## analysis_scripts
*contains post-processing files. Typically operate on log, .out, or trajectory files.*

  -  **last_avg_pproc_fold.py** takes a directory input and searches every subdirectory to create a csv from each lammps .out file, then creates a master csv that contains the last line and back half average from each individual .out it parses.
  -  **perov_struct_analysis.py** is a structural analysis script for HOIPs. Takes a lammps trajectory file (.dump) as input. Finds (among other things) lattice parameters, layer heights, organic tilts, and octahedral tilts. Outputs plots from averages and time series. 
