# lammps_pre_post_processing
Pre and post processing scripts for running LAMMPS simulations through slurm. Scripts are divided into subdirectories. Most are focused specifically on Hybrid Organic Inorganic Perovskites (HOIPs).

## structure_editors 
 *contains scripts that directly create or edit MD input files.*

  -  **triple-halide_perovskite_builder.py** creates lammps data files of CsPb(I_{1-x-y}Br_{y}Cl_{x})_3 cubic perovskites. It is intended as a standalone script, so the percentages for Br and Cl, cell_parameters, element charges, and masses are written in variables. The charges on the swapped atoms affect the charges on the neighboring atoms, split by the amount of nearest neighbors. 
  -  **halogen_positioning_HOIP.py** edits a master .car .mdf pair to create all the mono and di halogenated placement combinations on the organic rings. The placements and halogen charges can be easily edited.
  -  **MBA_chirality_alteration.py** Takes in a .car .mdf pair to create 1) the opposite chirality and 2) no methyl group versions. 

## slurm_scripts 
*contains scripts that are used to prep / submit jobs while on the HPC.*

  -  **hoip_array.sh** takes a folder of .car .mdf pairs and preps and runs lammps simulations for all of them. It converts each to a .data file using msi2lmp.exe, creates a designated folder, and submits the job via slurm.
  -  **batch_job_arr.script** Slurm array batch script to submit a list of lammps simulations at once. Used in conjunction with hoip_array.sh.

## analysis_scripts
*contains post-processing files. Typically operate on log, .out, or trajectory files.*

  -  **last_avg_pproc_fold.py** takes a directory input and searches every subdirectory to create a csv from each lammps .out file, then creates a master csv that contains the last line and back half average from each individual .out it parses.
  -  **perov_struct_analysis.py** is a structural analysis script for HOIPs. Takes a lammps trajectory file (.dump) as input. Finds (among other things) lattice parameters, layer heights, organic tilts, and octahedral tilts. Outputs plots from averages and time series. 

## misc_HPC
*contains files that are not directly related to HOIPs, but are adjacently useful.*

  -  **accelrys_license_check.bat** Checks the current users of the Accelrys Materials Studio floating licenses. Ran in command prompt.
  -  **node_list.sh** Prints out the used and available cpus per node when given a partition. Will highlight your usage green. To run on slurm system, '<path_to_script>node_list.sh <partition>'
