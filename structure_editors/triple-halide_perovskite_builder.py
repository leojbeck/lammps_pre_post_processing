# Triple-halide perovskite structure building
#            Leo Beck
#           2024-09-10
# This python script uses Atomic Simulation Environment (ASE)
# to build CsPbI3 perovskite structures, then convert some I to Br,
# then some more I to Cl, changing halogen and neighbor charges
# (goal structure: CsPb(I_{1-x-y}Br_{y}Cl_{x})_3

# This script (at time of creation) is fairly one-off, meant for 
# creating the above mentioned structures to run stability 
# simulations for a little joint project between Heinz and McGehee

# ______________________________ HOW TO USE ______________________________
"""
1) Define the percentage of Iodine atoms you want to change into Br & Cl (pBr, pCl)
2) Check / set the cell parameters. a=b=c=6.289 is for cubic CsPbI3
3) Define supercell size [#x, #y, #z]
4) List the forcefield filename and output name (with relative path)
5) Double check forcefield types, charges, and masses in the elem_props dictionary
6) Run!

"""

# Import packages
from ase import Atoms
from ase import build
from ase.visualize import view
from ase.neighborlist import NewPrimitiveNeighborList
from ase.neighborlist import NeighborList
from ase.io.lammpsdata import write_lammps_data

import numpy as np
from random import randrange
from random import sample as rv_sample

import matplotlib.pyplot as plt
from matplotlib import cm
import os

# ______________________________ USER INPUT SECTION ______________________________
# Define size, percentages, et c.
pBr = 0.30
pCl = 0.1
cell_par = [6.289, 6.289, 6.289, 90, 90, 90]
base_perov_elm = 'CsPbI3'
supercell_ABC  = np.array([3, 3, 3])
frc_filename = "cvff_iff_v1_5_G_poly_solv_ions_per_v55.frc"
output_filename = "CsPbI"+str(int(100*round((1-pBr-pCl),2)))+"Br"+str(int(100*pBr))+"Cl"+str(int(100*pCl))+".data"
#output_filename = "LUNAR-main"+os.path.sep+"test.data"

# Define dictionary of element properties: element, charge, mass
# Put in alphabetical order please! :)
elem_props = {
  "Br": ["Br", -0.88, 79.904],
  "Cl": ["Cl", -0.9, 35.453],
  "Cs": ["cs+", 1.0, 132.91],
  "I": ["Ipb1", -0.84, 126.9044],
  "Pb": ["pb2+", 1.52, 207.2]
}

# True / False for plotting charge density
want_to_plot = False
# True / False for writing nta file
write_nta = False
# ______________________________ ______________________________

# Other stuff that doesn't need to change as often
a = cell_par[0]
b = cell_par[1]
c = cell_par[2]
perov_locs = [(a/2, a/2, a/2), (0, 0, 0), (a/2, 0, 0), (0, a/2, 0), (0, 0, a/2)] 

# Function to make cleaner code later
def swap_atom(pre, post, list, elch, charges, supercell):
    for count, rv in enumerate(list):
        #print(rv)
        # Swap and change charge
        supercell[rv].symbol = elem_props[post][0]
        charges[rv] = elem_props[post][1]
        el_mass[rv] = elem_props[post][2]
        forces[rv]  = elem_props[post][0]
        # Find nearest neighbors, update charge
        # Find the nearest neighbors of the swapped atom
        indices, offsets = nl.get_neighbors(rv)
        
        # Distribute the charge difference between I and Br (-0.88 - (-0.84) = -0.04)
        # Split the difference among the neighbors
        charge_difference = (elch[pre][1] - elch[post][1]) / len(indices)

        # Update nearest neighbor charge
        for neighbor_index in indices:
            charges[neighbor_index] += charge_difference

# Function to map charges to colors using a colormap
def map_charges_to_colors(charges):
    norm = plt.Normalize(vmin=min(charges), vmax=max(charges))
    cmap = plt.get_cmap('coolwarm')  # Choose a colormap, 'coolwarm' shows a gradient from blue (negative) to red (positive)
    colors = cmap(norm(charges))
    return colors[:, :3], norm, cmap  # Only need RGB values

# Function to write masses to the lammps data file
def write_masses(output_filename, elem_props):
    with open(output_filename, 'r') as ld_file:
        lines = ld_file.readlines()
    atoms_start_index = next(i for i, line in enumerate(lines) if line.startswith("Atoms "))# == "Atoms # full")
    masses_end_index = atoms_start_index
    # Extract the sigma and epsilon for each atom type
    mass_lines = ["Masses\n\n"]
    for (index, (atom_type, a_vals)) in enumerate(elem_props.items()):
        #print(atom_type + str(a_vals[2]) + "\n")
        mass_lines.append(f"{index+1:<4}   {a_vals[2]:.10f} # {a_vals[0]}\n")
    # Add in the blank line afterwards
    mass_lines.append("\n")
    # Insert the Pair Coeffs section between the 'Masses' and 'Atoms' sections
    updated_lines = (
        lines[:masses_end_index] +
        mass_lines +
        lines[atoms_start_index:]
    )

    # Write the updated file back
    with open(output_filename, 'w', newline='\n') as lammps_file:
        lammps_file.writelines(updated_lines)

# Function to read IFF frc file, search for pair coeffs A & B
def search_frc(elem_props):
    pair_coeffs = {}  # Dictionary to store coefficients for each type

    with open(frc_filename, 'r') as frc_file:
        lines = frc_file.readlines()

    # Find start and end of the relevant section (nonbond)
    start_line = next(i for i, line in enumerate(lines) if line.startswith("> E = Aij"))
    end_line = next(i for i, line in enumerate(lines) if line.startswith("#bond_increments"))

    # Extract the relevant lines
    relevant_lines = lines[start_line + 1:end_line]

    # Parse each forcefield type
    for ff_type in elem_props.values():
        for line in relevant_lines:
            if f" {ff_type[0]} " in line:  # Match forcefield type in the line
                parts = line.split()
                try:
                    pair_coeffs[ff_type[0]] = (float(parts[-2]), float(parts[-1]))  # Extract A and B
                except ValueError:
                    print(f"Error parsing line: {line.strip()}")

    return pair_coeffs

def write_pair_coeffs(output_filename, elem_props, pair_coeffs):
    """
    Inserts the Pair Coeffs section in the correct location in the LAMMPS data file.
    """
    # Read the existing data file
    with open(output_filename, 'r') as lammps_file:
        lines = lammps_file.readlines()

    # Locate the 'Masses' section and the start of the 'Atoms' section
    masses_end_index = next(i for i, line in enumerate(lines) if line.strip() == "Masses") + len(elem_props) + 1
    atoms_start_index = next(i for i, line in enumerate(lines) if line.startswith("Atoms "))

    # Extract the sigma and epsilon for each atom type
    pair_coeff_lines = ["\nPair Coeffs\n\n"]
    for (index, (atom_type, a_vals)) in enumerate(elem_props.items()):
        if a_vals[0] in pair_coeffs:
            A, B = pair_coeffs[a_vals[0]]
            sigma, epsilon = sigeps(A, B)
            pair_coeff_lines.append(f"{index+1:<4}   {epsilon:.10f}   {sigma:.10f} # {a_vals[0]}\n")
    # Add in the blank line afterwards
    pair_coeff_lines.append("\n")
    # Insert the Pair Coeffs section between the 'Masses' and 'Atoms' sections
    updated_lines = (
        lines[:masses_end_index+1] +
        pair_coeff_lines +
        lines[atoms_start_index:]
    )

    # Write the updated file back
    with open(output_filename, 'w', newline='\n') as lammps_file:
        lammps_file.writelines(updated_lines)


# Function to convert IFF A & B to sigma epsilon
def sigeps(A, B):
    #for key, (A, B) in pair_coeffs.items():
    sigma = (2 * A / B) ** (1/6)
    eps   = A / (sigma ** 12)
    return sigma, eps
    #pair_coeffs[key] = (sigma, eps)

# Build the base CsPbI3 perovskite
perov = Atoms(base_perov_elm, perov_locs)
perov.set_cell(cell_par[0] * np.identity(3))
perov.set_pbc((True, True, True))
supercell = build.make_supercell(perov, np.eye(3) * supercell_ABC)

# Corrected approach: Use len(supercell) to determine the total number of atoms
num_atom_anion = 3 * np.sum(supercell_ABC**2)
num_atom_cell = len(supercell)  # Automatically get the number of atoms in the supercell

# Assign charges and masses______________________________
charges = []
el_mass = []
forces  = []
for atom in supercell:
    # Round the atom positions
    #print(atom.position)
    #atom.position = [np.round(float(i), 6) for i in atom.position]
    if atom.symbol in elem_props:
        #atom.symbol = elem_props["I"][0]
        charges.append(elem_props[atom.symbol][1])
        el_mass.append(elem_props[atom.symbol][2])
        forces.append(elem_props[atom.symbol][0])
    else:
        charges.append(0)  # Default charge for any other atom
        el_mass.append(0)
        forces.append("")

# View before swapping - should be ideal CsPbI3
#view(supercell)
print("Initial net charge: ", np.sum(charges))
print("Initial mass of system: ", np.sum(el_mass))


# Swap I for Br ______________________________
# 1. Identify all iodine (I) atoms in the supercell
iodine_indices = [i for i, atom in enumerate(supercell) if atom.symbol == 'I']
num_iodine_atoms = len(iodine_indices)
num_Br_swaps = round(num_iodine_atoms * pBr)
num_Cl_swaps = round(num_iodine_atoms * pCl)

# 2. Create randomized subset of the iodine list
I_Br_swaps = rv_sample(iodine_indices, num_Br_swaps)
print("Total atoms in supercell:", num_atom_cell)
print("Total iodine atoms to swap:", num_Br_swaps)

# 3. Neighbor detection using ASE NeighborList
cutoff_distance = 3.5  # Example cutoff distance for finding nearest neighbors (tune if necessary)
nl = NeighborList([cutoff_distance / 2] * len(supercell), skin=0.0, self_interaction=False)
nl.update(supercell)

# Call swap_atom function
swap_atom("I", "Br", I_Br_swaps, elem_props, charges, supercell)


# Swap I for Cl ______________________________
iodine_indices = [i for i, atom in enumerate(supercell) if atom.symbol == 'I']
I_Cl_swaps = rv_sample(iodine_indices, num_Cl_swaps)

# Call swap_atom function
swap_atom("I", "Cl", I_Cl_swaps, elem_props, charges, supercell)


# Set the charges
print("Net charge: ", np.sum(charges))
supercell.set_initial_charges(charges)
supercell.set_masses([round(i, 4) for i in el_mass])
# View altered supercell
view(supercell)

# Save as data file
#masses = True
write_lammps_data(output_filename, supercell, atom_style='full') # , masses=True

# Find the nonbond pair coeffs for each atom declared up top
pair_coeffs = search_frc(elem_props)
# Convert pair coefficients A and B to sigma and epsilon
#sigeps(pair_coeffs)

# Write the masses: masses=True in write_lammps_data does not work
write_masses(output_filename, elem_props)
# Insert the Pair Coeffs section into the data file
write_pair_coeffs(output_filename, elem_props, pair_coeffs)


# Write accompanying nta file - defines atom types to assign - not currently needed
if write_nta:
    nta_filename = output_filename.replace('.data', '.nta')
    with open(nta_filename, 'w') as nta_file:
        nta_file.write("New type assignment file for {output_filename} > \n")
        nta_file.write("style id\n")
        for index, forcefield_type in enumerate(forces):
            nta_file.write(f"{index + 1:<6}{forcefield_type}\n")  # Align to column 7

# 20240911 - Doesn't really work right now / give insight forcefield params

if want_to_plot:
    # Color map ______________________________
    # Map charges to colors
    atom_colors, norm, cmap = map_charges_to_colors(charges)
    fig, ax = plt.subplots()
    heatmap = ax.imshow([atom_colors], aspect='auto')

    # If you want to save a snapshot of the structure with charge color visualization:
    cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='Charge')
    plt.title('Charge Distribution')
    plt.show()
