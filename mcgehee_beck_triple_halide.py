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

# Define size, percentages, et c.
pBr = 0.30
pCl = 0.06
cell_par = [6.289, 6.289, 6.289, 90, 90, 90]
base_perov_elm = 'CsPbI3'
supercell_ABC  = np.array([2, 2, 2])
frc_filename = "cvff_iff_v1_5_G_poly_solv_ions_per_v55.frc"
# Define dictionary of element properties: element, charge, mass
# Put in alphabetical order please! :)
elem_props = {
  "Br": ["Br", -0.88, 79.904],
  "Cl": ["Cl", -0.9, 35.453],
  "cs+": ["Cs", 1.0, 132.91],
  "Ipb1": ["I", -0.84, 126.9044],
  "pb2+": ["Pb", 1.52, 207.2]
}
# Define dictionary of element charges
elem_charges = {
  "Cs": 1.0,
  "Pb": 1.52,
  "I": -0.84,
  "Br": -0.88,
  "Cl": -0.9
}

# Define dictionary of element masses
elem_masses = {
  "Cs": 132.91,
  "Pb": 207.2,
  "I": 126.9044,
  "Br": 79.904,
  "Cl": 35.453
}

# Define dictionary of forcefield types
ff_types = {
  "Cs": "cs+",
  "Pb": "pb2+",
  "I": "Ipb1",
  "Br": "Br",
  "Cl": "Cl"
}

# True / False for plotting charge density
want_to_plot = False

# Other stuff that doesn't need to change as often
a = cell_par[0]
b = cell_par[1]
c = cell_par[2]
perov_locs = [(a/2, a/2, a/2), (0, 0, 0), (a/2, 0, 0), (0, a/2, 0), (0, 0, a/2)] 

output_filename = "CsPbI"+str(int(100*round((1-pBr-pCl),2)))+"Br"+str(int(100*pBr))+"Cl"+str(int(100*pCl))+".data"

# Function to make cleaner code later
def swap_atom(pre, post, list, elch, charges, supercell):
    for count, rv in enumerate(list):
        #print(rv)
        # Swap and change charge
        supercell[rv].symbol = ff_types[post]
        charges[rv] = elem_charges[post]
        el_mass[rv] = round(elem_masses[post], 4)
        forces[rv]  = ff_types[post]
        # Find nearest neighbors, update charge
        # Find the nearest neighbors of the swapped atom
        indices, offsets = nl.get_neighbors(rv)
        
        # Distribute the charge difference between I and Br (-0.88 - (-0.84) = -0.04)
        # Split the difference among the neighbors
        charge_difference = round(((elch[pre] - elch[post]) / len(indices)), 2)  

        # Update nearest neighbor charge
        for neighbor_index in indices:
            charges[neighbor_index] += charge_difference

# Function to map charges to colors using a colormap
def map_charges_to_colors(charges):
    norm = plt.Normalize(vmin=min(charges), vmax=max(charges))
    cmap = plt.get_cmap('coolwarm')  # Choose a colormap, 'coolwarm' shows a gradient from blue (negative) to red (positive)
    colors = cmap(norm(charges))
    return colors[:, :3], norm, cmap  # Only need RGB values

# Function to read IFF frc file, search for pair coeffs A & B
def search_frc(ff_types):
    pair_coeffs = {}  # Dictionary to store coefficients for each type

    with open(frc_filename, 'r') as frc_file:
        lines = frc_file.readlines()

    # Find start and end of the relevant section
    start_line = next(i for i, line in enumerate(lines) if line.startswith("> E = Aij"))
    end_line = next(i for i, line in enumerate(lines) if line.startswith("#bond_increments"))

    # Extract the relevant lines
    relevant_lines = lines[start_line + 1:end_line]

    # Parse each forcefield type
    for ff_type in ff_types.values():
        for line in relevant_lines:
            if f" {ff_type} " in line:  # Match forcefield type in the line
                parts = line.split()
                try:
                    pair_coeffs[ff_type] = (float(parts[-2]), float(parts[-1]))  # Extract A and B
                except ValueError:
                    print(f"Error parsing line: {line.strip()}")

    return pair_coeffs

def write_pair_coeffs(output_filename, ff_types, pair_coeffs):
    """
    Inserts the Pair Coeffs section in the correct location in the LAMMPS data file.
    """
    # Read the existing data file
    with open(output_filename, 'r') as lammps_file:
        lines = lammps_file.readlines()

    # Locate the 'Masses' section and the start of the 'Atoms' section
    masses_end_index = next(i for i, line in enumerate(lines) if line.strip() == "Masses") + len(elem_props) + 1
    atoms_start_index = next(i for i, line in enumerate(lines) if line.strip() == "Atoms # full")

    # Extract the sigma and epsilon for each atom type
    pair_coeff_lines = ["\nPair Coeffs\n\n"]
    for (index, (atom_type, ff_type)) in enumerate(ff_types.items()):
        if ff_type in pair_coeffs:
            A, B = pair_coeffs[ff_type]
            sigma, epsilon = sigeps(A, B)
            pair_coeff_lines.append(f"{index+1:<4}   {epsilon:.10f}   {sigma:.10f} # {ff_type}\n")
    pair_coeff_lines.append("\n")
    # Insert the Pair Coeffs section between the 'Masses' and 'Atoms' sections
    updated_lines = (
        lines[:masses_end_index] +
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
    if atom.symbol == 'I':
        charges.append(elem_charges["I"])
        el_mass.append(elem_masses["I"])
        forces.append(ff_types["I"])
    elif atom.symbol == 'Br':
        charges.append(elem_charges["Br"])
        el_mass.append(elem_masses["Br"])
        forces.append(ff_types["Br"])
    elif atom.symbol == 'Cl':
        charges.append(elem_charges["Cl"])
        el_mass.append(elem_masses["Cl"])
        forces.append(ff_types["Cl"])
    # Assign charges for other atoms (Cs, Pb) as needed
    elif atom.symbol == 'Cs':
        charges.append(elem_charges["Cs"])  # Example charge for Cs
        el_mass.append(elem_masses["Cs"])
        forces.append(ff_types["Cs"])
    elif atom.symbol == 'Pb':
        charges.append(elem_charges["Pb"])  # Example charge for Pb
        el_mass.append(elem_masses["Pb"])
        forces.append(ff_types["Pb"])
    else:
        charges.append(0)  # Default charge for any other atom
        el_mass.append(0)
        forces.append("")

# View before swapping - should be ideal CsPbI3
view(supercell)
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
swap_atom("I", "Br", I_Br_swaps, elem_charges, charges, supercell)


# Swap I for Cl ______________________________
iodine_indices = [i for i, atom in enumerate(supercell) if atom.symbol == 'I']
I_Cl_swaps = rv_sample(iodine_indices, num_Cl_swaps)

# Call swap_atom function
swap_atom("I", "Cl", I_Cl_swaps, elem_charges, charges, supercell)


# Set the charges
print("Net charge: ", np.sum(charges))
supercell.set_initial_charges(charges)
supercell.set_masses([round(i, 4) for i in el_mass])
# View altered supercell
view(supercell)

# Save as data file
#masses = True
write_lammps_data(output_filename, supercell, atom_style='full', masses=True) # , masses=True

# Find the nonbond pair coeffs for each atom declared up top
pair_coeffs = search_frc(ff_types)
# Convert pair coefficients A and B to sigma and epsilon
#sigeps(pair_coeffs)

# Insert the Pair Coeffs section into the data file
write_pair_coeffs(output_filename, ff_types, pair_coeffs)


# Write accompanying nta file - defines atom types to assign
nta_filename = output_filename.replace('.data', '.nta')
with open(nta_filename, 'w') as nta_file:
    nta_file.write("New type assignment file for {output_filename} > \n")
    nta_file.write("style id\n")
    for index, forcefield_type in enumerate(forces):
        nta_file.write(f"{index + 1:<6}{forcefield_type}\n")  # Align to column 7

# 20240911 - Doesn't really work right now / give insign forcefield params


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