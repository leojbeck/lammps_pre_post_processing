# Leo Beck
# Created: March 12th, 2025 
# Latest edits: August 15th, 2025

# ---- Desired descriptors ----
# Plots of energy, volume, etc. over time ✅
# Inorganics  
#   layer heights - from Pb to Pb ✅
#   octahedra distortions (delta d) ✅
#   PbI lengths ✅, IPbI angle ✅, PbIPb angle 🟨
#   symmetry - 21 screw axis & inversion centers 🟥
# Organics
#   organic orientation (wrt 0 0 1 ✅, inter ring  🟨)
#   inter-halogen distances and angles  🟨
# Interface
#   N distance ✅ / orientation 🟨 wrt Pb I pockets 
#       prioritize I over Pb

##%pip install ase
# Load relevant packages

import numpy as np
import pandas as pd
from ovito.io import import_file
from ovito.data import DataCollection
import MDAnalysis as mda
from openpyxl import Workbook
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
import re
import timeit
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase.data import atomic_masses, chemical_symbols
import gc

clock = timeit.default_timer
plt.style.use('seaborn-v0_8')

# Function to calculate plane normal using Singular Value Decomposition (SVD)
def get_plane_normal(points):
    centroid = np.mean(points, axis=0)
    shifted = points - np.mean(points, axis=0)
    _, _, v = np.linalg.svd(shifted)  # SVD
    normal = v[-1]  # Last column is the normal vector
    return normal / np.linalg.norm(normal)  # Normalize

def angle_between(v1, v2):
            cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

# ------------------ PARAMETERS ------------------
pb_i_cutoff = 3.6  # Pb-I cutoff (modify as needed)
pb_pb_cutoff = 2 * pb_i_cutoff
n_i_cutoff = 3.5 # n4 - Ipb1 cutoff
n_cutoff = 2.7 # n4 - cg2 cutoff
ring_cutoff = 2.9 # cg2 - cg2 cutoff
layer_threshold = 3 # Pb layer height minimum cutoff
xyz_axis = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # Reference for tilt angles
layer_dir = 1 # 0 - x, 1 - y, 2 - z
supercell = np.array([4, 2, 3])

system = 'PbI6'
sys_names = ['2', '3', '4', '23', '24', '25', '26', '34', '35']
halogen = ['F', 'Cl', 'Br', 'I']
halogen = ['F']

type_dict = {
            "pb2+": ["Pb", 207.200000, 1, 207.2],
            "Ipb1": ["I",  126.904400, 2, 126.90447],
            "n4"  : ["N",  14.006700,  3, 14.007],
            "hn"  : ["H",  1.007970,   4, 1.008],
            "har" : ["H",  1.007970,   5, 1.008],
            "far" : ["F",  18.998400,  6, 18.998403163],
            "cnp" : ["C",  12.011150,  7, 12.011],
            "h1h" : ["H",  1.007970,   8, 1.008],
            "cg2" : ["C",  10.011150,  9, 12.011],
            "c3"  : ["C",  12.011150, 10, 12.011],
            "cge" : ["H",  1.000000,  11, 1.0],
            "clar": ["Cl", 35.453000, 12, 35.45],
            "brar": ["Br", 79.904000, 13, 79.904],
            "iar" : ["I", 126.904470, 14, 126.90447]
    }
for hal in halogen:
    for it, placement in enumerate(sys_names):
        basename = f"R_{placement}_{hal}MBA_{system}"+"_4x2x3"
        fpath = f"../../HOIS/HOIS_OUTPUT/{system}/4x2x3/{hal}MBA/"+basename+"/"
        dump_file =   fpath + basename + ".dump"  # LAMMPS dump file
        data_file =  fpath + basename + ".data"  # LAMMPS data file
        output_file = fpath + basename + ".xlsx"  # Excel output file
        output_csv = fpath + basename + "avgs.csv" # csv averages output




        # ------------------ LOAD DATA ------------------
        pipeline = import_file(dump_file)

        # Prepare output storage
        data_list = []

        # Forcefield type dictionary: element & mass (currently) do nothing
        # FF type: element, mass, # in dump

        # ---------------- Update type_dict by reading the lammps data file -----------------
        """
        Parses the 'Masses' section of a LAMMPS data file to extract atom type mapping.
        Returns a dictionary where keys are forcefield names and values are LAMMPS atom type numbers.
        """
        atom_dict = {}
        in_masses_section = False
        atom_type = 0
        with open(data_file, "r") as f:
            for line in f:
                line = line.strip()

                # Detect start of the Masses section
                if line.startswith("Masses"):
                    in_masses_section = True
                    continue  # Skip the header line

                # Stop parsing when we reach another section
                if in_masses_section and line == "":
                    if atom_type:
                        break

                # Extract atom type ID and forcefield name (ignoring the actual mass)
                elif in_masses_section:
                    parts = line.split()
                    atom_type = int(parts[0])  # First integer is LAMMPS type ID

                    # Extract the forcefield type from the comment (if present)
                    if len(parts) > 2 and parts[2].startswith("#"):
                        ff_type = parts[3]  # Forcefield name is after "#"
                        atom_dict[ff_type] = atom_type


        # Load LAMMPS data file to get mass-based type dictionary

        # Print the dynamically extracted type mapping
        print("✅ Extracted Atom Type Mapping from LAMMPS Data File: " + basename + ".data")
        for ff_type, atom_type in atom_dict.items():
            type_dict[ff_type][2] = atom_type
            #print(f"{ff_type}: {atom_type}")

        #print(type_dict)


        # ------------------ PREDEFINE STUFF ------------------
        print("Number of frames: ", pipeline.source.num_frames)
        pb_ys = [_*[] for _ in range(pipeline.source.num_frames)]
        tot_org_projs = [[_*[] for _ in range(pipeline.source.num_frames)] for k in range(0, 3)]
        #tot_org_projs = np.reshape(tot_org_projs, (3, pipeline.source.num_frames))
        print("Number of total organic projections: ",  len(tot_org_projs))
        print("Length of first organic projection: ", len(tot_org_projs[0]))

        # ------------------ LOOP THROUGH TIMESTEPS ------------------
        for frame in range(pipeline.source.num_frames):
            data: DataCollection = pipeline.compute(frame)

            # Get atomic positions and types
            positions = data.particles.positions.array
            types = data.particles.particle_types.array
            # Extract timestep from the LAMMPS dump file
            timestep = data.attributes['Timestep']  # Get the timestep from OVITO

            # Identify atom indices and positions ______________________
            # Extract atom IDs from dump file
            atom_ids = data.particles.identifiers.array  # Atom IDs from dump file

            # Sort all atomic positions by atom ID
            sorted_indices = np.argsort(atom_ids)  # Get sorting indices
            sorted_positions = positions[sorted_indices]  # Reorder positions
            sorted_types = types[sorted_indices]  # Reorder types

            # Get atom indices
            pb_indices = np.where(types == type_dict["pb2+"][2])
            i_indices = np.where(types == type_dict["Ipb1"][2])
            cnp_indices = np.where(sorted_types == type_dict["cnp"][2])
            organic_indices = np.where(types == type_dict["cg2"][2])  # Modify if needed
            n_indices = np.where(types == type_dict["n4"][2])
            # Get atom positions
            pb_positions = positions[pb_indices]
            i_positions = positions[i_indices]
            cnp_positions = positions[cnp_indices]
            organic_positions = positions[organic_indices]
            n_positions = positions[n_indices]


            # Compute Cell Parameters - divided by supercell nx ny nz
            cell = data.cell
            a, b, c = cell[0], cell[1], cell[2]
            lx, ly, lz = cell[0,0] / supercell[0], cell[1,1] / supercell[1], cell[2,2] / supercell[2]
            alpha, beta, gamma = angle_between(b,c), angle_between(a,c), angle_between(a,b), 

            # Get Pb heights for the layer heights _______________________
            # 0 for x, 1 for y, 2 for z
            pb_z_coords = positions[pb_indices, layer_dir].flatten()
            if pb_z_coords.size == 0:
                print(f"⚠️ Frame {frame}: No Pb atoms found!")
                layer_height = np.nan
            else:
                pb_z_sorted = np.sort(pb_z_coords)

            #layer_threshold = 1.0
            layers = [[pb_z_sorted[0]]]

            for z in pb_z_sorted[1:]:
                if abs(z - layers[-1][-1]) > layer_threshold:
                    layers.append([z])
                else:
                    layers[-1].append(z)

            mean_z_per_layer = [np.mean(layer) for layer in layers]
            layer_diffs = np.diff(sorted(mean_z_per_layer))
            layer_height = np.mean(layer_diffs) if len(layer_diffs) > 0 else np.nan

            pb_ys[frame] = pb_z_coords  # Store raw z positions (for histograms or debugging)

            tilt_angles = []
            i_tree = KDTree(i_positions)  # Nearest-neighbor search for I atoms
            pb_tree = KDTree(pb_positions)

            # ________________________________ Compute N-I Distances  ________________________________
            n_distances = []
            for n in n_positions:
                nearest_i_indices = i_tree.query_ball_point(n, n_i_cutoff)
                n_distances.extend(np.linalg.norm(i_positions[nearest_i_indices] - n, axis=1))

            min_n_i_distance = np.min(n_distances) if n_distances else np.nan


            # ________________________________ Compute Octahedral Tilt Angles ________________________________
            # ________________________________    Calculate Pb I distances    ________________________________
            # TODO: Add in I-Pb-I angles (just flip method)

            nn_distances = []
            pb_distances = []
            delta_ds = []
            # Loop through each Pb atom
            for i, pb1 in enumerate(pb_positions):
                # Find neighboring Pb atoms within cutoff (6.5 Å)
                neighbor_indices = pb_tree.query_ball_point(pb1, pb_pb_cutoff)
                i_neighbors_pb1 = i_tree.query_ball_point(pb1, pb_i_cutoff)
                # Calculate Pb I distances
                pb_distances.extend(np.linalg.norm(i_positions[i_neighbors_pb1] - pb1, axis=1))
                # delta d: intra-octahedral distortion
                d0 = np.mean(pb_distances)
                d_d = (1/6)*np.sum((pb_distances - d0)**2) / (d0**2)
                #print(d_d)
                delta_ds.append(d_d)

                # Iterate through the neighboring Pb atoms
                for j in neighbor_indices:
                    if j <= i:  # Avoid duplicate pairs
                        continue
                    
                    pb2 = pb_positions[j]  # Neighbor Pb atom

                    # Find iodine neighbors for both Pb atoms (within 3.5 Å)
                    i_neighbors_pb1 = set(i_neighbors_pb1)
                    i_neighbors_pb2 = set(i_tree.query_ball_point(pb2, pb_i_cutoff))

                    # Get common I between the two Pb
                    candidate_indices = set(i_neighbors_pb1) & set(i_neighbors_pb2)
                    #if len(candidate_indices < 1):
                        #candidate_indices = set(i_neighbors_pb1) | set(i_neighbors_pb2)


                    best_angle = None
                    max_angle = -1000

                    # Select iodine with shortest Pb-I-Pb distance
                    for i_idx in candidate_indices:
                        i_pos = i_positions[i_idx]
                        v1 = pb1 - i_pos
                        v2 = pb2 - i_pos

                        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                        pb_i_pb_angle = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

                        # Check if this is a new best
                        if pb_i_pb_angle > max_angle:
                            max_angle, best_angle = pb_i_pb_angle, pb_i_pb_angle

                    # Append best angle to list of angles
                    if best_angle is not None:
                            tilt_angles.append(pb_i_pb_angle)
                    #print(pb1, " ... ", pb2, " ... ", cos_theta)
                    #print()

            avg_pb_i_distance = np.mean(pb_distances) if pb_distances else np.nan
            avg_delta_d = np.mean(delta_ds) #if delta_ds else np.nan

            avg_pbi_tilt = np.nanmean(tilt_angles) # Handle empty case
            #avg_pbi_tilt = np.max(tilt_angles) # Handle empty case
            #print(tilt_angles)


            # ________________________________ Compute Organic Molecule Tilt (Center of Mass Angles) ________________________________
            # TODO: Compute inter-ring angles

            # Nearest-neighbor search for organic rings
            organic_tilt = [np.nan] * (3*len(cnp_positions))  # List for storing tilt of organic molecules
            organic_tilt = np.reshape(organic_tilt, (3, len(cnp_positions)))
            otree = KDTree(organic_positions)  # KDTree for cg2 atoms

            # Iterate over cnp since there is only 1 per organic - will grab cg2 in position 1
            for i, n in enumerate(cnp_positions):
                #print(i, " __ ", n)
                # Step 1: Find the nearest cg2 atom to nitrogen
                nearest_ring = otree.query(n, k=1)  # Get closest cg2 atom (returns distance, index)
                nearest_cg2_index = nearest_ring[1]  # Get index of closest cg2

                # Step 2: Find all cg2 atoms in the same ring (within ring_cutoff)
                nearest_ring_atoms = otree.query_ball_point(organic_positions[nearest_cg2_index], ring_cutoff)

                if len(nearest_ring_atoms) >= 3:  # Ensure we have 6 atoms (benzyl ring)
                    ring_positions = organic_positions[nearest_ring_atoms]  # Get their positions

                    # Step 3: Compute plane normal using SVD
                    ring_normal = get_plane_normal(ring_positions)  # Normal vector of the ring

                    # Step 4: Compute tilt angle with X, Y, Z-axis
                    for _, v in enumerate(xyz_axis):
                        ring_angle = angle_between(ring_normal, v)
                        ring_angle = 180 - ring_angle if ring_angle >= 90 else ring_angle
                        organic_tilt[_,i] = ring_angle  # Store tilt angle for this organic molecule

                # ---- Inter ring angles ----
                # Find the closest cnp
                nn_dist = 999
                for j, n2 in enumerate(cnp_positions):
                    if i == j: continue
                    if nn_dist > np.linalg.norm(n-n2): nn_dist = np.linalg.norm(n-n2)
                #print(nn_dist)
                nn_distances.append(nn_dist)


            # Get the average projected tilts for this timestep

            avg_org_tilt_yz = np.nanmean(organic_tilt[0,:])
            avg_org_tilt_xz = np.nanmean(organic_tilt[1,:])
            avg_org_tilt_xy = np.nanmean(organic_tilt[2,:])
            dev_org_tilt_yz = np.nanstd( organic_tilt[0,:])
            dev_org_tilt_xz = np.nanstd( organic_tilt[1,:])
            dev_org_tilt_xy = np.nanstd( organic_tilt[2,:])
            #print(organic_tilt)
            for _ in range(0, 3):
                tot_org_projs[_][frame] = organic_tilt[_,:]
            avg_nn = np.nanmean(nn_distances)
            dev_nn = np.nanstd(nn_distances)
            min_nn = np.min(nn_distances)
            max_nn = np.max(nn_distances)
            #print(dev_nn)


            # Store results
            row = [timestep, lx, ly, lz, alpha, beta, gamma, layer_height, avg_pbi_tilt, avg_delta_d, avg_pb_i_distance, min_n_i_distance, avg_nn, min_nn, max_nn, dev_nn, avg_org_tilt_yz, avg_org_tilt_xz, avg_org_tilt_xy, dev_org_tilt_yz, dev_org_tilt_xz, dev_org_tilt_xy] #+ organic_tilt
            data_list.append(row)
        print("Donezo!")

        # ------------------ SAVE TO EXCEL ------------------
        # Make Columns
        lattice_params = ["a (Å)", "b (Å)", "c (Å)", "alpha (°)", "beta (°)", "gamma (°)"]
        inorg_cols = ["Layer Height (Å)", "Octahedral Tilt (°)", "delta d", "Pb-I Distance (Å)", "N-I Distance (Å)", "Avg min N-N Distance (Å)", "Min N-N Dist (Å)", "Max N-N Dist (Å)","N-N StdDev"]
        organic_cols = [f"Avg Organic Tilt in {i} (°)" for i in ['x', 'y', 'z']]
        dev_cols = [f"StdDev Organic Tilt in {i} (°)" for i in ['x', 'y', 'z']]
        columns = ["Timestep"] + lattice_params + inorg_cols + organic_cols + dev_cols

        # Make and save df
        df = pd.DataFrame(data_list, columns=columns)
        df.to_excel(output_file, index=False)
        print(f"✅ Analysis complete! Results saved to {output_file}")

        # ------------------ List Time / Block Averages ------------------
        # Create an array / dataframe that has the back half averages
        bha_avgs = []
        # Append mean, median, and std dev to bha_avgs
        bha_avgs.append(df.iloc[round(df.shape[0]/2):-1, :].mean(axis=0))
        bha_avgs.append(df.iloc[round(df.shape[0]/2):-1, :].median(axis=0))
        bha_avgs.append(df.iloc[round(df.shape[0]/2):-1, :].std(axis=0))

        #Create dataframe, add in column for the basename (for combining csv later)
        df_bha = pd.DataFrame(bha_avgs, columns=columns)
        df_bha.insert(0, 'Basename', np.array([basename,basename,basename]).T)
        df_bha.insert(1, 'Halogen', np.array([hal,hal,hal]).T)
        df_bha.insert(2, 'Placement', np.array([placement,placement,placement]).T)
        df_bha.insert(3, 'System', np.array([system,system,system]).T)
        df_bha.insert(4, 'Stat', ['Mean', 'Median', 'StdDev'])


        # Save to csvs
        df_bha.to_csv(output_csv, index=False)
        print(f'Saved to {output_csv}!')

        # ------------------ MAKE PLOTS ------------------
        # Load the DataFrame from the Excel file
        df = pd.read_excel(output_file)

        # Define the time column
        time = df["Timestep"] * 1.0  # Convert to fs if needed

        # Loop through each column (except Timestep) and plot
        for column in df.columns[1:]:  # Skip "Timestep"
            plt.figure(figsize=(8, 5))
            plt.plot(time, df[column], marker='o', linestyle='-', label=column)

            # Formatting
            plt.title(f"{basename}  {column} vs Time", fontsize=14)
            plt.xlabel("Time (fs)", fontsize=12)
            plt.ylabel(column, fontsize=12)
            plt.legend()
            plt.grid(True)

            # Save the plot as an image
            plt.savefig(f"{fpath}{basename}_{column.replace(' ', '_')}_vs_Time.png", dpi=300, bbox_inches="tight", facecolor="w")
            plt.close()  # Close to prevent overlap issues

        print("✅ Plots saved successfully!")

        # ------------------ MAKE Organic PLOTS ------------------

        plt.figure(figsize=(8, 5))
        for column in organic_cols:  # Skip "Timestep"
            plt.plot(time, df[column], label=column)

        # Formatting
        plt.title(f"{basename} Organic Angle vs Time", fontsize=14)
        plt.xlabel("Time (fs)", fontsize=12)
        plt.ylabel("Organic Tilt (°)", fontsize=12)
        plt.legend(["yz plane", "xz plane", "xy plane"])
        plt.grid(True)
        plt.ylim([0, 90])

        # Save the plot as an image
        plt.savefig(f"{fpath}{basename}_Organic_tilt_vs_Time.png", dpi=300, bbox_inches="tight", facecolor="w")
        plt.close()
        #plt.close()  # Close to prevent overlap issues

        print("✅ Organic Tilt Plots saved successfully!")
        #plt.show()

        # ------------------ MAKE BINS ------------------
        # ------------------ LEAD POSITIONS ------------------
        plt.figure()
        plt.hist(pb_ys[-1], bins=61)

        xyz_labels = ['X', 'Y', 'Z']
        # Adding labels and title
        plt.xlabel(f'Pb {xyz_labels[layer_dir]} Height (A)')
        plt.ylabel('Frequency')
        plt.title(f'Histogram of Pb {xyz_labels[layer_dir]} in {basename}')

        plt.savefig(f"{fpath}{basename}_Pb_{xyz_labels[layer_dir]}_Height_at_End.png", dpi=300, bbox_inches="tight", facecolor="w")
        plt.close()
        # Display the plot
        #plt.show()

        # ------------------ MAKE BINS ------------------
        # ------------------ ring angles ------------------
        for i, k in enumerate(['x', 'y', 'z']):
            avg_tots = np.mean(tot_org_projs[i][-10:-1], axis=0)
            plt.figure()
            plt.hist(avg_tots, bins=41)

            # Adding labels and title
            plt.xlabel(f'Ring Angle against {k} plane')
            plt.ylabel('Frequency')
            plt.title(f'Histogram of {basename} Ring Projection onto the {k} Plane')

            plt.savefig(f"{fpath}{basename}_RingAngle_at_End_{k}.png", dpi=300, bbox_inches="tight", facecolor="w")
            plt.close()
        # Display the plot
        #plt.show()
        print(f'Done analyzing {basename} at time {clock()}!')

        # Garbage collection / memory clearing
        del [df, df_bha, pipeline, data]
        gc.collect()
