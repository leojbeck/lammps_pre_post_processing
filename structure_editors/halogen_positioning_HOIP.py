# Halogen Placement Script - v2
# Original work -   William Weise   -- 2025 07 15
# Rework - Leo Beck - most recently -- 2025 08 20

# This script takes a master car / mdf output from Materials Studio,
# then creates all of the combinations of mono & di halogenated variants
# of a specific halogen (change hal to change the halogen type)
# This v2 also adjusts the mdf, not just the car

import os
import re
from copy import deepcopy

# Define substitution patterns (positions on ring)
placements = [[2], [3], [4], [2, 3], [2, 4], [2, 5], [2, 6], [3, 4], [3, 5]]

# Charge constants
e_c_hyd = 0.25		# cg2 with hydrogen attached
e_c_hal = 0.6		# cg2 with halogen attached
e_hyd = 0.050		# hydrogen on carbon ring (har)
e_hal = -0.300		# halogen on carbon ring (curr. only far)

# Forcefield types to check for
ffs = ['cg2', 'har']
hyd_ff = ffs[1]

# List of halogens to swap in
halogens = ['F', 'Cl', 'Br', 'I']

# car 2 mdf index shift
c2m = 17 


def car_format_fixed(p):
	return (
		f"{p[0]:<5}"						# Atom name
		f"{float(p[1]):15.9f}"				# X
		f"{float(p[2]):15.9f}"				# Y
		f"{float(p[3]):15.9f}"				# Z
		f"{p[4]:>5}"						# XXXX
		f"{' '}"
		f"{int(p[5]):<4}"					# Mol #
		f"{' ':>3}"							# 6 spaces
		f"{p[6]:<8}"						# Atom type (har, clar)
		f"{p[7]:<4}"						# Element type (H, Cl)
		f"{float(p[8].lstrip('.')):<0.3f}"	# Charge
	)

def mdf_format_fixed(m):
	string = ''
	for i in range(12, len(m)):
		string += f"{m[i]:<{len(m[i]) + 1}}"
		
	#print(string)
	return (
		f"{m[0]:<14}"						# SubUnit:element
		f"{m[1]:>7}"						# Element
		f"  "
		f"{m[2]:<4}"						# Forcefield type
		f"{m[3]:>5}"						# charge group (?)
		f"{int(m[4]):>6}"					# isotope (0)
		f"  "
		f"{str(m[5]):<2}"					# formal charge (str from 1- or 1+)
		f"{float(m[6]):>10.4f}"				# charge
		f"{int(m[7]):>2}"					# switching atom
		f"{int(m[8]):>2}"					# oop_flag
		f"{int(m[9]):>2}"					# chirality flag
		f"{float(m[10]):>7.4f}"				# occupancy
		f"{float(m[11]):>8.4f}"				# xray_temp_factor
		f"{string:>{len(string) + 1}}"		# connections (12 +)
	)

def build_atom_index(lines):
	"""Map atom names to their line index in the .car file"""
	atom_index = {}
	
	for i, line in enumerate(lines):
		parts = line.split()
		if len(parts) >= 2:
			# Grab atom name and atom index
			atom_name = parts[0]
			atom_index[atom_name] = i
	return atom_index

def replace_atoms_with_ring_substitution(car_lines, mdf_lines, placement, hal):
	# Read car / mdf lines
	car_lines = deepcopy(car_lines)
	mdf_lines = deepcopy(mdf_lines)
	atom_index = build_atom_index(car_lines)
	
	# Iterate through lines
	for ind, line in enumerate(car_lines[:-2]):
		# Ignore the header (first 4 lines)
		if ind < 4:
			continue
		
		# Split up car line into parts
		parts = line.split()
		fftype = parts[6]
		
		# Get relevant mdf line, split into parts
		mline = mdf_lines[ind+c2m]
		m_parts = mline.split()
		#print(len(m_parts))
		if len(m_parts) < 13:
			m_parts.insert(12, '')
		
		# Check if forcefield type matches any of the designated ff types
		if fftype.startswith(tuple(ffs)):
			# Check if the trailing char (ring placement) is in placements
			if int(fftype[-1]) in placement:
				
				# Check whether current line matches first ffs entry
				if fftype.startswith(hyd_ff):
					# Substituted with Fluorine
					parts[0] = f"{hal}{parts[0][1:]}"
					parts[6] = f"{hal}{parts[6][1:]}".lower()
					parts[7] = hal
					parts[8] = f"{e_hal:.3f}"
					
					# mdf line edits
					# TODO: C & H hardcoded in here
					m_parts[0] = re.sub(r":[CH](\d+)", fr":{hal}\1", m_parts[0])
				else:
					# Leave as carbon
					parts[8] = f"{e_c_hal:.3f}"
					tnums = int(re.search(r'\d+$', m_parts[12]).group())
					m_parts[12] = f"{hal}{tnums}"
			# Remove last character from ff type
			parts[6] = parts[6][:-1]
			m_parts[2]= m_parts[2][:-1]
		
		# Formatting for car & mdf
		line = car_format_fixed(parts)
		car_lines[ind] = line
		mline = mdf_format_fixed(m_parts)
		mdf_lines[ind+c2m] = mline
		
		# Fix the mdf file (add 17 to the line index)
		
	return car_lines, mdf_lines

def process_car_file(carif, mdfif):
	# Open and read lines from car / mdf files
	with open(carif, 'r') as f:
		original_car_lines = [line.rstrip('\n') for line in f]
	
	with open(mdfif, 'r') as f:
		original_mdf_lines = [line.rstrip('\n') for line in f]
	
	
	# Grab path and filename to use path
	car_base_filename, car_ext = os.path.splitext(os.path.basename(carif))
	mdf_base_filename, mdf_ext = os.path.splitext(os.path.basename(mdfif))
	dirpath = os.path.dirname(carif)
	
	for hal in halogens:
		haldir = dirpath + "/" + hal + "MBA"
		if not os.path.isdir(haldir):
			os.mkdir(haldir)
			print(" ------- Subdirectory created: ", haldir)
		else: 
			print(" ------- Subdirectory found. Using: ", haldir)

		for placement in placements:
			modified_car_lines, modified_mdf_lines = replace_atoms_with_ring_substitution(deepcopy(original_car_lines), deepcopy(original_mdf_lines), placement, hal)
			
			# car output things
			suffix = '-'.join(map(str, placement))
			output_car_filename = f"R-{suffix}-{hal}MBA{car_ext}" if suffix else f"R-{hal}MBA{car_ext}"
			output_car_path = os.path.join(haldir, output_car_filename)
			
			with open(output_car_path, 'w') as f:
				f.write('\n'.join(modified_car_lines) + '\n')
			print(f"Generated: {output_car_path}")
			
			# mdf output things
			output_mdf_filename = f"R-{suffix}-{hal}MBA{mdf_ext}" if suffix else f"R-{hal}MBA{mdf_ext}"
			output_mdf_path = os.path.join(haldir, output_mdf_filename)
			
			with open(output_mdf_path, 'w') as f:
				f.write('\n'.join(modified_mdf_lines) + '\n')
			print(f"Generated: {output_mdf_path}")

if __name__ == "__main__":
	car_in_file = "../../HOIS/CRYSTAL_STRUCTURES/Rxn_energy/R-Y-XMBA_structures/R-Y-XMBA_I_amorphous.car"
	#input_file = input("Enter full path to the .car structure file (including the file name itself): ").strip()
	if os.path.isfile(car_in_file) and car_in_file.endswith(".car"):
		mdf_in_file = os.path.splitext(car_in_file)[0]+'.mdf'
		if os.path.isfile(mdf_in_file):
			print("------ mdf found: ", mdf_in_file, " ------")
			process_car_file(car_in_file, mdf_in_file)
		else:
			print("car exists, but mdf does not. Please include mdf of same basename in the directory.")
	else:
		print("Invalid file. Please provide a valid .car file path.")
