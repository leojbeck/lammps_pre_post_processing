# MBA backbone alteration - chirality swapper
# Originally adapted from halogen_positioning_script (2025 08 21)
# Leo Beck - most recently -- 2025 09 13

# This script takes a car / mdf output from Materials Studio,
# then edits the 'backbones' to change chirality, and (remove the methyl group)

# CURRENT ASSUMPTIONS
# c3 shows up in mdf as multiple of 8 
# h1h on cnp has been changed to h1hx


import os
import re
from copy import deepcopy


# Dictionary of backbone charges [0] and elements [1]
elem_dict = {
	'n4': [-0.4, 'N'],
	'c3': [-0.15, 'C'], #[-0.26, 0.56],
	'hn': [0.4, 'H'], #[-0.28, 0.58],
	'h1h': [0.05, 'H'], #[-0.24, 0.54]
	'cnp': [0.15, 'C'],
	}

# Forcefield types to check for
ffs = ['n4', 'c3', 'hn', 'h1h', 'cnp', 'h1hx']

# car 2 mdf index shift based on header lengths
c2m = 16
# Store c3 ids from mdf during chir_flag 1 to reference in chir_flag 2
c3_ids = {}

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
	# For loop to handle dynamic length from # bonds
	string = ''
	for i in range(12, len(m)):
		string += f"{m[i]:<{len(m[i]) + 1}}"
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

def car_swap(parts, n, new_id):
	#n = {0: 1, 1: 0, 2: 3, 3: 2}
	parts[0] = f"{elem_dict[ffs[n]][1]}{new_id}"
	parts[6] = f"{ffs[n]}".lower()
	parts[7] = elem_dict[ffs[n]][1]
	parts[8] = f"{elem_dict[ffs[n]][0]:.3f}"
	return parts

def mdf_swap(mp, n, new_id):
	mp[0] = re.sub(r":[A-Z]+(\d+)", fr":{elem_dict[ffs[n]][1]}{new_id}", mp[0])
	mp[1] = f"{elem_dict[ffs[n]][1]}"
	mp[2] = f"{ffs[n]}".lower()
	mp[6] = elem_dict[ffs[n]][0]
	return mp

def swap_backbone_atoms(car_lines, mdf_lines, chir_flag):
	# Read car / mdf lines
	car_lines = deepcopy(car_lines)
	mdf_lines = deepcopy(mdf_lines)
	atom_index = build_atom_index(car_lines)
	
	delind = 888
	start_ind = 0
	# Iterate through lines
	for ind, line in enumerate(car_lines[:-2]):
		# Ignore the header (first 4 lines)
		if ind <= 4:
			start_ind += 1
			continue
		
		# Split up car line into parts
		parts = line.split()
		fftype = parts[6]
		
		# Get relevant mdf line, split into parts
		mline = mdf_lines[ind+c2m]
		m_parts = mline.split()
		
		# Handle case of no bonds (likely from Ipb1)
		if len(m_parts) < 13:
			m_parts.insert(12, '')
		
		
		# Check if forcefield type matches any of the designated ff types
		if fftype.startswith(tuple(ffs)):
			# Check which entry line matches 
			# (nitrogen)
			match fftype:
				case 'n4': # Nitrogen n4
					if chir_flag == 1:
						# mdf line edits
						# mdf: map nitrogen to carbon (id * 8)
						atom_id = int(re.search(r"\d+$", m_parts[0]).group())
						new_id = atom_id * 8
						# Swap n4 to c3 
						car_swap(parts, 1, new_id)
						mdf_swap(m_parts, 1, new_id)
				case 'c3' : # Carbon c3
					if chir_flag == 1:
						# mdf: map carbon to nitrogen (id / 8)
						atom_id = int(re.search(r"\d+$", m_parts[0]).group())
						
						new_id = atom_id // 8   # integer division
						
						# Swap c3 to n4
						car_swap(parts, 0, new_id)
						mdf_swap(m_parts, 0, new_id)
						
						# Add to dict of c3 ids (for chir_flag = 2)
						c3_ids.update({atom_id: delind})
						delind += 1
					elif chir_flag == 2:
						new_id = int(re.search(r"\d+$", m_parts[0]).group())
						new_id = c3_ids[new_id]
						car_swap(parts, 3, new_id)
						mdf_swap(m_parts, 3, new_id)
						# Delete bonded hydrogens in mdf
						m_parts[12:] = [col for col in m_parts[12:] if not re.match(r"^H\d+$", col)]
						
				case 'hn': # Hydrogen (hn)
					# Swap hn to h1h
					if chir_flag == 1:
						new_id = parts[0][1:]
						car_swap(parts, 3, new_id)
						# Change bonded n to bonded h1h
						new_id = int(re.search(r"\d+$", m_parts[12]).group()) * 8
						m_parts[12] = f"C{new_id}"
						
				case 'h1h': # Hydrogen (h1h on c3)
					if chir_flag == 1:
						# Swap h1h to hn
						new_id = parts[0][1:]
						car_swap(parts, 2, new_id)
						new_id = int(re.search(r"\d+$", m_parts[12]).group()) // 8
						m_parts[12] = f"N{new_id}"
					elif chir_flag == 2:
						parts = None
						
				case 'cnp':
					if chir_flag == 2:
						# replace bonded C with H in MDF
						for i, bond in enumerate(m_parts[12:], start=12):
							# Look for :C<number>
							m = re.match(r"^C(\d+)$", bond) or re.match(r"^:C(\d+)$", bond)
							if m:
								c_id = int(m.group(1))
								if c_id in c3_ids:
									m_parts[i] = f"H{c3_ids[c_id]}"
						
				case _: # Hydrogen (h1h on cnp)
					# Remove last character from ff type
					parts[6] = parts[6][:-1]
					m_parts[2]= m_parts[2][:-1]
		
		# Formatting for car & mdf
		if parts is not None:
			line = car_format_fixed(parts)
			mline = mdf_format_fixed(m_parts)
		else:
			line = ""
			mline = ""
		
		car_lines[ind] = line
		mdf_lines[ind+c2m] = mline
		
	print(c3_ids)
	# Delete empty lines
	car_lines[start_ind:] = [l for l in car_lines[start_ind:] if l.strip()]
	mdf_lines[start_ind+c2m:] = [l for l in mdf_lines[start_ind+c2m:] if l.strip()]
	return car_lines, mdf_lines

def process_car_file(carif, mdfif, chir_flag):

	carfp = os.path.dirname(carif)
	mdffp = os.path.dirname(mdfif)
	# Grab path and filename to use path
	car_base_filename, car_ext = os.path.splitext(os.path.basename(carif))
	mdf_base_filename, mdf_ext = os.path.splitext(os.path.basename(mdfif))
	
	# Open and read lines from car / mdf files
	with open(carif, 'r') as f:
		original_car_lines = [line.rstrip('\n') for line in f]
	
	with open(mdfif, 'r') as f:
		original_mdf_lines = [line.rstrip('\n') for line in f]
	
	modified_car_lines, modified_mdf_lines = swap_backbone_atoms(deepcopy(original_car_lines), deepcopy(original_mdf_lines), chir_flag)
	
	# Change chirality char based on operation
	chirality = car_base_filename[0:2]
	if chir_flag == 1:
		if car_base_filename.startswith('R'):
			chirality = 'S' + str(chirality[1])
	elif chir_flag == 2:
		chirality = ''
	else:
		chirality = 'c' + chirality
	# car mdf output filenames
	output_car_filename = f"{chirality}{car_base_filename[2:]}{car_ext}"
	output_mdf_filename = f"{chirality}{mdf_base_filename[2:]}{mdf_ext}"
	output_car_path = os.path.join(carfp, output_car_filename)
	output_mdf_path = os.path.join(mdffp, output_mdf_filename)
	
	# Write new files
	with open(output_car_path, 'w') as f:
		f.write('\n'.join(modified_car_lines) + '\n')
	print(f"Generated: {output_car_path}")
	
	with open(output_mdf_path, 'w') as f:
		f.write('\n'.join(modified_mdf_lines) + '\n')
	print(f"Generated: {output_mdf_path}")

if __name__ == "__main__":
	car_in_file = "../../../HOIS/CRYSTAL_STRUCTURES/PbI4/master_files/R-X-YMBA_PbI4.car"
	#input_file = input("Enter full path to the .car structure file (including the file name itself): ").strip()
	if os.path.isfile(car_in_file) and car_in_file.endswith(".car"):
		mdf_in_file = os.path.splitext(car_in_file)[0]+'.mdf'
		if os.path.isfile(mdf_in_file):
			print("------ mdf found: ", mdf_in_file, " ------")
			for chir_flag in range(0,3):
				process_car_file(car_in_file, mdf_in_file, chir_flag)
		else:
			print("car exists, but mdf does not. Please include mdf of same basename in the directory.")
	else:
		print("Invalid file. Please provide a valid .car file path.")
