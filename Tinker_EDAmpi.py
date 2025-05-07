import subprocess
import os
import glob 
import time
import sys
import concurrent.futures
import threading
import statistics
import tempfile
import shutil


def get_header_and_repetitions(arc_file):
    """Reads the header line and counts its repetitions (frames)."""
    global line_to_search
    global repetitions
    
    line_to_search = None
    repetitions = None
    try:
        with open(arc_file, 'r', encoding='utf-8') as file:
            line_to_search = file.readline().strip()
    except FileNotFoundError:
        print(f"Error: The file {arc_file} was not found.")
        return None, None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None, None

    if line_to_search:
        print(f"Header line found: '{line_to_search}'")
        count = 0
        try:
            with open(arc_file, 'r', encoding='utf-8') as file:
                for line in file:
                    if line.strip() == line_to_search:
                        count += 1
            repetitions = count
        except FileNotFoundError:
            print(f"Error: The file {arc_file} was not found.")
            return None, None
        if count is not None:
             print(f"The header line repeats {count} times (frames).")
        return line_to_search, count
    else:
        print("Could not determine the header line.")
        return None, None

def read_tinker_xyz(arc_file):
    """
    Reads a TINKER ARC file and returns a list of atoms.
    """
    try:
        with open (arc_file, 'r') as file:
            try:
                n_atoms = int(file.readline().strip())
            except ValueError:
                print("Error: The first line of the file must be an integer representing the number of atoms.")
                return None
            
            n_atoms += 1
            atoms_list = []
            
            for _ in range(n_atoms):
                line = file.readline().strip()
                if not line:  # Skip empty lines
                    continue
                parts = line.split()
                if len(parts) < 7:  # Check if the line has enough data
                    continue
                atoms = {
                    'atom_num': int(parts[0]),
                    'atom_sym': parts[1],
                    'x': float(parts[2]),
                    'y': float(parts[3]),
                    'z': float(parts[4]),
                    'atom_type': int(parts[5]),
                    'connect': [int(num) for num in parts[6:]]                
                }
                atoms['connect'] = [num for num in atoms['connect'] if num != 0]
                atoms_list.append(atoms)
            return atoms_list
    except FileNotFoundError:
        print(f"Error: The file {arc_file} was not found.")
        return None

def atoms_to_residues(atoms_list):
    residues = []
    current_residue = []
    count_234 = 0
    for atom in atoms_list:
        current_residue.append(atom)
        
        # Assuming a new residue starts when the atom_sym is CA after N
        if len(current_residue) >= 3 and \
            current_residue[-3]['atom_sym'] == 'N' and \
            current_residue[-2]['atom_sym'] == 'CA' and \
            current_residue[-1]['atom_sym'] == 'C':
            
            residues.append(current_residue[:-3]) # Save the previous residue
            current_residue = current_residue[-3:] # Start a new residue 
            
        if atom['atom_type'] == 234:
            count_234 += 1
            if count_234 >= 2:
                break
        if atom['atom_type'] == 238:
            count_234 += 1
            if count_234 >= 2:
                break

    if current_residue:
        residues.append(current_residue)  # Append the last residue
    
    return residues

def atoms_to_nucleotide(atoms_list):
    nucleotides = []
    current_nucleotide = []
    
    for atom in atoms_list:
        if atom['atom_type'] == 349:
            break
        
        current_nucleotide.append(atom)

        # Assuming a new residue starts when the atom_sym is CA after N
        if len(current_nucleotide) >= 8 and \
            current_nucleotide[-8]['atom_sym'] == 'O' and \
            current_nucleotide[-7]['atom_sym'] == 'C' and \
            current_nucleotide[-6]['atom_sym'] == 'C' and \
            current_nucleotide[-5]['atom_sym'] == 'O' and \
            current_nucleotide[-4]['atom_sym'] == 'C' and \
            current_nucleotide[-3]['atom_sym'] == 'C' and \
            current_nucleotide[-2]['atom_sym'] == 'C' and \
            current_nucleotide[-1]['atom_sym'] == 'O' :
            
            nucleotides.append(current_nucleotide[:-8]) # Save the previous residue
            current_nucleotide = current_nucleotide[-8:] # Start a new residue 
        
    if current_nucleotide:
        nucleotides.append(current_nucleotide)  # Append the last residue
    
    return nucleotides

def atoms_to_water(atoms_list): 
    water = []
    water_molecule = []
    
    for atom in atoms_list:
        if atom['atom_type'] == 349:
            water_molecule = [atom]
            
            # Find connected atoms
            for connected_atom_num in atom['connect']:
                connected_atom = next((a for a in atoms_list if a['atom_num'] == connected_atom_num), None)
                if connected_atom:
                    water_molecule.append(connected_atom)
            
            # Check if we have a water molecule (e.g., at least one oxygen and some hydrogens)
            if len(water_molecule) > 1:
                water.append(water_molecule)
    
    return water


def analyze_residues(residues, prm_file):
    """
    Separates all of the possible pairs of the protein.
    """
    
    if residues:
        n = len(residues)
        
        # Function to analyze a single residue pair
        def analyze_pair(i, j):
            # Generate a unique filename
            pair_filename = f"pair_{i}_{j}_{threading.get_ident()}.arc"
            
            # Reset energy components for each pair
            all_energy_components = {}

            res1 = residues[i]
            res2 = residues[j]
            
            #Residue 1 atoms numbers
            res1_an1 = res1[0]['atom_num']
            res1_an2 = res1[-1]['atom_num']
            
            #Residue 2 atoms numbers
            res2_an1 = res2[0]['atom_num']
            res2_an2 = res2[-1]['atom_num']
            
            # Create a unique temp directory for this thread
            with tempfile.TemporaryDirectory(dir="./files") as tempdir:
                # Copy arc_file into tempdir with a unique name
                arc_basename = f"input_{threading.get_ident()}.arc"
                arc_copy = os.path.join(tempdir, arc_basename)
                # Use symlink if possible (faster than copy)
                try:
                    os.symlink(os.path.abspath(arc_file), arc_copy)
                except Exception:
                    shutil.copy2(arc_file, arc_copy)
                
                # Use only the basename in the input to archive
                command = "archive"
                if res1_an1 == 1 and res2_an1 != res1_an2 + 1:
                    input = f"{arc_basename}\n3\n-{res1_an2 + 1} -{res2_an1 - 1} -{res2_an2 + 1} -{line_to_search} 0\n1 {repetitions} 1\n\r\n"
                elif res1_an1 == 1 and res2_an1 == res1_an2 + 1:
                    input = f"{arc_basename}\n3\n -{res2_an2 + 1} -{line_to_search} 0\n1 {repetitions} 1\n\r\n"
                elif res1_an1 != 1 and res2_an1 != res1_an2 + 1:
                    input = f"{arc_basename}\n3\n-1 -{res1_an1 - 1} -{res1_an2 + 1} -{res2_an1 - 1} -{res2_an2 + 1} -{line_to_search} 1\n1 {repetitions} 1\n\r\n"
                elif res1_an1 != 1 and res2_an1 == res1_an2 + 1:
                    input = f"{arc_basename}\n3\n-1 -{res1_an1 - 1} -{res2_an2 + 1} -{line_to_search} 0\n1 {repetitions} 1\n\r\n"

                # Run archive in tempdir
                try:
                    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tempdir)
                    stdout, stderr = process.communicate(input=input.encode())
                    print(f"--- TINKER Output ---\n",stdout)
                    if stderr:
                        print(f"TINKER Error:\n", stderr)
                except FileNotFoundError:
                    print(f"Error: archive did not run correctly.")
                    return
                except Exception as e:
                    print(f"An error occurred: {e}")
                    return
                
                # Find the new output file in tempdir
                list_of_files = glob.glob(os.path.join(tempdir, '*'))
                newest_files = max(list_of_files, key=os.path.getctime)
                try:
                    os.rename(newest_files, pair_filename)
                except FileNotFoundError:
                    print(f"Error: The file '{newest_files}' was not found for renaming.")
                    return
                except OSError as e:
                    print(f"Error: Could not rename the file. {e}")
                    return
                
                #Connectivity check and delete 0 in the connect list
                # Read the file and process each atom's connect list
                try:
                    with open(pair_filename, 'r') as file:
                        lines = file.readlines()
                    
                    # Process each line (atom) in the file
                    for m in range(1, len(lines)):  # Skip header lines
                        parts = lines[m].split()
                        if len(parts) > 6:
                            connect_list = [int(num) for num in parts[6:]]
                            
                            # Remove 0 from the connect list
                            connect_list = [num for num in connect_list if num != 0]
                            
                            # Update the line with the modified connect list
                            lines[m] = ' '.join(parts[:6] + [str(num) for num in connect_list]) + '\n'
                    
                    # Write the modified content back to the file
                    with open(pair_filename, 'w') as file:
                        file.writelines(lines)
                        
                except FileNotFoundError:
                    print(f"Error: The file '{pair_filename}' was not found.")
                    return
                except Exception as e:
                    print(f"An error occurred while processing '{pair_filename}': {e}")
                    return
                    
                # Executes analyze.x process on the fly
                
                command = "analyze"
                a_res = pair_filename
                input = f"{a_res}\n{prm_file}\nE\n\r"
                try:
                    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    stdout, stderr = process.communicate(input=input)
                    
                    if stderr:
                        print(f"TINKER Error:\n", stderr)
                        print(f"DUDE WTH")
                        # Save the error message to a file
                        with open("error_log.txt", "a") as error_file:
                            error_file.write(f"TINKER Error for {pair_filename}:\n{stderr}\n")
                        sys.exit("Critical error in analyze, exiting.")
                    else:
                        # Parse the output
                        lines = stdout.splitlines()
                        energy_components = {}
                        for line in lines:
                            if "Intermolecular Energy :" in line:
                                if "Intermolecular Energy" not in energy_components:
                                    energy_components["Intermolecular Energy"] = []
                                energy_components["Intermolecular Energy"].append(float(line.split()[3]))
                            elif "Van der Waals" in line:
                                if "Van der Waals" not in energy_components:
                                    energy_components["Van der Waals"] = []
                                energy_components["Van der Waals"].append(float(line.split()[3]))
                            elif "Atomic Multipoles" in line:
                                if "Atomic Multipoles" not in energy_components:
                                    energy_components["Atomic Multipoles"] = []
                                energy_components["Atomic Multipoles"].append(float(line.split()[2]))
                            elif "Polarization " in line:
                                if "Polarization" not in energy_components:
                                    energy_components["Polarization"] = []
                                energy_components["Polarization"].append(float(line.split()[1]))
                        import statistics
                        # Prepare data for writing to the file
                        output_lines = []
                        
                        # Header
                        header_components = list(energy_components.keys())
                        header_line = " ".join([f"{c:<15}" for c in header_components])  # Adjust spacing as needed
                        output_lines.append(f"{header_line}\n")
                        
                        # Sub-header
                        subheader_line = " ".join(["AVG STD".center(15) for _ in header_components])
                        output_lines.append(f"{subheader_line}\n")
                        
                        # Data line
                        data_line = f"res {i} - res {j} "
                        for component in header_components:
                            values = energy_components[component]
                            n = len(values)
                            if n > 1:
                                avg = statistics.mean(values)
                                std_dev = statistics.stdev(values)
                            elif n == 1:
                                avg = values[0]
                                std_dev = 0.0  # Standard deviation is 0 if only one value
                            else:
                                avg = 0.0  # Handle the case where there are no values
                                std_dev = 0.0
                            data_line += f"{avg:<7.2f} {std_dev:<7.2f} "
                        output_lines.append(data_line + "\n")
                        
                        # Write to file
                        with open("energy_analysis.txt", "a") as outfile:  # "a" for append mode
                            outfile.writelines(output_lines)
                            # Delete the pair_{i}_{j}.arc file
                
                except FileNotFoundError:
                    print(f"Error: analyze did not run correctly.")
                except Exception as e:
                    print(f"An error occurred: {e}")
                
                try:
                    os.remove(pair_filename)
                    print(f"File {pair_filename} deleted successfully.")
                except FileNotFoundError:
                    print(f"Error: File {pair_filename} not found.")
                except Exception as e:
                    print(f"Error deleting file {pair_filename}: {e}")
        
        # Use a thread pool to parallelize the analysis
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            futures = []
            for i in range(1, n):
                for j in range(i + 1, n):
                    futures.append(executor.submit(analyze_pair, i, j))
            
            # Wait for all tasks to complete
            concurrent.futures.wait(futures)
    
    # Parallelize the nucleotides block as well
    if nucleotides:
        nu = len(nucleotides)
        def analyze_res_nuc_pair(i, j):
            # Reset energy components for each pair
            all_energy_components = {}

            res = residues[i]
            nuc = nucleotides[j]
            
            res1 = res[0]['atom_num']
            res2 = res[-1]['atom_num']
            
            nuc1 = nuc[0]['atom_num']
            nuc2 = nuc[-1]['atom_num']
            
            with tempfile.TemporaryDirectory(dir="./files") as tempdir:
                arc_basename = f"input_{threading.get_ident()}.arc"
                arc_copy = os.path.join(tempdir, arc_basename)
                # Use symlink if possible (faster than copy)
                try:
                    os.symlink(os.path.abspath(arc_file), arc_copy)
                except Exception:
                    shutil.copy2(arc_file, arc_copy)
                
                command = "archive"               
                if res1 == 1:
                    input = f"{arc_basename}\n3\n-{res2 + 1} -{nuc1- 1} -{nuc2 + 1} -{line_to_search} 1\n1 {repetitions} 1\n\r\n"
                elif res1 != 1:
                    input = f"{arc_basename}\n3\n-1 -{res1 - 1} -{res2 + 1} -{nuc1- 1} -{nuc2 + 1} -{line_to_search} 1\n1 {repetitions} 1\n\r\n"
                
                try:
                    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tempdir)
                    stdout, stderr = process.communicate(input=input.encode())
                    print(f"--- TINKER Output ---\n",stdout)
                    if stderr:
                        print(f"TINKER Error:\n", stderr)
                except FileNotFoundError:
                    print(f"Error: archive did not run correctly.")
                    return
                except Exception as e:
                    print(f"An error occurred: {e}")
                    return
                
                list_of_files = glob.glob(os.path.join(tempdir, '*'))
                newest_files = max(list_of_files, key=os.path.getctime)
                pair_filename = f"pair_res{i}_nuc{j}_{threading.get_ident()}.arc"
                try:
                    os.rename(newest_files, pair_filename)
                except FileNotFoundError:
                    print(f"Error: The file '{newest_files}' was not found for renaming.")
                    return
                except OSError as e:
                    print(f"Error: Could not rename the file. {e}")
                    return
                
                #Connectivity check and delete 0 in the connect list
                try:
                    with open(pair_filename, 'r') as file:
                        lines = file.readlines()
                    
                    # Process each line (atom) in the file
                    for m in range(1, len(lines)):  # Skip header lines
                        parts = lines[m].split()
                        if len(parts) > 6:
                            connect_list = [int(num) for num in parts[6:]]
                            connect_list = [num for num in connect_list if num != 0]
                            lines[m] = ' '.join(parts[:6] + [str(num) for num in connect_list]) + '\n'
                    with open(pair_filename, 'w') as file:
                        file.writelines(lines)
                except FileNotFoundError:
                    print(f"Error: The file '{pair_filename}' was not found.")
                    return
                except Exception as e:
                    print(f"An error occurred while processing '{pair_filename}': {e}")
                    return
                
                # Executes analyze.x process on the fly
                command = "analyze"
                a_res = pair_filename
                input = f"{a_res}\n{prm_file}\nE\n\r"
                try:
                    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    stdout, stderr = process.communicate(input=input)
                    
                    if stderr:
                        print(f"TINKER Error:\n", stderr)
                        print(f"DUDE WTH")
                        with open("error_log.txt", "a") as error_file:
                            error_file.write(f"TINKER Error for {pair_filename}:\n{stderr}\n")
                        sys.exit("Critical error in analyze, exiting.")
                    else:
                        lines = stdout.splitlines()
                        energy_components = {}
                        for line in lines:
                            if "Intermolecular Energy :" in line:
                                if "Intermolecular Energy" not in energy_components:
                                    energy_components["Intermolecular Energy"] = []
                                energy_components["Intermolecular Energy"].append(float(line.split()[3]))
                            elif "Van der Waals" in line:
                                if "Van der Waals" not in energy_components:
                                    energy_components["Van der Waals"] = []
                                energy_components["Van der Waals"].append(float(line.split()[3]))
                            elif "Atomic Multipoles" in line:
                                if "Atomic Multipoles" not in energy_components:
                                    energy_components["Atomic Multipoles"] = []
                                energy_components["Atomic Multipoles"].append(float(line.split()[2]))
                            elif "Polarization " in line:
                                if "Polarization" not in energy_components:
                                    energy_components["Polarization"] = []
                                energy_components["Polarization"].append(float(line.split()[1]))
                        # Prepare data for writing to the file
                        output_lines = []
                        header_components = list(energy_components.keys())
                        header_line = " ".join([f"{c:<15}" for c in header_components])
                        output_lines.append(f"{header_line}\n")
                        subheader_line = " ".join(["AVG STD".center(15) for _ in header_components])
                        output_lines.append(f"{subheader_line}\n")
                        data_line = f"res {i} - nuc {j} "
                        for component in header_components:
                            values = energy_components[component]
                            nvals = len(values)
                            if nvals > 1:
                                avg = statistics.mean(values)
                                std_dev = statistics.stdev(values)
                            elif nvals == 1:
                                avg = values[0]
                                std_dev = 0.0
                            else:
                                avg = 0.0
                                std_dev = 0.0
                            data_line += f"{avg:<7.2f} {std_dev:<7.2f} "
                        output_lines.append(data_line + "\n")
                        with open("energy_analysis.txt", "a") as outfile:
                            outfile.writelines(output_lines)
                except FileNotFoundError:
                    print(f"Error: analyze did not run correctly.")
                except Exception as e:
                    print(f"An error occurred: {e}")
                try:
                    os.remove(pair_filename)
                    print(f"File {pair_filename} deleted successfully.")
                except FileNotFoundError:
                    print(f"Error: File {pair_filename} not found.")
                except Exception as e:
                    print(f"Error deleting file {pair_filename}: {e}")

        # Use a thread pool to parallelize the residue-nucleotide analysis
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            futures = []
            for i in range(1, n):
                for j in range(1, nu):
                    futures.append(executor.submit(analyze_res_nuc_pair, i, j))
            concurrent.futures.wait(futures)
    else:
        print("Nothing to be processed was found.")
        return None

    
class Main:
    all_energy_components = {}

arc_file = "./files/test.arc"
prm_file = "./files/amoebabio18.prm"

get_header_and_repetitions(arc_file)
atoms_list = read_tinker_xyz(arc_file)
residues = atoms_to_residues(atoms_list)
nucleotides = atoms_to_nucleotide(atoms_list)
analyze_residues(residues, prm_file)



