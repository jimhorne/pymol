# -*- coding: utf-8 -*-
"""
ksCreated on Fri Apr  8 11:37:30 2022

@author: Jim

This script requires DSSP to be installed in windows or linux.

Can be downloaded from https://anaconda.org/speleo3/dssp (Windows)
or https://anaconda.org/salilab/dssp (Linux and OSX)

This has been updated from the original script made as part of the BamA-SurA project.
That script created a 2d array of filepaths and parent directories to link the name of the
AF2 prediction folder with the path to the file. This then feeds the prediction name and the path
to the PDBParser object to get the structure.

The colours assigned to each SS were determined by a lab-wide poll asking the question 
"What colour do you associate with each SS element" in the Radford Lab, Uni Leeds,
c. 2022. The most popular choices are used here.

This whole script takes 0.5 - 0.9s per entry.

2024-08-20 Working Points
--
"""

"""
DSSP codes
H = α-helix
B = residue in isolated β-bridge
E = extended strand, participates in β ladder
G = 3-helix (310 helix)
I = 5 helix (π-helix)
T = hydrogen bonded turn
S = bend
- = A blank in the DSSP secondary structure determination stands for loop or irregular.
Loops and irregular elements are often, very incorrectly, called "random coil" or "coil"

Tuple Index Value
0 = DSSP index
1 = Amino acid
2 = Secondary structure
3 = Relative ASA
4 = Phi
5 = Psi
6 = NH–>O_1_relidx
7 = NH–>O_1_energy
8 = O–>NH_1_relidx
9 = O–>NH_1_energy
10 = NH–>O_2_relidx
11 = NH–>O_2_energy
12 = O–>NH_2_relidx
13 = O–>NH_2_energy
"""

import time
start_time = time.time()

import os
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.DSSP import DSSP
import pymol
from pymol import cmd


def PDB_to_DSSP_and_PDB_obj(strucID, fname):
    """
    Requires PDBParser from Bio.PDB and DSSP from Bio.PDB.DSSP to be
    imported. Structure object follows a Structure/Model/Chain/Residue/Atom
    hierarchy architecture
    """
    p = PDBParser()
    structure = p.get_structure(strucID, fname)
    model = structure[0]
    return DSSP(model, fname, dssp='mkdssp'), model

def dssp_code_to_map(dsspObj):
    """
    Takes DSSP object secondary structure elements, maps key to dict values,
    and returns a list object. Requires PDBParser from Bio.PDB and DSSP from
    Bio.PDB.DSSP to be imported.
    """
    listOut = []
    for i, data in enumerate(dsspObj.keys()):
        a_key = list(dsspObj.keys())[i]
        listOut.append(structureMap.get(dsspObj[a_key][2]))
    return listOut
    #print('2D structure at each residue mapped to dictionary code and list created.')

def createFnameArray(globPattern, parentDir=False):
    """
    Requires glob and os to be imported. Takes a glob wildcard expression and returns a 2D numpy
    array containing the path+filename in one column and parent directory of file in
    another.
    """
    pathnames = []
    fnames = []
    subdirs = []
    if parentDir == True:
        for name in glob.glob(globPattern):
            pathnames.append(name)
            subdirs.append(os.path.basename(os.path.dirname(name)))
        return np.column_stack((pathnames, subdirs))
    else:
        for name in glob.glob(globPattern):
            pathnames.append(name)
            fnames.append(name)
        return np.column_stack((pathnames, fnames))

def plot_color_flag(structureMap):
    """
    Plots a list of DSSP assignments for protein sequence as a 1D color flag 
    according to whether residue a-helix, beta-sheet, or disordered. Plots data as a 
    color bar plot. Requires matplotlib.patches import Rectangle, matplotlib.collections 
    import PatchCollection, and matplotlib.colors import ListedColormap
    Args:
        structureMap (list): a 1D list of DSSP assignments were each index represents
        a single residue in the sequence
    """

    cmap = ListedColormap(['r', 'b', '0.8']) # set the colours of the 3 possible values in the list
    # Set width of figure proportional to length of sequence
    fig = plt.figure(figsize=((len(structureMap)/100), 1))
    ax = fig.add_axes([0, 0, 1, 1])
    # Turn of all axis labels
    ax.set_axis_off()
    # create a collection with rectangle for each residue
    col = PatchCollection([Rectangle((y, 0), 1, 1) for y in range(len(structureMap))])
    # set data, colormap and color limits
    col.set_array(np.array(structureMap))
    col.set_cmap(cmap)
    ax.add_collection(col)
    ax.set_ylim(0,1)
    ax.set_xlim(0, len(structureMap))
    # Show figure immediately and close after plotting
    plt.show()
    plt.close(fig)

def plot_histogram_sizes(input_data, key):
    """
    Takes in data and plots a histogram based on values for a specific
    key or column.
    Args:
        input_data (list or str): list of dictionaries where there is a key value in common you
        want to plot, or a string representing a path to a file to load into memory
        key (str): the key (dict) or column (file) you want to access to plot from
    """
    if type(input_data) == list:
        # Extract values from list using key - assumes list of dictionaries
        print('Plotting histogram from list of dictionaries.')
        values = [dic[key] for dic in input_data]
    
    else:
        print(f'Loading data from local file ({input_data}) and plotting.')
        df = load_spreadsheet(input_data)
        values = list(df[key])
    
    # Set the figure size
    plt.figure(figsize=(10, 6))
    
    # Plot a histogram
    plt.hist(values, bins=50, color='blue', edgecolor='black') # choose num bins or statistical rule forbinsize (auto, scott, sturges, fd, doane, rice, sqrt)

    # Set y-axis and x-axis range
    #plt.ylim(0, 50)
    #plt.xlim(0, 800)

    # Add labels and title
    plt.xlabel('Height in y (A)')
    plt.ylabel('Frequency')
    plt.title('Extent of AF2 predicted structure')

    # Show the plot
    #plt.show()
    
    # Save the plot to a file
    plt.savefig('../structural_bioinformatics/histogram_plot.png')

def load_spreadsheet(file_path, xlsx_file=False):
    '''Load data in a spreadsheet format (either TSV or XLSX/XLS) into memory as
    a pandas dataframe. Expects a header row. Checks extension to guess spreadsheet
    type. Assumes tab-seaprated rather than comma-separated file
    '''
    try:
        if 'xls' in file_path[-4:]:
            df = pd.read_excel(file_path, header=0)
        else:
            df = pd.read_csv(file_path, sep='\t', header=0)
        
        print(f'Check dataframe headers loaded in are correct: {df.columns.tolist()}')
        print(f'Check dataframe size is correct (excl. header row): {df.shape}')
        return df
    
    except Exception as e:
        print(f"Failed to process file {file_path}: {e}")

def trim_and_align_pymol(reference, sample, ctrim, ntrim, counter='1', align='align'):
    """
    Loads Pymol, loads PDB files, trims residues, aligns to reference PDB,
    saves trimmed and aligned file
    Args:
        reference (str): path to reference pdb file you want to align TO
        sample (str): path to sample file you want to align with ref
    """
    # Initialize PyMOL
    pymol.finish_launching(['pymol', '-cq']) # launch pymol in CLI mode quietly without GUI
        
    # Load the PDB files
    cmd.load(sample, 'sample')
    cmd.load(reference, 'reference')
    
    # Create model objects from loaded pymol objects
    sam_model = cmd.get_model('sample')
    ref_model = cmd.get_model('reference')
    
    # Calculate seq length directly from model
    len_sam = max([int(atom.resi) for atom in sam_model.atom])
    len_ref = max([int(atom.resi) for atom in ref_model.atom]) # need to convert to int (otherwise str)

    # Select and remove final Y residues from sample using DSSP calculated
    # length of disordered C-terminus
    Y = ctrim   # Number of residues to trim C-terminally
    cmd.remove(f'sample and resi {len_sam - Y + 1}-{len_sam}') # need to add 1 to select correct range  
    cmd.remove('reference and resi 465-493') # manually set for reference  

    # Select and remove first X residues
    X = ntrim	# Number of residues to trim N-terminally
    cmd.remove(f'sample and resi 1-{X}')
    cmd.remove('reference and resi 1-22') # manually set for reference

    # Align the sample structure to the reference and return the rmsd
    if align == 'align':
        rmsd = cmd.align('sample', 'reference')
    elif align == 'super':
        rmsd = cmd.super('sample', 'reference')
    else:
        print('Unrecognised alignment type. Exiting function early.')
        return

    # Save the aligned (and trimmed) structures to new PDB files
    parts = sample.split('-')
    sample_name = f'../structural_bioinformatics/TolC_pdbs/35pc_aligned_to_tolC_withSuper/sample_{parts[1]}_trim_aligned.pdb'
    cmd.save(sample_name, 'sample')
    if counter < 2:
        cmd.save('../structural_bioinformatics/TolC_pdbs/35pc_aligned_to_tolC_withSuper/reference_trim_aligned.pdb', 'reference')
    
    # Optional: Quit PyMOL or reinitialize
    cmd.reinitialize()
    #pymol.cmd.quit()
    
    return sample_name, rmsd

def calculate_xyz_extent_pdb(sample_name, verbose=False):
    """
    Loads aligned PDB file output from trim_and_align_pymol and calculates the extent
    in X, Y, Z on the coordinate system of the PDB file. This is why the structure must
    be aligned to a reference which has been manually aligned.
    """
    # Load the PDB file
    parser = PDBParser()
    structure = parser.get_structure('sample', sample_name)

    # Extract all atom coordinates
    all_atoms = [atom.get_coord() for atom in structure.get_atoms()]

    # Convert to a NumPy array for easier manipulation
    coords = np.array(all_atoms)

    # Calculate the 'height' along axes
    x_min = np.min(coords[:, 0])
    x_max = np.max(coords[:, 0])
    x_extent = x_max - x_min

    y_min = np.min(coords[:, 1])
    y_max = np.max(coords[:, 1])
    y_extent = y_max - y_min

    z_min = np.min(coords[:, 2])
    z_max = np.max(coords[:, 2])
    z_extent = z_max - z_min

    if verbose == True:
        print(f"The x-extent of the protein is: {x_extent:.3f} Å")
        print(f"The y-extent of the protein is: {y_extent:.3f} Å")
        print(f"The z-extent of the protein is: {z_extent:.3f} Å")
    
    return {'x': round(x_extent, 3), 'y': round(y_extent, 3), 'z': round(z_extent, 3)}

def export_dict_list_to_spreadsheet(dict_list):
    """
    Converts a list of dictionaries into a spreadsheet with the keys as headers
    """
    df = pd.DataFrame(dict_list)
    df.to_csv('../structural_bioinformatics/output.tsv', sep='\t', index=False)

def get_pdb_sequence(chain):
    ppb = PPBuilder()
    peptides = ppb.build_peptides(chain)
    
    # This should be redundant if your pdb chain is a single polypeptide stretch
    sequence = "".join([str(peptide.get_sequence()) for peptide in peptides])
                       
    return sequence

# ====================================
# Running analysis
# ====================================

# Create dictionary to map specific DSSP 2ndary structure elements to integer value
structureMap = {'H': 1, 'G': 1, 'I': 1, 'E': 2, 'B': 4, 'T': 4, 'S': 4, '-': 4}

# Pattern match for getting files, pathnames, and subdirectories - uses GLOB wildcard patterns not regex
# This matches everything in a specific depth from PWD containing certain string
#globPattern = r"AlphaFoldModels\*\*_1.pdb"
# This is used to match only files in PWD that contain AF2DB naming convention
globPattern = r'../structural_bioinformatics/TolC_pdbs/35pc_base/AF*.pdb'

# Create 2D numpy array of file path and file name
# If you want filepath and parent dir you parentDir=True
fnameArray = createFnameArray(globPattern)

# Specify PDB locations and names
reference_pdb = r'../structural_bioinformatics/reference_oriented_AF-P02930_tolCec_topump.pdb'

# Alignment method
alignment = 'super'

# Multi entry analysis
#confirmation = input(f'\n\n\nThe reference PDB is {reference_pdb} and the alignment method will be {alignment}.\n\n\n'
#                     'If this is correct please type "y".')
dict_list = []
counter = 0
try:
    for i, data in enumerate(fnameArray):
        counter += 1
        # Check if this entry has already been calculated before and skip if so
        basedir = os.path.dirname(data[0])
        alignTrimFile = f'sample_{data[1].split("-")[1]}_trim_aligned.pdb'
        final_file_check = os.path.join(basedir, alignTrimFile)
        if os.path.exists(final_file_check):
            continue
        # Main loop
        temp_dict = {}
        dsspObj, model = PDB_to_DSSP_and_PDB_obj(data[1], data[0])
        structure2Dmap = dssp_code_to_map(dsspObj)
        print(data)
        temp_dict['pdbPath'] = data[0]
        temp_dict['pdbName'] = data[1]
        temp_dict['len'] = len(get_pdb_sequence(model['A']))
        for i, value in enumerate(reversed(structure2Dmap), start=0): # iterates backwards (make it 1-indexed for 'true' index)
            if value != 4:
                cterm_disorder_len = len(structure2Dmap) - (len(structure2Dmap) - i) # this number is inclusive
                temp_dict['ctermdis'] = cterm_disorder_len
                print(f'Breaking... final {cterm_disorder_len} residues are disordered')
                break
            else:
                pass # can remove else statement but leaving here for readability
        for i, value in enumerate(structure2Dmap, start=0):
            if value != 4:
                nterm_disorder_len = i
                temp_dict['ntermdis'] = nterm_disorder_len
                print(f'Breaking... starting {nterm_disorder_len} residues are disordered')
                break
            else:
                pass
        #plot_color_flag(structure2Dmap)
        sam_name, rmsd = trim_and_align_pymol(reference_pdb, data[0], cterm_disorder_len, nterm_disorder_len, counter=counter, align=alignment)
        temp_dict['rmsd'] = float(rmsd[0])
        xyz_size = calculate_xyz_extent_pdb(sam_name)
        temp_dict.update(xyz_size)
        dict_list.append(temp_dict)
        print(f'Processed file {counter}.')
    export_dict_list_to_spreadsheet(dict_list)
    plot_histogram_sizes(dict_list, 'y')
    
except:
    print("Loop failed early. Saving data and plotting current results.")
    export_dict_list_to_spreadsheet(dict_list)
    plot_histogram_sizes(dict_list, 'y')



# =============================================================================
# # If have already performed analysis and wish to replot
# #precalc_output_file = r'C:\Users\james\OneDrive\Cambridge Postdoc Project Docs\Projects\OMF-PAP_length_distributions\interpro_35pc_resclust_tolc\output_interpro35pc.tsv'
# precalc_output_file = r'C:\Users\james\OneDrive\Cambridge Postdoc Project Docs\Projects\OMF-PAP_length_distributions\interpro_40pc_reclust_pap\output_interproPAP40pc_aligned_to_acrA_withSuper.tsv'
# plot_histogram_sizes(precalc_output_file, 'y')
# =============================================================================


# =============================================================================
# # Single entry analysis / debugging
# entry = 'AF-A0A0M2V3T8-F1-model_v4.pdb'
# dict_list = []
# counter = 0
# counter += 1
# temp_dict = {}
# dsspObj, model = PDB_to_DSSP_and_PDB_obj(data[1], data[0])
# structure2Dmap = dssp_code_to_map(dsspObj)
# print(entry)
# for i, value in enumerate(reversed(structure2Dmap), start=0): # iterates backwards (make it 1-indexed for 'true' index)
#     if value != 4:
#         cterm_disorder_len = len(structure2Dmap) - (len(structure2Dmap) - i) # this number is inclusive
#         print(f'Breaking... final {cterm_disorder_len} residues are disordered')
#         break
#     else:
#         pass # can remove else statement but leaving here for readability
# for i, value in enumerate(structure2Dmap, start=0):
#     if value != 4:
#         nterm_disorder_len = i
#         print(f'Breaking... starting {nterm_disorder_len} residues are disordered')
#         break
#     else:
#         pass
# #plot_color_flag(structure2Dmap)
# sam_name = trim_and_align_pymol(reference_pdb, entry, cterm_disorder_len)
# xyz_size = calculate_xyz_extent_pdb(sam_name)
# temp_dict.update(xyz_size)
# dict_list.append(temp_dict)
# 
# =============================================================================

end_time = time.time()
execution_time = end_time - start_time
print('\nExecution time:', execution_time, "seconds")
#print('\nTime per entry:', execution_time / counter, "seconds per entry")