# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 14:38:55 2024

@author: james
"""

import csv
import os
import pymol

# Function to open a subset of models
def open_subset_models(tsv_file):  
    models = []
    
    # Read the TSV file and process it
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            models.append((float(row['y']), row['pdbName']))
    
    # Sort list based on chosen metric
    models.sort(key=lambda x: x[0])

    # Open models
    for _, pdb_name in models[8515:8619:1]:
        mirror_dir = r'C:\Users\james\OneDriveOverflow\Scripting\StructuralBioinformatics' # comment out this line if files in pwd
        alignTrimFile = f'TolC_pdbs/35pc_aligned_to_tolC_withAlign/sample_{pdb_name.split("-")[1]}_trim_aligned.pdb'
        if mirror_dir:
            alignTrimFile = os.path.join(mirror_dir, alignTrimFile)
        pymol.cmd.load(alignTrimFile)

# Example usage: Replace 'your_file.tsv' with the path to your TSV file
open_subset_models(r'C:\Users\james\OneDrive\Cambridge Postdoc Project Docs\Projects\OMF-PAP_length_distributions\interpro_35pc_resclust_tolc\output_interpro35pc.tsv')
#open_subset_models(r'C:\Users\james\OneDrive\Cambridge Postdoc Project Docs\Projects\OMF-PAP_length_distributions\interpro_40pc_reclust_pap\output_interproPAP40pc_aligned_to_acrA_withSuper.tsv')


# Automatically zoom to the loaded structures
pymol.cmd.zoom()

