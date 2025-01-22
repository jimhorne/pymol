# PyMOL script to align all loaded objects to the reference structure
from pymol import cmd

# Name of the reference object
#reference_name = 'reference_oriented_AF-P02930_tolCec'
reference_name = 'fold_2024_09_09_14_52_acra6acrb2_atp2_model_0'

# Get a list of all loaded objects in PyMOL
object_list = cmd.get_object_list()

# Specify chain to align against
chain = 'A'

# Align each object to the reference, skipping the reference itself
for obj in object_list:
    if obj != reference_name:
        # Perform alignment
        if chain and chain.isalpha(): # if chain variable has been assigned and is a letter 
            rmsd = cmd.align(f'{obj} & c. {chain}', f'{reference_name} & c. {chain}')
        else:
            rmsd = cmd.align(obj, reference_name)
        # Print the results
        print(f'Aligned {obj} to {reference_name} with RMSD: {rmsd[0]:.4f} Å')

# =============================================================================
# # Align each object to the reference using only the first 75 residues
# for obj in object_list:
#     if obj != reference_name:
#         residues_to_keep = []
#         
#         # Iterate over residues in the object
#         resi_count = 0
#         for model in cmd.get_model(obj).atom:
#             resi = model.resi
#             if resi not in residues_to_keep:
#                 residues_to_keep.append(resi)
#                 resi_count += 1
# 
#             # Stop if we've selected 150 residues
#             if resi_count >= 86:
#                 break
#         
#         # Join the list into a string for PyMOL selection
#         resi_selection = '+'.join(residues_to_keep)
#         rmsd = cmd.super(f'/{obj} and resi {resi_selection}', f'{reference_name} and resi 24-91')
#         
#         # Print the results
#         print(f'Aligned {obj} to {reference_name} with RMSD: {rmsd[0]:.4f} Å')
# =============================================================================


# Optionally, save the aligned structures
# for obj in object_list:
#     if obj != reference_name:
#         cmd.save(f'{obj}_aligned.pdb', obj)

print("Alignment completed for all objects.")
