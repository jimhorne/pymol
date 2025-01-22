from pymol import cmd

# Name of the old and new objects
old_object = 'multi_state'
new_object = 'new_object'

# List or range of states to copy
states_to_copy = range(1, 24)  # This corresponds to states 1 to 23

# Loop over each state and create the new object
for i, state in enumerate(states_to_copy):
    target_state = i + 1  # States in PyMOL start from 1, not 0
    if target_state == 1:
        cmd.create(new_object, old_object, state, target_state)
    else:
        cmd.create(new_object, old_object, state, target_state)
