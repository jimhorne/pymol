# Get a list of all objects in the session
all_objects = cmd.get_names('objects')

# Get the current visibility status of all objects as dict
# Each value is a list where the index0 = enabled/disabled view, index1 = ?, index2 = ?, index3 = ?
view_dict = viewing.get_vis()

# Loop through all objects
for obj in all_objects:
    # Check if the object is visible
    visibility_flag = view_dict[obj][0]
    print(f'Obj {obj} is set to {visibility_flag}')
    # If object not visible, delete
    if visibility_flag == 0:
        cmd.delete(obj)