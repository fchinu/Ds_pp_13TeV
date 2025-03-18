import os
import shutil

def delete_directories_without_log(root_folder):
    # Loop through all subdirectories in the root folder
    for dirpath, dirnames, filenames in os.walk(root_folder, topdown=False):
        # Check if 'charmHadr.txt' is not in the list of files
        if 'charmHadr.txt' not in filenames and dirpath != root_folder and "JobPbPb" in dirpath:
            pass
            print(f"Deleting directory: {dirpath}")
            # Remove the directory and its contents
            shutil.rmtree(dirpath)
        else:
            # # check if the charmHadr.txt is empty
            # if 'charmHadr.txt' in filenames:
            #     if os.stat(os.path.join(dirpath, 'charmHadr.txt')).st_size == 0:
            #         print(f"Deleting directory: {dirpath}")
            #         # Remove the directory and its contents
            #         shutil.rmtree(dirpath)
            # check if "--------" is not in the file (end of file)
            if 'charmHadr.txt' in filenames:
                with open(os.path.join(dirpath, 'charmHadr.txt'), 'r') as f:
                    if "--------" not in f.read():
                        print(f"Deleting directory: {dirpath}")
                        # Remove the directory and its contents
                        shutil.rmtree(dirpath)
# Specify the root folder to start from
root_folder = '/home/fchinu/Run3/Ds_pp_13TeV/PYTHIA_Simulations/pbpb'

# Call the function to delete directories without 'pythiaAA.log'
delete_directories_without_log(root_folder)