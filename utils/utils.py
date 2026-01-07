from dataclasses import dataclass, field
from pathlib import Path
import shutil
import os
import sys

def check_dir(dir):

	if not os.path.exists(dir):
		print(f"\033[32m{dir} does not exist, it will be created\033[0m")
		os.makedirs(dir)
	else:
		print(f"\033[33m{dir} already exists, it will be overwritten\033[0m")
		shutil.rmtree(dir)
		os.makedirs(dir)

	return

def logger(message, level='INFO'):
	"""
	Function to log messages with different levels.
	Args:
		message (str): The message to log.
		level (str): The level of the message ('INFO', 'WARNING', 'ERROR').
	"""
	message = f"[{level}] {message}"
	if level == 'INFO':
		print(f"\033[32m{message}\033[0m")
	elif level == 'WARNING':
		print(f"\033[33m{message}\033[0m")
	elif level == 'ERROR':
		print(f"\033[31m{message}\033[0m")
		sys.exit(1)
	elif level == 'COMMAND':
		print(f"\033[35m{message}\033[0m")
	elif level == 'DEBUG':
		print(f"\033[34m{message}\033[0m")
	elif level == 'PAUSE':
		input(f"\033[36m{message}\n{level}: Press Enter to continue.\033[0m")
	else:
		print(f"\033[37m{message}\033[0m")  # Default to white for unknown levels

def make_dir_root_file(directory, file):
    if not file.GetDirectory(directory):
        file.mkdir(directory)
        logger(f"Created directory {directory} in file {file.GetName()}", level='WARNING')
    else:
        logger(f"Directory {directory} already exists in file {file.GetName()}", level='WARNING')


def enforce_list(var):
	"""
	Ensure that the input variable is a list.

	Parameters:
	- var (any): The input variable.

	Returns:
	- list: The input variable as a list.
	"""
	if not isinstance(var, list):
		var = [var]
	return var

@dataclass
class BinInfoBase:
    """The template - contains all the logic."""
    mins: list[float]
    maxs: list[float]
    edges: list[float] = field(init=False, repr=True)

    def __post_init__(self):
        if len(self.mins) != len(self.maxs):
            raise ValueError("mins and maxs must have the same length.")
        self.edges = self.mins + ([self.maxs[-1]] if self.maxs else [])

@dataclass
class CentralityInfo(BinInfoBase):
    """Identical to BinInfoBase, but explicitly for Centrality."""
    pass

@dataclass
class PtInfo(BinInfoBase):
    """Identical to BinInfoBase, but explicitly for Pt."""
    pass


def pt_folders_exist(folder_name, pt_cuts):
    """
    Check if pt folder structure (as created in preprocessing) exists.
      	It should contain folders named as 'pt_<pt_min>_<pt_max>'.

    Args:
    - folder_name (str): The name of the folder.
    - cutset (list): List of cutsets to check for pt folders.

    Returns:
    - bool: True if all pt folders exist, False otherwise.
    """
    folder_pt_mins = []  # list of pt mins in the folders
    folder_pt_maxs = []  # list of pt maxs in the folders

    folder_path = Path(folder_name)
    for folder in folder_path.iterdir():
        if folder.is_dir() and not folder.name.startswith("pt_"):
            continue
        suffix = folder.name[3:]  # remove 'pt_' prefix
        pt_min_str, pt_max_str = suffix.split('_')[0:2]
        folder_pt_mins.append(int(pt_min_str) / 10.0)
        folder_pt_maxs.append(int(pt_max_str) / 10.0)

    folder_pt_intervals = set(zip(folder_pt_mins, folder_pt_maxs))
    for pt_min, pt_max in zip(pt_cuts[0], pt_cuts[1]):
        if (pt_min, pt_max) not in folder_pt_intervals:
            return False
    return True


def get_root_files_in_directory(directory):
    """
    Get all ROOT files in the specified directory.

    Args:
    - directory (str): The directory to search for ROOT files.

    Returns:
    - root_files (list): A list of ROOT file paths.
    """
    path = Path(directory)
    root_files = [str(f) for f in sorted(path.glob('*.root'))]
    return root_files
