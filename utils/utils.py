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


