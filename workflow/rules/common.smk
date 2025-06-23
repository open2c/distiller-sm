# Copied from https://github.com/mirnylab/distiller-sm/blob/master/_distiller_common.py
import os, pathlib
import numpy as np
import shlex

def argstring_to_dict(argstring):
    """
    Convert a command line argument string into a dictionary.
    
    Args:
        argstring (str): The command line argument string.
        
    Returns:
        dict: A dictionary with argument names as keys and their values.
        
    Note:
        - Arguments that start with '-' are considered keys.
        - If an argument has no value, it is set to True.
        - Values can be separated by spaces and will be joined into a single string.
        - If the final argument is a value without a preceding key, it will be included as a value for the last key (issue in case of specified input as last argument).
    """
    args = shlex.split(argstring)
    keys = np.where([arg.startswith('-') for arg in args])[0]
    vals = np.where([not arg.startswith('-') for arg in args])[0]
    args_arrs = [arr for arr in np.split(args, keys) if arr.size > 0]
    argdict = {str(arr[0]): (' '.join(arr[1:]) if arr.size > 1 else True) for arr in args_arrs }
    return argdict

def needs_downloading(fastq_files, side):
    if len(fastq_files) == 1 and fastq_files[0].startswith("sra:"):
        return True
    elif (
        fastq_files[side].startswith("sra:")
        or fastq_files[side].startswith("http://")
        or fastq_files[side].startswith("ftp://")
    ):
        return True
    else:
        return False


def parse_fastq_folder(root_folder_path):
    library_run_fastqs = {}
    fastq_root_folder = pathlib.Path(root_folder_path)
    for fastq_file in sorted(fastq_root_folder.glob("**/*")):
        if not fastq_file.is_file():
            continue
        split_path = fastq_file.relative_to(fastq_root_folder).parts
        if len(split_path) != 3:
            raise Exception(
                "The fastq folder has a non-standard structure! "
                "Please, make sure that the fastq folders only contains a folder "
                "per library, each containing a folder per run, each _only_ containing "
                "two fastq(.gz) files."
            )
        library, run, fastq = split_path

        library_run_fastqs.setdefault(library, {}).setdefault(run, []).append(
            str(fastq_file.absolute())
        )

    return library_run_fastqs


def check_fastq_dict_structure(library_run_fastqs):
    for library, run_dict in library_run_fastqs.items():
        if not isinstance(run_dict, dict):
            return False
        for run, fastq_files in run_dict.items():
            if not isinstance(fastq_files, list) or (len(fastq_files) > 2):
                return False
    return True


def organize_fastqs(config):
    if isinstance(config["input"]["raw_reads_paths"], str):
        library_run_fastqs = parse_fastq_folder(config["input"]["raw_reads_paths"])
    elif isinstance(config["input"]["raw_reads_paths"], dict):
        if not check_fastq_dict_structure(config["input"]["raw_reads_paths"]):
            raise Exception(
                "An unknown format for library_fastqs! Please provide it as either "
                'a path to the folder structured as "library/run/fastqs" or '
                "a dictionary specifying the project structure."
            )
        library_run_fastqs = config["input"]["raw_reads_paths"]
    else:
        raise Exception(
            "An unknown format for library_fastqs! Please provide it as either "
            "a path to the folder with the structure library/run/fastqs or "
            "a dictionary specifying the project structure."
        )

    return library_run_fastqs
