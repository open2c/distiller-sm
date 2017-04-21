import os, pathlib


def parse_fastq_folder(root_folder_path):
    library_run_fastqs = {}
    fastq_root_folder = pathlib.Path(root_folder_path)
    for fastq_file in sorted(fastq_root_folder.glob('**/*')):
        if not fastq_file.is_file():
            continue
        split_path = fastq_file.relative_to(fastq_root_folder).parts
        if len(split_path) != 3:
            raise Exception(
                'The fastq folder has a non-standard structure! '
                'Please, make sure that the fastq folders only contains a folder '
                'per library, each containing a folder per run, each _only_ containing '  
                'two fastq(.gz) files.')
        library, run, fastq = split_path

        library_run_fastqs.setdefault(library, {}).setdefault(run, []).append(
            str(fastq_file.absolute()))

    return library_run_fastqs


def check_fastq_dict_structure(library_run_fastqs):
    for library, run_dict in library_run_fastqs.items():
        if not isinstance(run_dict, dict):
            return False
        for run, fastq_files in run_dict.items():
            if (not isinstance(fastq_list, list) or (len(fastq_list) != 2)):
                return False
    return True


def organize_fastqs(config):
    if isinstance(config['fastq_paths'], str):
        library_run_fastqs = parse_fastq_folder(config['fastq_paths'])
    elif isinstance(config['fastqs'], dict):
        if not check_fastq_dict_structure(config['fastq_paths']):
            raise Exception(
                'An unknown format for library_fastqs! Please provide it as either '
                'a path to the folder structured as "library/run/fastqs" or '
                'a dictionary specifying the project structure.')
        library_run_fastqs = config['fastq_paths']
    else:
        raise Exception(
            'An unknown format for library_fastqs! Please provide it as either '
            'a path to the folder with the structure library/run/fastqs or '
            'a dictionary specifying the project structure.')

    return library_run_fastqs 

