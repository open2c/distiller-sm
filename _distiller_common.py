import os, pathlib


def iflatten_tree(tree, root=None):
    for k, v in tree.items():
        if isinstance(v, dict):
            for z in iflatten_tree(v, root+[k] if root is not None else [k]):
                yield z
        else:
            yield (tuple(root+[k] if root is not None else [k]), v)


def flatten_tree(tree):
    return dict(list(iflatten_tree(tree)))


def organize_fastqs(config):
    if isinstance(config['fastq'], str):
        FLAT_TREE = {}
        for p in (pathlib.Path('.')/config['fastq']).glob('**/*'):
            if p.is_file():
                k = p.parts[1:-1]
                FLAT_TREE[k] = FLAT_TREE.get(k, []) + [str(p)]

        FLAT_TREE = {k:sorted(v) for k,v in FLAT_TREE.items()}
    else:
        FLAT_TREE = flatten_tree(config['fastq'])

    RUN_TO_FASTQS = {config['library_name_sep'].join(k):v 
                    for k,v in FLAT_TREE.items()}

    RUN_FULL_NAMES = list(RUN_TO_FASTQS.keys())

    LIBRARY_TO_FASTQS = {exp:[config['library_name_sep'].join(k)
                                 for k in FLAT_TREE if k[0] == exp]
                            for exp in [k[0] for k in FLAT_TREE] }

    LIBRARY_NAMES = list(LIBRARY_TO_FASTQS.keys())

    return RUN_TO_FASTQS, RUN_FULL_NAMES, LIBRARY_TO_FASTQS, LIBRARY_NAMES 

