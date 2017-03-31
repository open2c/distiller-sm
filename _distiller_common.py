def iflatten_tree(tree, root=None):
    for k, v in tree.items():
        if isinstance(v, dict):
            for z in iflatten_tree(v, root+[k] if root is not None else [k]):
                yield z
        else:
            yield (tuple(root+[k] if root is not None else [k]), v)

def flatten_tree(tree):
    return dict(list(iflatten_tree(tree)))
