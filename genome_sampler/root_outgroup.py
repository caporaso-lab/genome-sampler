from skbio import TreeNode


# Should node be a string?
def root_outgroup(tree: TreeNode, node: str) -> TreeNode:
    # Since the outgroup node passed will always be a tip it will always have
    # exactly one neighbor.
    # TODO: Maybe we make an explicit check and error on more than one
    # neighbor?
    return tree.root_at(tree.find(node).neighbors[0])
