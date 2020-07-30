from skbio import TreeNode


def root_outgroup(tree: TreeNode, node: str,
                  remove_outgroup: bool = True) -> TreeNode:
    fake_node = _insert_fake_parent(tree.find(node))
    rooted = tree.root_at(fake_node)

    new_subtree = rooted.find(fake_node.name)
    new_subtree.name = None

    return new_subtree if remove_outgroup else rooted


def _insert_fake_parent(node):
    real_parent = node.parent
    # We don't actually care what this node is called because it will become an
    # unnamed inner node. This name was chosen to avoid collisions
    fake_node = TreeNode(name='~~!ROOT_BRANCH!~~')
    for idx, child in enumerate(real_parent.children):
        if child is node:
            break

    del real_parent.children[idx]
    real_parent.children.append(fake_node)
    fake_node.children.append(node)
    node.parent = fake_node
    fake_node.parent = real_parent
    return fake_node
