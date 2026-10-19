#!/usr/bin/env python3
"""Find the labelled internal-node name that is the MRCA of a set of tips, in a
Newick tree whose internal nodes are already named (e.g. tree_doctor --name-ancestors,
or a .mod TREE: line). Pure stdlib. Prints the MRCA node label to stdout.
Usage: find_mrca.py <newick_or_mod_file> tip1 tip2 [tip3 ...]
"""
import sys, re

def read_newick(path):
    txt = open(path).read()
    # if a .mod file, grab the TREE: line
    m = re.search(r'TREE:\s*(.+)', txt)
    nwk = m.group(1) if m else txt
    nwk = nwk.strip().rstrip(';').strip()
    return nwk

class Node:
    __slots__ = ("name", "children", "tips")
    def __init__(self): self.name=None; self.children=[]; self.tips=set()

def parse(nwk):
    pos = 0
    def parse_clade():
        nonlocal pos
        node = Node()
        if nwk[pos] == '(':
            pos += 1
            while True:
                node.children.append(parse_clade())
                if nwk[pos] == ',':
                    pos += 1; continue
                if nwk[pos] == ')':
                    pos += 1; break
        # read label (name) up to , ) : ;
        m = re.match(r'[^,():;]*', nwk[pos:])
        lbl = m.group(0); pos += len(lbl)
        # strip branch length
        node.name = lbl.split(':')[0] if lbl else None
        # skip a branch length attached after a name-less internal node
        if pos < len(nwk) and nwk[pos] == ':':
            m2 = re.match(r':[0-9.eE+-]+', nwk[pos:]); pos += len(m2.group(0))
        if not node.children:
            node.tips = {node.name}
        else:
            for c in node.children: node.tips |= c.tips
        return node
    return parse_clade()

def find_mrca(node, targets):
    # smallest node whose tips superset targets
    if not targets <= node.tips:
        return None
    for c in node.children:
        r = find_mrca(c, targets)
        if r is not None:
            return r
    return node

def main():
    path = sys.argv[1]; tips = set(sys.argv[2:])
    root = parse(read_newick(path))
    missing = tips - root.tips
    if missing:
        sys.stderr.write(f"WARN tips not in tree: {sorted(missing)}\n")
        tips &= root.tips
    node = find_mrca(root, tips)
    if node is None or node.name is None:
        sys.stderr.write("ERROR: MRCA not found or unnamed (did you --name-ancestors?)\n")
        sys.exit(2)
    print(node.name)

if __name__ == "__main__":
    main()
