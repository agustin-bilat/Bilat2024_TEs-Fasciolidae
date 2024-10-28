#!/usr/bin/env python3

import os
import numpy as np
from scipy import stats
from ete3 import Tree
from Bio import SeqIO
import argparse

def get_ancestors(tree):
    return list(set(leaf.up for leaf in tree.iter_leaves()))

def get_terminals(treefile):
    try:
        tetree = Tree(treefile)
    except:
        print(f'{treefile} could not be read')
        return
    ancestors = get_ancestors(tetree)
    erroneous = []
    for node in ancestors:
        pair = node.children
        dists = {pair[0].dist: pair[0], pair[1].dist: pair[1]}
        mind, maxd = min(dists.keys()), max(dists.keys())
        if (mind == 0 and maxd > 0.1):
            shorter = dists[mind]
            if shorter.name != '':
                erroneous.append(shorter.name)
        elif mind > 0 and maxd/mind > 100:
            shorter = dists[mind]
            if shorter.name != '':
                erroneous.append(shorter.name)
    prune_list = []
    for leaf in tetree.iter_leaves():
        if leaf.name not in erroneous:
            prune_list.append(leaf.name)
    tetree.prune(prune_list)
    return {leaf.name: leaf.dist for leaf in tetree.iter_leaves()}

def summarise_trees(treedir, outfilename,species):
    with open(outfilename, 'w') as outfile:
        outfile.write('species\ttename\tmedlen\tnumbranches\n')
        for filename in os.listdir(treedir):
            if not filename.endswith('.nwk'):
                continue
            tename = filename.split('.')[0]
            terminals = get_terminals(f'{treedir}/{filename}')
            if not terminals:
                continue
            elif len(terminals) < 5:
                continue
            blens = list(terminals.values())
            bmed = np.median(blens)
            numbranches = len(blens)

            line = [species,tename, f"{bmed:.3f}", numbranches]
            line = '\t'.join([str(i) for i in line])
            outfile.write(f'{line}\n')

def main(treedir,species):
	outfilename = f'{species}_trees_summary.txt' # Use the name of either the species or assembly as outfile prefix
	summarise_trees(treedir, outfilename, species)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="For each tree get the median value of the lengths of terminal branches")
    parser.add_argument(
        '-d', '--directory',
        type=str,
        default=os.getcwd(),  # Use current working directory as default
        help='Path to the directory containing tree files in newick format'
    )
    parser.add_argument(
        '-s', '--species',
        type=str,
        required=True,  # Make this argument mandatory
        help='Prefix name of the species or assembly'
    )
    
    args = parser.parse_args()
    main(args.directory,args.species)
