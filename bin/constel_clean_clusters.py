###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################


#----------------------------Imports-------------------------------------------

import argparse
import textwrap
import sys
import numpy as np

from soma import aims

# ---------------------------Command-line--------------------------------------

doc = """Merge small regions into bigger ones according to connectivity
profile"""


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(doc))

    parser.add_argument(
        'individual_clustering',
        type=str,
        help='List of clustering with different numbers of clusters')
    parser.add_argument(
        'reduced_individual_matrix',
        type=str,
        help='Reduced individual matrix')
    parser.add_argument(
        'group_clustering',
        type=str,
        help='Group clustering')
    parser.add_argument(
        'atlas_matrix',
        type=str,
        help='Atlas matrix')
    parser.add_argument(
        'min_size',
        type=int,
        help='Minimal size of regions')
    parser.add_argument(
        'cleaned_clustering',
        type=str,
        help='Cleaned clustering')

    args = parser.parse_args(argv)

    return parser, args

# ----------------------------Main program-------------------------------------


def main():
    """
    """

    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    individual_clustering = aims.read(args.individual_clustering)
    group_clustering = aims.read(args.group_clustering)
    cleaned_clustering = aims.TimeTexture_S16()

    # read the individual reduced matrix, shape (region vertices, basins)
    vol = aims.read(args.reduced_individual_matrix)
    rmatrix = np.asarray(vol)[:, :, 0, 0]
    if rmatrix.shape[0] < rmatrix.shape[1]:
        rmatrix = rmatrix.T

    # read the group reduced matrix, shape (region vertices, basins)
    atlas_matrix = aims.read(args.atlas_matrix)
    amatrix = np.asarray(atlas_matrix)[:, :, 0, 0]
    if amatrix.shape[0] < amatrix.shape[1]:
        amatrix = amatrix.T
    for k in range(len(individual_clustering)):
        clustering = individual_clustering[k].np
        gclustering = group_clustering[k].np
        labels, counts = np.unique(clustering[clustering != 0],
                                   return_counts=True)
        merged_labels = [0]
        print(labels, counts)
        min = np.min(counts)
        while min < args.min_size:
            label = labels[np.argmin(counts)]
            merged_labels.append(label)
            # clusters = clustering[clustering != 0]
            indices = np.where(clustering == label)[0]
            # Research new labels
            gclustering = gclustering[gclustering != 0]
            amat = amatrix[~np.isin(gclustering, merged_labels)]
            gclusters = gclustering[~np.isin(gclustering, merged_labels)]
            for i, idx in enumerate(indices):
                distance = np.sum((rmatrix[i] - amat)**2, axis=1)
                new_label = gclusters[np.argmin(distance)]
                clustering[idx] = new_label
            labels, counts = np.unique(clustering[clustering != 0],
                                       return_counts=True)
            print(labels, counts)
            min = np.min(counts)
        cleaned_clustering[k].assign(clustering)
    aims.write(cleaned_clustering, args.cleaned_clustering)


if __name__ == "__main__":
    main()
