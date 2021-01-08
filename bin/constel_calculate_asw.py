#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################


#----------------------------Imports-------------------------------------------

from __future__ import print_function

# system module
from __future__ import absolute_import
import sys
import json
import numpy
import argparse
import textwrap

from six.moves import range

from sklearn.metrics import silhouette_score

# soma module
from soma import aims


#----------------------------Functions-----------------------------------------

# Command parameters
doc = """
Average Silhouette Width
~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the average silhouette width in order to find the optimal number of
clusters.
"""


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(doc))

    parser.add_argument(
        "matrix", type=str,
        help="Individual reduced matrix.")
    parser.add_argument(
        "clusters", type=str,
        help="Clustering")
    parser.add_argument(
        "kmax", type=int,
        help="The max number of clusters.")
    parser.add_argument(
        "silhouette_file", type=str,
        help="JSON file storing silhouette average scores for each clustering")
    parser.add_argument(
        "--kmin", type=int, default=2,
        help="Minimum number of clusters taken into account.")

    # parsing arguments
    args = parser.parse_args(argv)

    # verification of kmin and kmax values
    if args.kmin > args.kmax:
        raise ValueError("kmin should be lower than kmax")
    elif args.kmin < 2:
        raise ValueError('There should be at least 2 clusters (kmin >= 2)')

    return parser, args

# ----------------------------Main program-------------------------------------


def main():
    """
    """

    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    dict_clusters = {}

    aims_red_mat = aims.read(args.matrix)
    reduced_matrix = numpy.asarray(aims_red_mat)[:, :, 0, 0]

    aims_clusters = aims.read(args.clusters)
    cluster_list = numpy.asarray(aims_clusters)

    for k in range(args.kmin, args.kmax+1):
        clusterid = cluster_list[k-2][numpy.where(cluster_list[k-2] != 0)]
        if reduced_matrix.shape[0] != len(clusterid):
            reduced_matrix = reduced_matrix.T
        dict_clusters[k] = str(silhouette_score(reduced_matrix, clusterid))
    print(dict_clusters)

    with open(args.silhouette_file, "w") as outfile:
        json.dump(dict_clusters, outfile)


if __name__ == "__main__":
    main()
