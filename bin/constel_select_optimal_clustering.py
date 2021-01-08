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
import json
import numpy as np

from soma import aims

# ---------------------------Command-line--------------------------------------

doc = """Select the optimal number of clusters according to a reference
file"""


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
        'region',
        type=str,
        help='Region of interest')
    parser.add_argument(
        'optimal_cluster_numbers',
        type=str,
        help='The optimal numbers of clusters to apply')
    parser.add_argument(
        'optimal_clustering',
        type=str,
        help='Single optimal clustering')

    args = parser.parse_args(argv)

    return parser, args

# ----------------------------Main program-------------------------------------


def main():
    """
    """

    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    with open(args.optimal_cluster_numbers, 'r') as f:
        dict_optimal = json.load(f)

    opt_nb_clusters = int(dict_optimal[args.region])

    ind_clustering = aims.read(args.individual_clustering)
    individual_clustering = np.asarray(ind_clustering)

    # Select the optimal clustering
    opt_clustering = individual_clustering[opt_nb_clusters - 2]

    optimal_clustering = aims.TimeTexture_S16()
    optimal_clustering[0].assign(opt_clustering)

    aims.write(optimal_clustering, args.optimal_clustering)


if __name__ == "__main__":
    main()
