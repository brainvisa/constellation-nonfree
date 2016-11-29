#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

"""
This script does the following:
*

Main dependencies: PyAims library

Author: Sandrine Lefranc, 2015
"""


#----------------------------Imports-------------------------------------------


# system module
from __future__ import print_function
import argparse
import textwrap
import numpy
import sys

# soma
from soma import aims

# constel
from constel.lib.clustering.clusterstools import nearest_neighbour_profiles


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Compute the nearest neighbours profiles.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        "reduce_matrix",
        help="Individual reduced matrix.")
    parser.add_argument(
        "atlas_matrix",
        help="Group reduced matrix.")
    parser.add_argument(
        "atlas_classes",
        help="Clustering of the cortical region.")
    parser.add_argument(
        "individual_classes",
        help="Individual clustering of the cortical region defined from atlas."
        )

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)

    # read the individual reduced matrix, shape (region vertices, basins)
    rmatrix = aims.read(args.reduce_matrix)
    rmat = numpy.asarray(rmatrix)[:, :, 0, 0].T

    # read the group reduced matrix, shape (region vertices, basins)
    atlas = aims.read(args.atlas_matrix)
    amat = numpy.asarray(atlas)[:, :, 0, 0]

    # read the clustering texture of the cortical region
    # (e.g. gyrus of FreeSurfer)
    atlas_clusters = aims.read(args.atlas_classes)

    # define the nearest neighbour profiles
    indiv_clusters = atlas_clusters.__class__()
    for i in range(len(atlas_clusters)):
        aclusters = atlas_clusters[i].arraydata()
        indices = numpy.where(aclusters != 0)[0]
        res_classes = nearest_neighbour_profiles(
            rmat,
            amat,
            aclusters,
            indices)
        indiv_clusters[i].assign(res_classes)
    aims.write(indiv_clusters, args.individual_classes)

if __name__ == "__main__":
    main()
