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
from __future__ import absolute_import
import argparse
import textwrap
import numpy
import sys
import os

from matplotlib import pyplot

# soma
from soma import aims

# constel
from soma.aims.texturetools import meshDiceIndex
from six.moves import range


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            BLABLABLA
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        "individual_clustering", type=str,
        help="Individual reduced matrix.")
    parser.add_argument(
        "atlas_clustering", type=str,
        help="Group reduced matrix.")
    parser.add_argument(
        "white_mesh", type=str,
        help="")
    parser.add_argument(
        "individual_time_step", type=int,
        help="")
    parser.add_argument(
        "atlas_time_step", type=int,
        help="")
    parser.add_argument(
        "cortical_region", type=str,
        help="")
    parser.add_argument(
        "ofile", type=str,
        help="Save the stem plot in a PDF file.")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)

    # read the individual clustering
    indiv_clusters = aims.read(args.individual_clustering)

    # read the group clustering
    atlas_clusters = aims.read(args.atlas_clustering)

    # read white mesh
    mesh = aims.read(args.white_mesh)

    # Compute the dice index between two clustering texture
    dice, areas1, areas2, inter = meshDiceIndex(
        mesh,
        indiv_clusters,
        atlas_clusters,
        timestep1=args.individual_time_step,
        timestep2=args.atlas_time_step)
    print(dice)

    nb_clusters = args.atlas_time_step + 2

    # The stem plot
    f = pyplot.figure()
    markerLines, stemLines, baseLines = pyplot.stem(
        list(range(1, nb_clusters + 1)),
        dice[1:])
    pyplot.setp(markerLines,
                color="cyan",
                markersize=10,
                markeredgecolor="blue",
                markeredgewidth=0.5)
    pyplot.setp(stemLines,
                color="green",
                linewidth=1,
                linestyle="dashed")
    pyplot.setp(baseLines,
                color="red",
                linewidth=2,
                linestyle="solid")
    pyplot.margins(0.1, 0.1)
    pyplot.ylabel("Dice index", fontsize=18)
    pyplot.xlabel("Clusters (k)", fontsize=18)
    pyplot.title(
        "Parcellation with {0} clusters for the cortical region: '{1}'".format(
            nb_clusters,
            args.cortical_region), fontsize=18)
    f.savefig(args.ofile)


if __name__ == "__main__":
    main()

