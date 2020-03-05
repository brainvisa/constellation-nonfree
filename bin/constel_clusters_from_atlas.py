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

# soma
from soma import aims

# constel
from constel.lib.clustering.clusterstools import nearest_neighbour_profiles
from six.moves import range


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
    parser.add_argument(
        "individual_patch_texture",
        nargs='?',
        help="Individual patch texture for the cortical region. Default: use atlas one."
        )
    parser.add_argument(
        "individual_patch_num",
        nargs='?',
        type=int,
        help="Individual patch number identifier in texture. Default: use atlas."
        )

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)

    # read the individual reduced matrix, shape (region vertices, basins)
    rmatrix = aims.read(args.reduce_matrix)
    rmat = numpy.asarray(rmatrix)[:, :, 0, 0]
    if rmat.shape[0] < rmat.shape[1]:
        rmat = numpy.transpose(rmat)
    print("rmat: ", rmat.shape)

    # read the group reduced matrix, shape (region vertices, basins)
    atlas = aims.read(args.atlas_matrix)
    amat = numpy.asarray(atlas)[:, :, 0, 0]
    if amat.shape[0] < amat.shape[1]:
        amat = numpy.transpose(amat)
    print("amat: ", amat.shape)

    # read the clustering texture of the cortical region
    # (e.g. gyrus of FreeSurfer)
    atlas_clusters = aims.read(args.atlas_classes)

    iindices = None
    if args.individual_patch_texture is not None:
        if args.individual_patch_num is None:
            raise ValueError(
                'individual_patch_num must be specified along with '
                'individual_patch_texture')
        ipatch = aims.read(args.individual_patch_texture)
        iindices = numpy.where(ipatch[0].arraydata()
                               == args.individual_patch_num)[0]
        ntex = len(ipatch[0])
        if len(iindices) != rmat.shape[0]:
            raise ValueError(
                'individual matrix and texture do not match. Matrix shape is: '
                '%d rows, while ind. texture is %s'
                % (rmat.shape[0], len(iindices)))
    else:
        ntex = len(atlas_clusters[0])

    # define the nearest neighbour profiles
    indiv_clusters = atlas_clusters.__class__()
    for i in range(len(atlas_clusters)):
        aclusters = atlas_clusters[i].arraydata()
        indices = numpy.where(aclusters != 0)[0]
        if iindices is None:
            icindices = indices
        else:
            icindices = iindices
        res_classes = nearest_neighbour_profiles(
            rmat,
            amat,
            aclusters,
            indices,
            icindices,
            ntex)
        indiv_clusters[i].assign(res_classes)
    aims.write(indiv_clusters, args.individual_classes)

if __name__ == "__main__":
    main()
