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
* calculate the average mask from the individual connectivity matrices
* according to two criteria:
          - a vertex is valid if it pass the cthreshold
          - 50 percent of subjects must possess the valid vertex

Main dependencies: PyAims library

Author: Sandrine Lefranc
"""


#----------------------------Imports-------------------------------------------


# python system module
from __future__ import print_function
from __future__ import absolute_import
import sys
import json
import numpy
import argparse
import textwrap

# aims module
from soma import aims


#----------------------------Functions-----------------------------------------


def mylist(string):
    return json.loads(string)


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            calculate an average mask from the individual connectivity profiles
            according to two criteria (arbitrary):
                - a vertex is valid if it pass the threshold
                - 50 percent of subjects must possess the valid vertex
            -------------------------------------------------------------------
            """))

    #adding arguments
    parser.add_argument(
        "profiles", type=mylist,
        help="list of individual connectivity profiles")
    parser.add_argument("mask", type=str,
                        help="output file: average mask")
    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # Load the arguments of parser (delete script name: sys.arg[0])
    arguments = (json.dumps(eval(sys.argv[1])), sys.argv[2])
    parser, args = parse_args(arguments)

    # load all the names of the individual connectivity profiles
    fnames = args.profiles

    for idx, fname in enumerate(fnames):
        # read the file with aims
        aimsprofile = aims.read(fname)

        # transform the aims file to numpy array
        np_profile = aimsprofile[0].arraydata()

        # initialize a numpy array, the same size as the individual profile
        if idx == 0:
            nvertex = np_profile.shape[0]
            fvertices = numpy.zeros((nvertex, ), dtype=numpy.int16)

        # first criteria: a vertex is valid if it pass the cthreshold
        cthreshold = np_profile.max() / 100.
        valid_vert = np_profile >= cthreshold
        fvertices[valid_vert] += 1

    # second criteria: 50 percent of subjects must possess the valid vertex
    # 1: keep the vertex in the mask
    sthreshold = round(len(fnames) / 2.)
    fvertices[numpy.where(fvertices < sthreshold)] = 0
    fvertices[numpy.where(fvertices >= sthreshold)] = 1

    # create a time texture object by assigning the valid vertices
    tex = aims.TimeTexture_S16()
    tex[0].assign(fvertices)

    # write the file (mask) on the disk
    aims.write(tex, args.mask)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
