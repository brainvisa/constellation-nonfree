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
* normalize the individual connectivity profiles
*
*

Main dependencies: PyAims library

Author: Sandrine Lefranc
"""


#----------------------------Imports-------------------------------------------


# python system module
import sys
import numpy
import argparse
import textwrap

# aims module
from soma import aims


#----------------------------Functions-----------------------------------------


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument("mask", type=str, help="")
    parser.add_argument("meanprofile", type=str, help="")
    parser.add_argument("normprofile", type=str, help="")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1]), sys.argv[2]), sys.argv[3])
    parser, args = parse_args(arguments)

    # read the files with aims
    aims_mask = aims.read(args.mask)
    aims_profile = aims.read(args.meanprofile)
    
    # transform the aims files to numpy array
    np_profile = aims_profile[0].arraydata()
    np_mask = aims_mask[0].arraydata()
    
    # create the connectivity profile overlap the mask
    np_profile[numpy.where(np_mask == 0)] = 0

    # normalize of the profile by its max of connections
    if np_profile.max() > 0:
        np_profile *= 1. / np_profile.max()

    # create a time texture object by assigning the mean profile
    aimstex = aims.TimeTexture_FLOAT()
    aimstex[0].assign(np_profile)

    # write the file on the disk
    aims.write(aimstex, args.normprofile)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
