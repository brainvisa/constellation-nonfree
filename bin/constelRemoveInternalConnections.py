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


# python system modules
import numpy
import textwrap
import argparse
import sys
import json

# constel modules
from constel.lib.texturetools import management_internal_connections
from constel.lib.texturetools import normalize_profile

# soma module
from soma import aims


#----------------------------Functions-----------------------------------------


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Remove the patch internal connections of the cortical connectivity
            profile. Then, the profile is normalized by its total number of
            connections.
            -------------------------------------------------------------------
            """))

    #adding arguments
    parser.add_argument("label", type=str,
                        help="a number of ROI")
    parser.add_argument("profile", type=str,
                        help="a coritcal connectivity profile")
    parser.add_argument("gyriseg", type=str,
                        help="a labeling of gyri cortical segmentation")
    parser.add_argument("normprofile", type=str,
                        help="normalized connectivity profile")
    parser.add_argument("-c", action="store_true",
                        dest="internal_connections")
    parser.add_argument("-q", action="store_false",
                        dest="internal_connections")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # Load the arguments of parser (delete script name: sys.arg[0])
    print 
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)
    
    # load files
    aims_mask = aims.read(args.gyriseg)
    numpy_mask = aims_mask[0].arraydata()
    aims_profile = aims.read(args.profile)
    numpy_profile = aims_profile[0].arraydata()
    
    new_profile = management_internal_connections(
        args.label, numpy_mask, numpy_profile,
        args.internal_connections)

    norm_profile = normalize_profile(new_profile)
    
    tex = aims.TimeTexture_S16()
    tex[0].assign(norm_profile)
    
    aims.write(tex, args.normprofile)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
