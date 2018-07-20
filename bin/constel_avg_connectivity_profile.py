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
* Create and write a group profile

Main dependencies: PyAims library

Author: Sandrine Lefranc, 2015
"""


#----------------------------Imports-------------------------------------------


# python system module
import sys
import json
import argparse
import textwrap
import shutil
import soma.subprocess

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
            Create a group profile. All individual profiles are added and each
            final value is then divided by the number of profile.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        "normprofiles", type=mylist, help="list of individual profiles")
    parser.add_argument(
        "meanprofile", type=str, help="output: group of profile")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    """Execute the command to create and write a group profile on the disk.
    """
    # Load the arguments of parser (delete script name: sys.arg[0])
    arguments = (json.dumps(eval(sys.argv[1])), sys.argv[2])
    parser, args = parse_args(arguments)

    # apply a formula to a set of homogeneous textures
    for idx, tex_filename in enumerate(args.normprofiles):
        if idx == 0:
            shutil.copy(str(tex_filename), str(args.meanprofile))
        else:
            cmd = ["cartoLinearComb.py",
                   "-f", "I1 + I2",
                   "-i", tex_filename,
                   "-i", args.meanprofile,
                   "-o", args.meanprofile]
            soma.subprocess.check_call(cmd)

    # read mean profile
    meanprofile = aims.read(args.meanprofile)

    # mean of the values
    for i in xrange(meanprofile.nItem()):
        val = meanprofile[0][i]
        meanprofile[0][i] = val/len(args.normprofiles)

    # write the file on the disk
    aims.write(meanprofile, args.meanprofile)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
