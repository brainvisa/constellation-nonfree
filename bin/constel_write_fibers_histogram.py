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
"""


#----------------------------Imports-------------------------------------------

# python system module
from __future__ import absolute_import
import sys
import argparse
import textwrap
import soma.subprocess
import matplotlib.pyplot as plt

# aims module
from soma import aims


#----------------------------Functions-----------------------------------------

def parse_args(argv):
    """Parses the given list of arguments."""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Write the fibers lengths in a text file.
            -------------------------------------------------------------------
            """))

    parser.add_argument("fiber_tracts", type=str, help="")
    parser.add_argument("lengths_filename", type=str, help="")

    return parser, parser.parse_args(argv)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)

    # Give the options of the command
    cmd = ["constelFibersLengths",
           "-i", args.fiber_tracts,
           "-o", args.lengths_filename,
           "-verbose"]

    # Execute the command
    soma.subprocess.check_call(cmd)

    # Give the number of curves in the bundle file
    bundles = aims.read(args.fiber_tracts)
    curves_count = bundles["curves_count"]

    # Give the number of lengths in the text file
    # Create a list with all fiber lengths
    lengths = []
    with open(args.lengths_filename, "r") as lengthsfile:
        n = 0
        for line in lengthsfile:
            n += 1
            lengths.append(int(line.rstrip()))

    # Check
    if curves_count != n:
        raise ValueError(
            "The number of fibers is different of the number of lengths.")

    minlength = min(lengths)
    maxlength = max(lengths)

    plt.hist(sorted(lengths),
             range=(0, 300),
             bins=300,
             color='yellow',
             edgecolor='red')
    plt.xlabel("Fiber lengths (mm)")
    plt.ylabel("Number of fiber tracts")
    title = "{0} fiber tracts, {1} < lengths < {2} mm.".format(
        int(curves_count), minlength, maxlength)
    plt.title(title)

    filename = args.lengths_filename.split('.')[0] + "_histogram"
    plt.savefig(filename, dpi=350)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
