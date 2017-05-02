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
import sys
import argparse
import textwrap
import subprocess
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

            -------------------------------------------------------------------
            """))

    parser.add_argument("clusters", type=str, help="")
    parser.add_argument("white_mesh", type=str, help="")
    parser.add_argument("kmax", type=int, help="")

    return parser, parser.parse_args(argv)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)

    # mesh information
    mesh = args.white_mesh
    vertices = numpy.asarray(mesh.vertex())
    polygons = numpy.asarray(mesh.polygon())
    
    ftable = []
    # calculate the number of vertex by cluster
    for k in range(2, args.kmax + 1):
        # clusters information
        parcellation = aims.read(args.clusters)
        clusters = parcellation[k - 2].arraydata()
        labels = numpy.unique(clusters)
        # number of vertices
        measures = []
        for label in labels.nonzero()[0]:
            num_vert = len(numpy.where(clusters == label)[0])
            measures.append(num_vert)
        ftable.append(measures)
    

#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
