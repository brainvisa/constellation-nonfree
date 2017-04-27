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
*
*

Main dependencies: PyAims library

Author: Sandrine Lefranc
"""


#----------------------------Imports-------------------------------------------


# python system module
from __future__ import print_function
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
    parser.add_argument("watershed", type=str, help="")
    parser.add_argument("filterwatershed", type=str, help="")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)

    # Low connections to gyrus : filtered watershed with "minVertex_nb"
    minVertex_nb = 20
    
    # read the file with aims
    basins_tex = aims.read(args.watershed)
    
    # transform the aims files to numpy array
    basinTex_ar = basins_tex[0].arraydata()
    
    # recover the list of labels that are applied
    basins_labels = numpy.unique(basinTex_ar)
    
    # criteria: replace the basins with less than 20 vertices by zero
    for basin_label in basins_labels:
        if numpy.where(basinTex_ar == basin_label)[0].size < minVertex_nb:
            basinTex_ar[basinTex_ar == basin_label] = 0
    
    # keep only the labels different of zero
    otex_kept_labels = numpy.unique(basinTex_ar)
    if list(otex_kept_labels).count(0) != 0:
        tmp = list(otex_kept_labels)
        tmp.remove(0)
        otex_kept_labels = numpy.array(tmp)
    
    # renumbers all the labels from 1
    otex = basinTex_ar
    for i in range(len(otex_kept_labels)):
        current_label = otex_kept_labels[i]
        otex[basinTex_ar == current_label] = i + 1

    # create a time texture object
    aimstex = aims.TimeTexture_S16()
    aimstex[0].assign(otex)

    #filteredBasins = remove_labels(basins_tex, labelsToRemove_list)
    
    # write the file on the disk
    aims.write(aimstex, args.filterwatershed)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
