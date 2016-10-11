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


# python system modules
from __future__ import print_function
import sys
import numpy
import optparse

# soma
from soma import aims

# constel
from constel.lib.utils.matrixtools import resize_matrix


#----------------------------Functions-----------------------------------------


def parseOpts(argv):
    desc = """usage: %prog [options] filename"
    Calculate concatenated or averaged matrix.
    """

    parser = optparse.OptionParser(desc)

    parser.add_option("-m", "--matrices", 
                      dest="list_matrices",
                      action="append", 
                      help="list of reduced matrices")
    parser.add_option("-o", "--output", 
                      dest="matrix", 
                      help="output matrix (averaged or concatenated)")
    parser.add_option("-s", "--study", 
                      dest="study",
                      help="Choice 'avg' or 'concat'")

    return parser, parser.parse_args(argv)


#----------------------------Main program--------------------------------------


def main():
    parser, (options, args) = parseOpts(sys.argv)
 
    # two cases: 
    # (1) concatenated matrix: M(vertices_patch_concatenated, basins)
    # (2) averaged matrix: M(vertices_patch, basins)
    if options.study == 'concat':
        shapex = None
        shapey = 0
        for matrix in options.list_matrices:
            print(matrix)
            reduced_matrix = aims.read(matrix)
            reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]
            if shapex is None:
                shapex = reduced_matrix.shape[0]
            shapey += reduced_matrix.shape[1]
        # generate a list of asarray (data_type: DOUBLE)
        l = numpy.zeros((shapex, shapey), dtype=numpy.float32)
        shifty = 0
        for matrix in options.list_matrices:
            reduced_matrix = aims.read(matrix)
            reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]
            l[:, shifty:(shifty+reduced_matrix.shape[1])] = reduced_matrix
            shifty += reduced_matrix.shape[1]
        print("1: ", l.shape)
        l = resize_matrix(aims.Volume(l))
        print("2: ", l.getSize())
        aims.write(l, options.matrix)
    else: #avg
        # generate a list of asarray
        list_matrices = []
        for matrix in options.list_matrices:
            reduced_matrix = aims.read(matrix)
            reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]
            list_matrices.append(reduced_matrix.astype('float32'))
            print(reduced_matrix.shape)
        sum_matrix = [sum(i) for i in zip(*list_matrices)]
        averaged_matrix = numpy.array(sum_matrix) / len(list_matrices)
        print(averaged_matrix.shape)
        averaged_matrix = resize_matrix(aims.Volume(averaged_matrix))
        aims.write(averaged_matrix, options.matrix)

if __name__ == "__main__":
    main()
