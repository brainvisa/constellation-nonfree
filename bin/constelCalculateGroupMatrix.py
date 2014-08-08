#!/usr/bin/env python

# python system modules
import sys
import numpy
import optparse

# soma
from soma import aims

# Constel
from constel.lib.connmatrix.connmatrixtools import matrix_converter


def parseOpts(argv):
    desc = """usage: %prog [options] filename"
    Calculate concatenated or averaged matrix.
    """

    parser = optparse.OptionParser(desc)

    parser.add_option("-m", "--matrices", 
                      dest="list_matrices",
                      action="append", 
                      help="list of matrices")
    parser.add_option("-o", "--output", 
                      dest="matrix", 
                      help="output matrix (averaged or concatenated)")
    parser.add_option("-s", "--study", 
                      dest="study",
                      help="Choice 'avg' or 'concat'")

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)
 
    # generate a list of asarray
    list_matrices = []
    for matrix in options.list_matrices:
        reduced_matrix = aims.read(matrix)
        reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]
        list_matrices.append(numpy.transpose(reduced_matrix))

    # two cases: 
    # (1) concatenated matrix
    # (2) averaged matrix
    if options.study == 'concat':
        concatenated_matrix = numpy.concatenate(list_matrices)
        options.matrix = matrix_converter(concatenated_matrix)
    else:
        sum_matrix = [sum(i) for i in zip(*list_matrices)]
        averaged_matrix = numpy.array(sum_matrix) / len(list_matrices)
        options.matrix = matrix_converter(averaged_matrix)

if __name__ == "__main__":
    main()
