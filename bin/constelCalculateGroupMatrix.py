#!/usr/bin/env python

# python system modules
import sys
import numpy
import optparse

# soma
from soma import aims

# constel
from constel.lib.connmatrix.connmatrixtools import resize_matrix

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
        list_matrices.append(reduced_matrix.astype('float32'))
    print list_matrices
        
    # two cases: 
    # (1) concatenated matrix: M(target_regions, group_vertices_patch)
    # (2) averaged matrix: M(target_regions, vertices_patch)
    if options.study == 'concat':
        concatenated_matrix = numpy.concatenate(list_matrices, axis=1)
        concatenated_matrix = resize_matrix(aims.Volume(concatenated_matrix))
        aims.write(concatenated_matrix, options.matrix)
    else: #avg
        sum_matrix = [sum(i) for i in zip(*list_matrices)]
        averaged_matrix = numpy.array(sum_matrix) / len(list_matrices)
        print averaged_matrix.shape
        averaged_matrix = resize_matrix(aims.Volume(averaged_matrix))
        aims.write(averaged_matrix, options.matrix)

 
#    if options.study == 'concat':
#        shapex = None
#        shapey = 0
#        for matrix in options.list_matrices:
#            reduced_matrix = aims.read(matrix)
#            reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]
#            if shapex is None:
#                shapex = reduced_matrix.shape[0]
#            shapey += reduced_matrix.shape[1]
#        # generate a list of asarray (data_type: DOUBLE)
#        l = numpy.zeros((shapex, shapey), dtype=numpy.float16)
#        shifty = 0
#        for matrix in options.list_matrices:
#            reduced_matrix = aims.read(matrix)
#            reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]
#            l[:, shifty:(shifty+reduced_matrix.shape[1])] = reduced_matrix
#            shifty += reduced_matrix.shape[1]
#
#    # two cases: 
#    # (1) concatenated matrix: M(target_regions, group_vertices_patch)
#    # (2) averaged matrix: M(target_regions, vertices_patch)
#    
##        concatenated_matrix = numpy.concatenate(list_matrices, axis=1)
#        concatenated_matrix = resize_matrix(aims.Volume(l, dtype='int16'))
#        aims.write(concatenated_matrix, options.matrix)
#    else: #avg
#        list_matrices = []
#        for matrix in options.list_matrices:
#            reduced_matrix = aims.read(matrix)
#            reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]
#            list_matrices.append(reduced_matrix)
#        sum_matrix = [sum(i) for i in zip(*list_matrices)]
#        averaged_matrix = numpy.array(sum_matrix) / len(list_matrices)
#        averaged_matrix = resize_matrix(aims.Volume(averaged_matrix))
#        aims.write(averaged_matrix, options.matrix)

if __name__ == "__main__":
    main()
