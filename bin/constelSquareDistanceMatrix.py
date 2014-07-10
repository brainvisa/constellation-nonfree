#!/usr/bin/env python

# python system modules
import numpy as np
import optparse
import sys

# scipy
from scipy.spatial.distance import pdist

def parseOpts(argv):
    desc = """Generate a square distance matrix in parallel."""
    parser = optparse.OptionParser(desc)
    parser.add_option('-m', '--matrix', dest='matrix',
                      help='input reduced connectivity matrix')   
    parser.add_option('-s', '--startiter', dest='starting_iteration', 
                      type='int', help='starting of iteration')                   
    parser.add_option('-t', '--nbiter', dest='nb_iteration', type='int',
                      help='number of iteration by job')                    
    parser.add_option('-k', '--idk', dest='id_distmat', type='int',
                      help='position in the distance matrix')
    parser.add_option('-i', '--input', dest='indir', metavar='FILE',
                      help='file to construct the distance matrix')
    parser.add_option('-d', '--method', dest='method', default='scipy',
                      help='three methods to generate the distance matrix:'
                          '(1) scipy --> risk of Memory Error,' 
                          '(2) square --> partition of square distance matrix,'
                          '(3) vector --> partition of distance matrix as vector'
                          '[default: %default]')
    return parser, parser.parse_args(argv)

def main():
    parser, (options, args) = parseOpts(sys.argv)
    
    # matrix loading
    matrix = np.load(options.matrix)
    m, n = matrix.shape

    #######################################################################
    #        Three methods to generate the distance matrix
    #          (1) scipy: generate all matrix, risk of Memory Error 
    #          (2) square: generate a portion of square distance matrix
    #          (3) vector: generate a portion of distance matrix as vector
    #######################################################################

    if options.method == 'scipy':
        # generate a triangular distance matrix as vector, with scipy
        # computes the Euclidean distance
        # WARNING: risk of Memory Error
        distance_matrix = pdist(matrix, metric='euclidean')
    elif options.method == 'square':
        # generate the distance matrix bounded by
        # starting_iteration and  max iteration required
        square_distance_matrix = np.zeros((options.nb_iteration,m))
        ii = 0
        for i in xrange(options.starting_iteration, options.nb_iteration + 
                    options.starting_iteration):
            for j in xrange(m):
                square_distance_matrix[ii,j] = np.sqrt(
                                       sum((matrix[i, :] - matrix[j, :]) * 
                                       (matrix[i, :] - matrix[j, :])))
            ii += 1
    else:
        # generate a triangular distance matrix as vector
        # bounded by starting_iteration and  max iteration required
        # advantage over the scipy: noy Memory Error
        size = 0
        count = 0
        for i in xrange(options.starting_iteration, options.nb_iteration + 
                    options.starting_iteration):
            size += m - i - 1  
        square_distance_matrix = np.zeros((size,))
        for i in xrange(options.starting_iteration, options.nb_iteration + 
                    options.starting_iteration):
            for j in xrange(i + 1, m):
                square_distance_matrix[count] = (np.sqrt(
                                       sum((matrix[i, :] - matrix[j, :]) * 
                                       (matrix[i, :] - matrix[j, :]))))
                count +=1
                
    # open the file to construct the square distance matrix    
    f = open(options.indir, 'r+')
    
    #######################################################################
    # complete the square distance matrix bounded by starting_iteration and 
    # max iteration required (It's not that big a deal)
    #######################################################################
    if options.method == 'scipy':
        # initialization of the square distance matrix
        square_distance_matrix = np.zeros((m, m), 'd')
        
        # upper/right triangular matrix
        k = int(options.id_distmat)
        for i in xrange(options.starting_iteration, options.nb_iteration + 
                        options.starting_iteration):
            for j in xrange(i + 1, m):
                square_distance_matrix[i][j] = distance_matrix[k]
                k += 1
                
        # lower/left triangular matrix
        kk = 0
        for jj in xrange(m):
            for ii in xrange(jj, m):
                if jj == ii:
                    pass
                elif ii <= (options.starting_iteration - 1):
                    kk += 1
                else:
                    square_distance_matrix[ii][jj] = distance_matrix[kk]
                    kk += 1    
        del distance_matrix
        
    # initial position to complete the distance matrix in the file 
    if options.method == 'square':
        # square distance matrix
        f.seek(options.starting_iteration * m * 8)
    else:
        # distance matrix as vector
        k = m - options.starting_iteration
        small_matrix = k*(k-1)/2
        large_matrix = m*(m-1)/2
        f.seek( (large_matrix - small_matrix) * 8)
        print (large_matrix - small_matrix)
    
    # write the distance matrix bounded by starting_iteration and 
    # max iteration required
    if options.method == 'scipy':
        f.write(square_distance_matrix[options.starting_iteration:options.nb_iteration + 
                options.starting_iteration][:].ravel())
    else:
        f.write(square_distance_matrix.ravel())
    
    print square_distance_matrix.ravel()

    # remove the square distance matrix to free the memory
    del square_distance_matrix
    f.close()
                    
if __name__ == "__main__":
    main()