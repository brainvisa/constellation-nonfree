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
* compute squared Euclidean distance matrices (SEDM)
* evaluate various ways of computing SEDM

Main dependencies: numpy, scipy

Author: Sandrine Lefranc, 2015
"""

#----------------------------Imports-------------------------------------------


# python modules
import numpy

# scipy module
import scipy


#----------------------------Functions-----------------------------------------


def compute_SEDM(matrix, method=0):
    """Compute squared Euclidean distance matrices.

    Parameters
    ----------
    matrix:
    method: int
        method 0: a naïve approach
        method 1: avoiding square roots
        method 2: avoiding repeated inner products
        method 3: avoiding for loops
        method 4: resorting to build-in functions
        method 5: using inline C code

    Return
    ------
    dist_mat:
    """
    if method == 0:
        # determine dimensions of data matrix
        m, n = matrix.shape
        # initialize SEDM
        dist_mat = numpy.zeros((n, n))
        # iterate over upper triangle of dist_mat
        for i in range(n):
            for j in range(i + 1, n):
                dist_mat[i, j] = numpy.linalg(matrix[:, i] - matrix[:, j])**2
                dist_mat[j, i] = dist_mat[i, j]
    elif method == 1:
        # determine dimensions of data matrix
        m, n = matrix.shape
        # initialize SEDM
        dist_mat = numpy.zeros((n, n))
        # iterate over upper triangle of dist_mat
        for i in range(n):
            for j in range(i + 1, n):
                d = matrix[:, i] - matrix[:, j]
                dist_mat[i, j] = numpy.dot(d, d)
                dist_mat[j, i] = dist_mat[i, j]
    elif method == 2:
        # determine dimensions of data matrix
        m, n = matrix.shape
        # compute Gram matrix
        G = numpy.dot(matrix.T, matrix)
        # initialize SEDM
        dist_mat = numpy.zeros((n, n))
        # only iterate over upper triangle
        for i in range(n):
            for j in range(i + 1, n):
                # make use of |a-b|^2 = a'a + b'b -2a'b
                dist_mat[i, j] = G[i, i] - 2*G[i, j] + G[j, j]
                dist_mat[j, i] = dist_mat[i, j]
    elif method == 3:
        # determine dimensions of data matrix
        m, n = matrix.shape
        # compute Gram matrix
        G = numpy.dot(matrix.T, matrix)
        # compute matrix H
        H = numpy.tile(numpy.diag(G), (n, 1))
        dist_mat = H + H.T - 2*G
    elif method == 4:
        V = scipy.spatial.distance.pdist(matrix.T, "sqeuclidean")
        dist_mat = scipy.spatial.distance.squareform(V)
    elif method == 5:
        # determine dimensions of data matrix
        m, n = matrix.shape
        # compute Gram matrix
        G = numpy.dot(matrix.T, matrix)
        # initialize SEDM
        dist_mat = numpy.zeros((n, n))
        code = \
            """
            int i, j;
            for(i=0; i<n; i++){
                for(j=i+1; j<n; j++){
                    dist_mat(i, j) = G(i, i) - 2*G(i, j) + G(i, j);
                    dist_mat(j, i) = dist_mat(i, j)
            """
        dist_mat = scipy.weave.inline(
            code, ["G", "dist_mat", "n"],
            type_converters=scipy.weave.converters.blitz, compiler="gcc")
    return dist_mat

if __name__ == "__main__":
    import timeit
    n = 8
    m = [50, 100, 200, 400, 800, 1600, 3200, 6400]
    for i in m:
        matrix = numpy.matrix((n, i))
        print(timeit.timeit("compute_SEDM(matrix, method=0)",
                            setup="from __main__ import test"))
        print(timeit.timeit("compute_SEDM(matrix, method=1)",
                            setup="from __main__ import test"))
        print(timeit.timeit("compute_SEDM(matrix, method=2)",
                            setup="from __main__ import test"))
        print(timeit.timeit("compute_SEDM(matrix, method=3)",
                            setup="from __main__ import test"))
        print(timeit.timeit("compute_SEDM(matrix, method=4)",
                            setup="from __main__ import test"))
        print(timeit.timeit("compute_SEDM(matrix, method=5)",
                            setup="from __main__ import test"))