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
* implementation of the kmedoids clustering algorithm

Main dependencies: numpy

Author: Sandrine Lefranc, 2015
"""

#----------------------------Imports-------------------------------------------


# python modules
import numpy


#----------------------------Functions-----------------------------------------


def kmedoids(distance_matrix, k, tmax=100):
    """
    """
    # determine dimensions of distance matrix
    m, n = distance_matrix.shape

    # randomly initialize an array of k medoid indices
    indices = numpy.sort(numpy.random.choice(n, k))

    # create a copy of the array of medoid indices
    new_indices = numpy.copy(indices)

    # initialize a dictionary to represent clusters
    clusters = {}

    for t in xrange(tmax):
        # determine clusters (arrays of data indices)
        J = numpy.argmin(distance_matrix[:, indices], axis=1)
        for kappa in range(k):
            clusters[kappa] = numpy.where(J == kappa)[0]

        # update cluster medoids
        for kappa in range(k):
            J = numpy.mean(distance_matrix[
                numpy.ix_(clusters[kappa], clusters[kappa])], axis=1)
            j = numpy.argmin(J)
            new_indices[kappa] = clusters[kappa][j]
        numpy.sort(new_indices)

        # check for convergence
        if numpy.array_equal(indices, new_indices):
            break

        indices = numpy.copy(new_indices)
    else:
        # final update of cluster memberships
        J = numpy.argmin(distance_matrix[:, indices], axis=1)
        for kappa in range(k):
            clusters[kappa] = numpy.where(J == kappa)[0]

    return indices, clusters
