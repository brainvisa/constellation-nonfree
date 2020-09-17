#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################


# -------------------- imports ------------------------------------------------


# python system module
from __future__ import print_function
from __future__ import absolute_import
import argparse
import numpy
import sys

# aims module
from soma import aims

# Pycluster library
from Pycluster import kmedoids

# constel module
from constel.lib.utils.texturetools import texture_time

# SciPy library
from scipy.spatial.distance import pdist, squareform
from six.moves import range


# -------------------- definitions of functions -------------------------------


def parseArgs(argv):
    """Command-line options and arguments
    """
    parser = argparse.ArgumentParser(prog="constelIntraSubjectClustering",
        description="""For each subject, clustering was performed on the
        reduced connectivity matrix. The classical kmedoids algorithm has been
        used, with the euclidean distance between profiles as dissimilarity
        measure.""")
    parser.add_argument("matrix", type=str,
        help="a reduced connectivity matrix")
    parser.add_argument("patch", type=int,
        help="the patch number")
    parser.add_argument("gyri_segmentation", type=str,
        help="the gyri segmentation of the cortical surface")
    parser.add_argument("mesh", type=str,
        help="the mesh of the cortical surface")
    parser.add_argument("kmax", type=int,
        help="the number of clusters in the patterns")
    parser.add_argument("clustering_time", type=str,
        help="the result pattern")

    return parser, parser.parse_args(argv)


def rearrange_numbering(square_dmat, kmax, n_iter_k=1000):
    """Rearrange the numbering of clusters id.

    Parameters
    ----------
    square_dmat: array (mandatory)
        square-form distance matrix
    n_iter_k: int (optional, by default 1000)
        the number of times the k-medoids clustering algorithm is performed,
        each time with a different random initial condition
    kmax: int (mandatory)
        the number of clusters k

    Return
    ------
    clusterID: array
        an array containing the number of the cluster to which each item was
        assigned
    """
    clusterID = []
    for K in range(2, kmax + 1):
        clusterid, err, nfound = kmedoids(distance, K, NiterK)
        labels = 1
        for i in numpy.unique(clusterid):
            clusterid[clusterid == i] = labels
            labels += 1
        print('clusterid', numpy.unique(clusterid), 'for K =', K)
        clusterID.append(clusterid)
    return clusterID


def compute_distance_matrix(matrix, metric="euclidean", transpose=True):
    """Compute the square-form distance matrix.

    Parameters
    ----------
    matrix: ndarray (mandatory)
        an m by n array of m original observations in an n-dimensional space
    metric: str (optional, by default "euclidean")
        the metric can be "braycurtis", "canberra", "chebyshev", "cityblock",
        "correlation", "cosine", "dice", "euclidean", "hamming", "jaccard",
        "kulsinski", "mahalanobis", "matching", "minkowski", "rogerstanimoto",
        "russellrao", "seuclidean", "sokalmichener", "sokalsneath", "yule"
        "sqeuclidean".
    transpose: boolean (optional, by default True)
        the shape of the distance matrix must be (n_sample, n_dim)

    Return
    ------
    distance matrix: two-dimensional array
        the matrix containing the distances between the elements
        square-form distance matrix of shape (n_sample, n_sample)
    """
    if transpose:
        mat = numpy.transpose(numpy.asarray(matrix)[:,:,0,0])
    else:
        mat = numpy.asarray(matrix)[:,:,0,0]
    distance_matrix = squareform(pdist(mat, metric))
    return distance_matrix


# -------------------- main program -------------------------------------------


if __name__ == "__main__":
    parser, args = parseArgs(sys.argv)

    # load the files
    mesh = aims.read(args.mesh)
    gyri_seg = aims.read(args.gyri_segmentation)
    reduced_matrix = aims.read(args.matrix)

    # generate the distance matrix (euclidean distance)
    distance_matrix = compute_distance_matrix(reduced_matrix)

    # rearrange the numbering of the clusters
    clusterid = rearrange_numbering(distance_matrix, args.kmax)

    # extract the vertices of the patch
    vertices_patch = numpy.where(gyri_seg[0].arraydata() == args.patch)[0]

    # extract the total number of vertices of the mesh
    nb_vertices = mesh.vertex().size()

    # generate the clustering texture
    clusters = texture_time(kmax, clusterid, nb_vertices, vertices_patch, 1)

    # write the clustering result
    aims.write(clusters, args.clustering_time)
