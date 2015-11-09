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

Main dependencies: PyAims library

Author: Sandrine Lefranc, 2015
"""

#----------------------------Imports-------------------------------------------


# python system module
import csv
import numpy
import itertools
import exceptions

#soma
from soma import aims


#----------------------------Functions-----------------------------------------


def resize_matrix(matrix, rows_length=100.0, cols_length=80.0):
    """Resize a matrix.

    The visual effect may be easier.

    Parameters
    ----------
        matrix: (aims.Volume)
        rows_length: fixed value
        cols_length: fixed value

    Return
    ------
        matrix: (aims.Volume)
    """
    if matrix.getSizeX() == matrix.getSizeY():
        matrix.header()['voxel_size'] = aims.vector_FLOAT(
            [rows_length / matrix.getSizeX(),
             rows_length / matrix.getSizeY(), 1, 1])
    else:
        matrix.header()['voxel_size'] = aims.vector_FLOAT(
            [rows_length / matrix.getSizeX(),
             cols_length / matrix.getSizeY(), 1, 1])

    return matrix


def permutation_resampling(features):
    """Resampling by permutation of the features

    Reference distribution will be sampled by permutations of the original
    one. The original sample represents the population.

    Parameters
    ----------
        features: matrix of size M(n_samples, n_dimension)
                  n_samples = features.shape[0]
                  n_dimension = features.shape[1]

    Returns
    -------
        new_features: new features matrix
    """
    n_samples, n_dimension = features.shape
    new_features = numpy.random.permutation(
        features[:, 0].reshape((n_samples, 1)))
    for f in range(1, n_dimension):
        new_features = numpy.hstack(
            (new_features, numpy.random.permutation(
                features[:, f].reshape((n_samples, 1)))))
    return new_features


def bootstrap_resampling(features):
    """Bootstrap resampling process

    The original sample represents the population from which it was drawn.
    So resamples from this sample represent what we would get if we took many
    samples from the population. The bootstrap distribution of a statistic,
    based on many resamples, represents the sampling distribution of the
    statistic, based on many samples.

    Parameters
    ----------
        features: discrete uniform matrix of size M(n_samples, n_dimension)
                  n_samples = features.shape[0]
                  n_dimension = features.shape[1]

    Returns
    -------
        new_features: new features matrix
        integers_vector: list of indices resampling
        survivors: list of items that survived the resampling
    """
    n_samples = features.shape[0]

    # Return random integers from the "discrete uniform" distribution in the
    # closed interval [low, high]: results are from [1, low]
    integers_vector = numpy.random.random_integers(n_samples, size=n_samples)
    integers_vector = integers_vector - 1

    # integers_vector becomes part of the features
    new_features = features[integers_vector[0], :]
    for i in range(1, integers_vector.size):
        new_features = numpy.vstack(
            (new_features, features[integers_vector[i], :]))

    # generate a list of items that survived the resampling
    survivors = numpy.array([])
    for j in integers_vector:
        if ((numpy.where(survivors == j)[0]).size == 0):
            survivors = numpy.append(survivors, numpy.array([j]))

    return new_features, integers_vector, survivors


def generate_uniform_matrix(features):
    """Generate samples from uniform distribution

    Samples with same shape than original samples.

    Parameters
    ----------
        features: matrix of size M(n_samples, n_dimension)
                  n_samples = features.shape[0]
                  n_dimension = features.shape[1]

    Return
    ------
        uniform_matrix: uniform matrix of size M(n_samples, n_dimension)
    """
    n_samples, n_dimension = features.shape
    vector_a = (features[:, 0].max() - features[:, 0].min()) * \
        numpy.random.random_sample((n_samples, 1)) + features[:, 0].min()
    vector_b = (features[:, 1].max() - features[:, 1].min()) * \
        numpy.random.random_sample((n_samples, 1)) + features[:, 1].min()
    uniform_matrix = numpy.hstack((vector_a, vector_b))
    for i in range(2, n_dimension):
        ni = (features[:, i].max() - features[:, i].min()) * \
            numpy.random.random_sample((n_samples, 1)) + features[:, i].min()
        uniform_matrix = numpy.hstack((uniform_matrix, ni))

    return uniform_matrix


def orderDataMatrix(mat, labels):
    """
    inputs:
            mat: matrix of observations, shape = (n,p), n observations, p attributes
            labels: numpy array of shape (n) corresponding to the labels of observations

    output:   order_mat: ordered matrix, of shape (n,p)
            labels: ordered labels
    """
    (n, p) = mat.shape
    obs_nb = labels.size
    if n != obs_nb:
        raise exceptions.ValueError(
            "matrix dimensions and labels size are not compatible")
    labels_argsort = labels.argsort()
    order_mat = numpy.zeros((n, p), dtype=numpy.float32)
    for i in xrange(n):
        order_mat[i, :] = mat[labels_argsort[i], :]
    sortLabels = labels.copy()
    sortLabels.sort()
    return order_mat, sortLabels


def euclidianDistance(v1, v2):
    """
    input:
          v1, v2: two vectors, shape (1,p) (numpy arrays)
    output:
          dist(v1,v2)
    """
    dist = ((v1 - v2) ** 2).sum()
    return numpy.sqrt(dist)


def euclidianDistanceMatrix(matrix):
    (n, p) = matrix.shape
    euclidian_dist_matrix = numpy.zeros((n, n), dtype=numpy.float)
    for i in xrange(n):
        v1 = matrix[i]
        dist_value = 0
        euclidian_dist_matrix[i][i] = 0
        for j in xrange(0, i):
            v2 = matrix[j]
            dist_value = euclidianDistance(v1, v2)
            euclidian_dist_matrix[i][j] = dist_value
            euclidian_dist_matrix[j][i] = dist_value
    return euclidian_dist_matrix


def compute_mclusters_by_nbasins_matrix(reducedmatrix, parcels,
                                        timestep=0, mode="meanOfProfiles"):
    """
    Compute the mean connectivity profile of each parcels and create the associated matrix:
    inputs:
      reducedmatrix: size (patchVertex_nb, targets_nb)
      vertices_patch: size (patchVertex_nb): index of the patch vertex in the parcels
      parcels: labeled aims Time Texture_S16, parcellation of the patch
      timestep: if the input parcels has several time steps, indicates the chosen step. (each step corresponding to a parcellation of the patch, time_step = 0, parcellation into 2 clusters)
      mode: "mean": normalized by the number of parcels vertices
         or "sum": not normalized
    outputs:
      matrix: size (parcels_nb, targets_nb)
    """
    (n, p) = reducedmatrix.shape
    print "reducedmatrix.shape", reducedmatrix.shape
    labels = numpy.zeros((n), dtype=int)

    vertices_patch = []

    for i in xrange(parcels[0].nItem()):
        if parcels[0][i] != 0:
            vertices_patch.append(i)

    vertices_patch = numpy.array(vertices_patch)
    vertices_patch = vertices_patch.tolist()
    for i in xrange(n):  # iteration on patch Vertex
        index_i = numpy.uint32(vertices_patch[i])
        labels[i] = parcels[timestep][index_i]

    labels_unique = numpy.unique(labels).tolist()

    matrix = numpy.zeros((len(labels_unique), p), dtype=numpy.float32)
    labelCount = 0
    for i in xrange(len(labels_unique)):
        label = labels_unique[i]
        label_connMatrixToTargets_array = reducedmatrix[
            numpy.where(labels == label)]
        if mode == "meanOfProfiles":
            matrix[labelCount] = label_connMatrixToTargets_array.mean(axis=0)
        elif mode == "sumOfProfiles":
            matrix[labelCount] = label_connMatrixToTargets_array.sum(axis=0)
        else:
            raise exceptions.ValueError(
                "the mode must be \"mean\" or \"sum\" ")
        labelCount += 1

    return matrix


def write_matrix2csv(matrix, csvfilename):
    """
    """
    fieldnames = []
    (n, p) = matrix.shape
    print n, p
    c = n
    for i in xrange(p+1):
        if i == 0:
            fieldnames.append("Cluster")
        else:
            fieldnames.append("Basin " + str(i))
    
    with open(csvfilename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        m = 0
        while c<(n+1) and c>0:
            dict = {}
            count = 0
            for element in fieldnames:
                if count == 0:
                    dict[str(element)] = m + 1
                else:
                    dict[str(element)] = matrix[m][count - 1]
                    print dict
                count += 1
            writer.writerow(dict)   
            c -= 1
            m += 1
        

def contingency_matrix(labels1, labels2):
    classes = list(set(labels1))
    n = len(classes)
    contingency_matrix = numpy.array([zip(
        labels1, labels2).count(x) for x in itertools.product(
        classes, repeat=2)]).reshape(n, n)
    #contingency_matrix = numpy.bincount(
        #n * (labels1 - 1) + (labels2 - 1), minlength = n * n).reshape(n, n)
    return contingency_matrix


def partialWhiten(features):
    """
    Normalize the data (sigma = 1) per category.
    """
    dist = features[:, 0:6]
    direc = features[:, 6:12]
    sdist = numpy.std(dist)
    sidrec = numpy.std(direc)
    direc = direc / sidrec
    dist = dist / sdist
    white = numpy.hstack((dist, direc))
    return white
