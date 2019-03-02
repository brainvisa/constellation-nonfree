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
Matrix utilities.
"""

# ---------------------------Imports-------------------------------------------

# System module
import csv
import numpy
import itertools

# Soma module
from soma import aims

# ---------------------------Functions-----------------------------------------


def calculate_percentage(matrix):
    """Calculate a percentage on the rows of a matrix.

    Parameters
    ----------
        matrix: numpy.array (mandatory)
            The connectivity matrix.

    Returns
    -------
        percent_matrix: numpy.array
            The table with the percentages of each line of the matrix.
    """
    # tuple of array dimensions
    ndmatrix = matrix.shape

    # initialyze a matrix with the same size as input matrix
    percent_matrix = numpy.zeros(
        (ndmatrix[0], ndmatrix[1]), dtype=numpy.float32)

    # work on the rows of the matrix
    for i in range(ndmatrix[0]):
        # each element of the rows is divided by total number
        for idx, value in enumerate(matrix[i]):
            percent = (value / sum(matrix[i])) * 100
            percent_matrix[i][idx] = percent

    return percent_matrix


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


def order_data_matrix(mat, labels):
    """
    Parameters
    ----------
        mat:
            matrix of observations, shape = (n,p), n observations, p attributes
        labels: numpy array
            of shape (n) corresponding to the labels of observations

    Returns
    -------
        order_mat: ordered matrix, of shape (n,p)
            labels: ordered labels
    """
    (n, p) = mat.shape
    obs_nb = labels.size
    if n != obs_nb:
        raise ValueError(
            "matrix dimensions and labels size are not compatible")
    labels_argsort = labels.argsort()
    order_mat = numpy.zeros((n, p), dtype=numpy.float32)
    for i in xrange(n):
        order_mat[i, :] = mat[labels_argsort[i], :]
    sortLabels = labels.copy()
    sortLabels.sort()
    return order_mat, sortLabels


def euclidian_distance(v1, v2):
    """
    Parameters
    ----------
        v1, v2: two vectors, shape (1,p) (numpy arrays)
    Returns
    -------
        dist(v1,v2)
    """
    dist = ((v1 - v2) ** 2).sum()
    return numpy.sqrt(dist)


def euclidian_distance_matrix(matrix):
    """
    """
    (n, p) = matrix.shape
    euclidian_dist_matrix = numpy.zeros((n, n), dtype=numpy.float)
    for i in xrange(n):
        v1 = matrix[i]
        dist_value = 0
        euclidian_dist_matrix[i][i] = 0
        for j in xrange(0, i):
            v2 = matrix[j]
            dist_value = euclidian_distance(v1, v2)
            euclidian_dist_matrix[i][j] = dist_value
            euclidian_dist_matrix[j][i] = dist_value
    return euclidian_dist_matrix


def compute_mclusters_by_nbasins_matrix(reducedmatrix, clusters, mode="mean"):
    """
    Compute the mean connectivity profile of each clusters and create the
    associated matrix.

    Parameters
    ----------
      reducedmatrix: numpy.array (mandatory)
          matrix of size M(vertex_ROI, vertex_targets)
      clusters: numpy.array (mandatory)
          connectivity-based parcelation of a ROI, labelled by clusters
      mode:
          "mean": normalized by the number of parcels vertices
          "sum": not normalized

    Return
    ------
      matrix: size (parcels_nb, targets_nb)
    """
    # if need, transpose the matrix to always work with
    # M(vertex_ROI, vertex_targets)
    if reducedmatrix.shape[0] < reducedmatrix.shape[1]:
        reducedmatrix = reducedmatrix.T

    # delete the zeros of the list to keep only the ROI vertices
    vertex_ROI = [i for i in clusters if i != 0]

    # return an error if the number of vertices is different between the rows
    # of the reduced matrix and the number of ROI vertices
    if len(vertex_ROI) != reducedmatrix.shape[0]:
        raise ValueError("reduced matrix and clustering must be the same size")

    # give the various labels
    labels = numpy.unique(vertex_ROI)

    # initialyze a matrix M(labels, targets)
    clusters_matrix = numpy.zeros(
        (len(labels), reducedmatrix.shape[1]), dtype=numpy.float32)

    # construct the clusters matrix M(labels, targets)
    for i, label in enumerate(labels):
        # give the matrix corresponding to idx
        cluster_matrix = reducedmatrix[numpy.where(vertex_ROI == label)]
        # mean or sum the profiles (rows) of the matrix
        if mode == "mean":
            clusters_matrix[i] = cluster_matrix.mean(axis=0)
        elif mode == "sum":
            clusters_matrix[i] = cluster_matrix.sum(axis=0)
        else:
            raise("The mode must be 'mean' or 'sum'.")

    return clusters_matrix


def write_matrix2csv(matrix, csvfilename):
    """Write a CSV file.

    Parameters
    ----------
        matrix: numpy.array (mandatory)
            matrix of size M(labels, vertex_targets)
            labels defines the clusters of the parcellation
        csvfilename: str (mandatory)
            directory and name of the CSV file
    """
    # tuple of array dimensions
    ndmatrix = matrix.shape

    # initialyze a list to add titles in columns of the csv file
    fieldnames = []

    # define the titles in columns of the csv file
    for i in range(ndmatrix[1] + 1):
        if i == 0:
            fieldnames.append("Cluster")
        else:
            fieldnames.append("Basin " + str(i))

    c = ndmatrix[0]
    # open a CSV file to write the results
    with open(csvfilename, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        m = 0
        # complete the table with the values of the matrix (values in percent)
        while c < (ndmatrix[0] + 1) and c > 0:
            dict = {}
            count = 0
            for element in fieldnames:
                if count == 0:
                    dict[str(element)] = m + 1
                else:
                    dict[str(element)] = matrix[m][count - 1]
                count += 1
            writer.writerow(dict)
            c -= 1
            m += 1


def contingency_matrix(labels1, labels2):
    """
    """
    classes = list(set(labels1))
    n = len(classes)
    contingency_matrix = numpy.array([zip(
        labels1, labels2).count(x) for x in itertools.product(
        classes, repeat=2)]).reshape(n, n)
    # contingency_matrix = numpy.bincount(
    # n * (labels1 - 1) + (labels2 - 1), minlength = n * n).reshape(n, n)
    return contingency_matrix


def partialWhiten(features):
    """Normalize the data (sigma = 1) per category.

    Parameters
    ----------
        features

    Returns
    -------
        white
    """
    dist = features[:, 0:6]
    direc = features[:, 6:12]
    sdist = numpy.std(dist)
    sidrec = numpy.std(direc)
    direc = direc / sidrec
    dist = dist / sdist
    white = numpy.hstack((dist, direc))
    return white


def save_normalization(filename):
    """Save the normalization.

    Parameters
    ----------
    filename: str (mandatory)
    """
    matrix = aims.read(filename)
    mat = numpy.array(matrix)[:, :, 0, 0]
    rows = mat.shape[0]
    norm_rows = []
    for i in range(rows):
        line_i = mat[i]
        norm_i = numpy.linalg.norm(line_i)
        norm_rows.append(norm_i)
    name1 = filename.split('.')[0]
    name2 = filename.split('.')[1]
    norm_name = name1 + name2 + "_normalization"
    numpy.save(norm_name, norm_rows)


def replace_negative_values(matrix, value=0.0):
    """ Replace the negative values by 0.0 (by default).
    """
    mat = aims.read(matrix)
    mat.muteToDense()
    numpy.asarray(
        mat.denseMatrix())[numpy.asarray(mat.denseMatrix()) < 0.0] = 0.0
    aims.write(mat, matrix)
