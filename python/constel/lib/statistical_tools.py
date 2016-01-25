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

Author: Sandrine Lefranc, 2014
"""

#----------------------------Imports-------------------------------------------


# system module
import numpy

# Matplotlib module
import matplotlib.pyplot as plt


#----------------------------Functions-----------------------------------------


def calculate_percentage(matrix):
    """Calculate a percentage on the rows of a matrix.

    Parameters
    ----------
        matrix: numpy.array (mandatory)

    Return
    ------
        percent_matrix: numpy.array
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


def sum_list(list_elements, start, end):
    """ Calculate the sum of the list elements from start/end list.

    Parameters
    ----------
    list_elements:
        give a list of elements
    start:
        first item in the list
    end:
        last item in the list

    Return
    ------
    s:
        sum of the list item
    """
    s = 0
    for k in xrange(start, end + 1):
        s += list_elements[k]
    return s


def draw_line(x, y, marker1, xm, ym, marker2):
    """ Draw a broken line passing through the points whose list of
    coordinates is given.

    Parameters
    ----------
    x:
        list of abscissa fir the input data
    y:
        list of ordinate for the input data
    marker1:
        marker and color of the input curve
    xm:
        list of abscissa for the moving average
    ym:
        list of ordinate for the movinf average
    marker2:
        marker and color of the moving average curve

    Return
    ------
    plot the curves on the same graphic
    """
    plt.title('Moving Average')
    plt.xlabel('Coordinates on isomap axis')
    plt.ylabel('Subjects')
    plt.plot(x, y, marker1, xm, ym, marker2)
    return plt.show()


def moving_average(x, y, k):
    """ Calculate the moving average.

    Parameters
    ----------
    x:
        list of abscissa
    y:
        list of ordinate
    k:
        moving average of order k

    Return
    ------
    xm, ym:
        values coordinates of the moving average
    """
    xm = []
    ym = []
    n = len(x)
    if k % 2 == 0:
        p = k // 2
        for j in xrange(p, n - p):
            ym += [(0.5 * y[j - p] + sum_list(y, j - p + 1, j + p - 1)
                   + 0.5 * y[j + p]) / k]
            xm += [x(j)]
    else:
        p = (k - 1) // 2
        for j in xrange(p, n - p):
            ym += [(sum_list(y, j - p, j + p)) / k]
            xm += [x[j]]
    return xm, ym


def overlap(x, y):
    i = 0
    j = 0
    c = 0
    len_x = len(x)
    len_y = len(y)
    while i < len_x and j < len_y:
        if x[i] > y[j]:
            j += 1
        elif x[i] < y[j]:
            i += 1
        else:  # x[i] == y[j]
            c += 1
            i += 1
            j += 1
    return c
