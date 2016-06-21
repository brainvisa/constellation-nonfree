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
* calculate the percentage on the rows of a matrix.

Main dependencies:

Author: Sandrine Lefranc, 2014
"""

#----------------------------Imports-------------------------------------------


# system module
import numpy


#----------------------------Functions-----------------------------------------


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
