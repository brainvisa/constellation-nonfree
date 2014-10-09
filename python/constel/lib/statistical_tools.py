#! /usr/bin/env python
"""
Smoothing technique: moving average.
"""
# Created on Thu Oct  9 10:30:34 2014
# @author: sl236442

# Matplotlib module
import matplotlib.pyplot as plt


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
    plt.ylabel('Coordinates on isomap axis')
    plt.xlabel('Subjects')
    plt.plot(x, y, marker1, xm, ym, marker2)
    return plt.show()


def moving_average(x, y, k):
    """ Calculate the moving average.

    Parameters
    ----------
    x:
        list of obscissa
    y:
        list of ordinate
    k:
        moving average of order k

    Return
    ------
    mavg:
        values of the moving average
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