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
* check if two lists are the same number of elements.
* check no empty list.

Main dependencies:

Author: Sandrine Lefranc, 2014
"""

#----------------------------Imports-------------------------------------------


from __future__ import print_function
# system module
import numpy


#----------------------------Function------------------------------------------


def frame(data):
    """
    
    Parameters
    ----------
        data: (mandatory)
    Returns
    -------
        new:
    """
    x = data.shape[1] + 2
    y = data.shape[0] + 2
    new = [x * [0]] * y
    new = numpy.array(new)
    
    for i in range(1, y - 1):
        for j in range(1, x - 1):
            new[i][j] = data[i - 1][j - 1]
            
    return new 


def dilatation(img):
    """Dilatation of the image.
    
    Parameters
    ----------
        img: (mandatory)

    Returns
    -------
        dilate
    """
    dilate = [0] * (img.shape[1] - 2) * (img.shape[0] - 2)
    h = 0
    for i in range(1, img.shape[0] - 1):
        for j in range(1, img.shape[1] - 1):
            # max value (3*3) for the point(x, y)
            dilate[h] = max(
                [img[i - 1][j - 1], img[i][j - 1], img[i + 1][j - 1],
                 img[i - 1][j], img[i][j], img[i + 1][j], img[i - 1][j + 1],
                 img[i][j + 1], img[i + 1][j + 1]]) 
            h += 1
    
    dilate = numpy.array(dilate)
    dilate = numpy.reshape(dilate,(img.shape[0] - 2, img.shape[1] - 2)) 
    return dilate


def erosion(img):
    """Erosion of the image.
    
    Parameters
    ----------
        img: (mandatory)

    Returns
    -------
        erode
    """
    erode = [0] * (img.shape[1] - 2) * (img.shape[0] - 2)
    h = 0
    for i in range(1, img.shape[0] - 1):
        for j in range(1,img.shape[1] - 1):
            erode[h]=min(
                [img[i - 1][j - 1], img[i][j - 1], img[i + 1][j - 1],
                 img[i - 1][j], img[i][j], img[i + 1][j], img[i - 1][j + 1],
                 img[i][j + 1], img[i + 1][j + 1]]) 
            h += 1
    
    erode = numpy.array(erode)
    erode = numpy.reshape(erode, (img.shape[0] - 2, img.shape[1] - 2))
    return erode


def verts_to_bbox(verts):
    """Determine the min and max of each axis.
    
    Parameters
    ----------
        verts: (mandatory)

    Returns
    -------
        min_x: float
        min_y: float
        min_z: float
        max_x: float
        max_y: float
        max_z: float
    """
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    zs = [v[2] for v in verts]
    min_x = min(xs)
    min_y = min(ys)
    min_z = min(zs)
    max_x = max(xs)
    max_y = max(ys)
    max_z = max(zs)
    return min_x, min_y, min_z, max_x, max_y, max_z


def transform_mesh_to_volume(white_meshes):
    """Apply the meshes on a binary grid.

    Parameters
    ----------
    meshes: (AimsTimeSurface_3_VOID)

    Return
    ------
    """
    # asarray of submeshes to determine min/max of xi, yi and zi
    asarray_submesh = numpy.concatenate(white_meshes, axis=0)
    min_x, min_y, min_z = numpy.min(asarray_submesh, axis=0)
    max_x, max_y, max_z = numpy.max(asarray_submesh, axis=0)

    # bounding box: [(max_x, max_y, max_z), (min_x, min_y, min_z)]
    # floor on min and ceil on max
    a_min = [numpy.floor(min_x), numpy.floor(min_y), numpy.floor(min_z)]
    a_max = [numpy.ceil(max_x), numpy.ceil(max_y), numpy.ceil(max_z)]
    a = [numpy.asarray(a_min), numpy.asarray(a_max)]

    # voxel size
    spacing = (1., 1., 1.)

    # amax - amin
    size = (a[1] - a[0]) + 1
    t = numpy.zeros((size), dtype=numpy.single)
    print(t.shape)

    # to obtain binary matrix
    for p in asarray_submesh:
        p -= (a[0])
        p = numpy.round(p)
        t[p[0], p[1], p[2]] = 1
    
    print(t[numpy.where(t > 0)])
    
    # 2D
    t = numpy.max(t, axis=2)  

    # add frame
    t_frame = frame(t)
    # erosion
    erose = erosion(t_frame)
    # dilatation
    dilate = dilatation(t_frame)
    # opening
    opening = dilatation(erose)
    # closing    
    closing = erosion(dilate)
    
    return t, t_frame, dilate, erose, opening, closing
