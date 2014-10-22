#!/usr/bin/env python
"""
@author: sl236442
"""

import numpy
import matplotlib.pyplot as plt


def verts_to_bbox(verts):
    """Determine the min and max of each axis.
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
    meshes: (AimsTimeSurface_3)

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
    print t.shape

    # to obtain binary matrix
    for p in asarray_submesh:
        p -= (a[0])
        p = numpy.round(p)
        t[p[0], p[1], p[2]] = 1

    print t[numpy.where(t > 0)]

    fig, ax = plt.subplots()
    ax.imshow(numpy.max(t, axis=2), cmap=plt.cm.gray, interpolation='nearest')
    plt.show()