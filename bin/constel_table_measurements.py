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
"""


#----------------------------Imports-------------------------------------------

# python system module
from __future__ import absolute_import
import sys
import csv
import numpy
import argparse
import textwrap

# aims module
from soma import aims

# constel module
from constel.lib.utils.meshtools import mesh_surface_area
from six.moves import range


#----------------------------Functions-----------------------------------------

def parse_args(argv):
    """Parses the given list of arguments."""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            -------------------------------------------------------------------
            """))

    parser.add_argument("clusters", type=str, help="")
    parser.add_argument("white_mesh", type=str, help="")
    parser.add_argument("kmax", type=int, help="")
    parser.add_argument("matrix", type=str, help="")
    parser.add_argument("stats_file", type=str, help="")

    return parser, parser.parse_args(argv)


def mean(data):
    return sum(data, 0.0) / len(data)


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)
    
    # mesh
    mesh = aims.read(args.white_mesh)
    vertices = numpy.asarray(mesh.vertex())
    polygons = numpy.asarray(mesh.polygon())
    ftable = {}
    for k in range(2, args.kmax + 1):
        measures_table = []
        # read the parcellation at k clusters
        parcellation = aims.read(args.clusters)
        clusters = parcellation[k - 2].arraydata()
        # labels in the parcellation
        labels = numpy.unique(clusters)
        #######################################################################
        # number of vertex
        #######################################################################
        measures = []
        for label in labels.nonzero()[0]:
            num_vert = len(numpy.where(clusters == label)[0])
            measures.append(num_vert)
        #######################################################################
        # surface
        #######################################################################
        clusters_areas = numpy.zeros((numpy.max(clusters) + 1, ))
        verts = vertices[polygons]
        for idx, polygon in enumerate(polygons):
            # indexing to define two vector arrays from triangle vertices
            a = verts[:, 0, :][idx] - verts[:, 1, :][idx]
            b = verts[:, 0, :][idx] - verts[:, 2, :][idx]
            # area of triangle in 3D = 1/2 * Euclidean norm of cross product
            triangle_area = numpy.linalg.norm(numpy.cross(a, b), 2).sum() / 2.
            l1 = clusters[polygon[0]]
            clusters_areas[l1] += triangle_area / 3
            l1 = clusters[polygon[1]]
            clusters_areas[l1] += triangle_area / 3
            l1 = clusters[polygon[2]]
            clusters_areas[l1] += triangle_area / 3
        clusters_areas = clusters_areas[1:] # delete the background
        #######################################################################
        # moyenne
        #######################################################################
        mat = aims.read(args.matrix)
        m = numpy.asarray(mat)[:, :, 0, 0]
        if m.shape[0] > m.shape[1]:
            m = numpy.transpose(numpy.asarray(mat)[:, :, 0, 0])
        smat = numpy.sum(mat, axis=1)
        
            
        #m = mean(data)
        #######################################################################
        # variance
        #######################################################################
        #variance = mean([(x - m) ** 2 for x in data])
        #######################################################################
        # ecart type
        #######################################################################
        #variance(data)**0.5
        for i in range(len(measures)):
            measures_table.append([measures[i], clusters_areas[i]])
        ftable[k - 1] = measures_table

    with open(args.stats_file, "wb") as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=",")
        row = []
        row.extend([["Number of vertex", "Surface", "Variance", "Ecart type"]])
        csv_writer.writerow(row)
        for lid in sorted(ftable.keys()):
            row = [lid]
            row.extend(["{0}".format(elem) for elem in ftable[lid]])
            csv_writer.writerow(row)

if __name__ == "__main__":
    main()

