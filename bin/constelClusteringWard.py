#! /usr/bin/env python
#############################################################################
# This software and supporting documentation are distributed by
#      CEA/NeuroSpin, Batiment 145,
#      91191 Gif-sur-Yvette cedex
#      France
# This software is governed by the CeCILL license version 2 under
# French law and abiding by the rules of distribution of free software.
# You can  use, modify and/or redistribute the software under the
# terms of the CeCILL license version 2 as circulated by CEA, CNRS
# and INRIA at the following URL "http://www.cecill.info".
#############################################################################

# python system modules
import sys
import numpy as np
from optparse import OptionParser

# soma
from soma import aims

# constel modules
from constel.lib.clustering.clusterstools import ward_method
from constel.lib.texturetools import texture_time


def parseOpts(argv):
    description = "Clustering Ward: one texture per subject"

    parser = OptionParser(description)

    parser.add_option(
        '-k', '--kmax', dest='kmax', type='int', help='Number of clusters')
    parser.add_option('-l', '--label', dest='patch',
                      help='Region of interest\'s,' 
                           'between 1 to parcelRegionsNb')
    parser.add_option('-x', '--matrix', dest='cmat', 
                      help='input reduced connectivity matrix')
    parser.add_option(
        '-c', '--distmatrix', dest='dmat', help='input distance matrix file')
    parser.add_option('-m', '--avgmesh', dest='avg_mesh',
                      help='Input average mesh to the correct group')
    parser.add_option('-g', '--gyritex', dest='gyri_tex', action='append', 
                      help='List of gyri segmentation (inputs)')
    parser.add_option(
        '-o', '--outputdir', dest='odir', help='output directory')
    parser.add_option('-t', '--textime', dest='tex_time',
                      action='append', help='List of time texture (outputs)')

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)

    # Load the reduced matrix
    matrix = aims.read(options.cmat)
    cmat = np.asarray(matrix)[:, :, 0, 0]

    # Compute linkage
    clusterid = ward_method(
        options.dmat, cmat.shape[0], options.odir, options.kmax)

    min_index = 0
    avg_mesh = aims.read(options.avg_mesh)
    nb_vertices = avg_mesh.vertex().size()

    # generate all individual textures time
    for i in xrange(len(options.gyri_tex)):
        tex = aims.read(options.gyri_tex[i])
        vertices_patch = []
        for j in xrange(tex[0].nItem()):
            if tex[0][j] == int(options.patch):
                vertices_patch.append(j)
        nb_vertex = len(vertices_patch)
        max_index = min_index + nb_vertex
        clusters = texture_time(options.kmax, clusterid, vertices_patch, 
                                nb_vertices, 2, min_index, max_index)
        min_index += nb_vertex
        aims.write(clusters, str(options.tex_time[i]))

if __name__ == "__main__":
    main()
