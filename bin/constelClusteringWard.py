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
import numpy
from optparse import OptionParser

# soma
from soma import aims

# constel modules
from constel.lib.clustering.clusterstools import ward_method
from constel.lib.utils.texturetools import texture_time


def parseOpts(argv):
    description = """ Ward method minimizes the total within-cluster variance.
    At each step the pair of clusters with minimum cluster distance are merged.
    To implement this method, at each step find the pair of clusters that
    leads to minimum increase in total within-cluster variance after
    merging.

    Usage = %prog [options] filename
    """

    parser = OptionParser(description)

    parser.add_option(
        '-k', '--kmax', dest='kmax', type='int', help='Number of clusters')
    parser.add_option('-l', '--label', dest='patch', type='int',
                      help='Region of interest, between 1 to parcelRegionsNb')
    parser.add_option('-x', '--dendrogram', dest='dendrogram',
                      help='input reduced connectivity matrix')
    parser.add_option(
        '-c', '--distmatrix', dest='dmat', help='input distance matrix file')
    parser.add_option('-m', '--avgmesh', dest='avg_mesh',
                      help='Input average mesh to the correct group')
    parser.add_option('-g', '--gyritex', dest='gyri_tex', action='append',
                      help='List of gyri segmentation (inputs)')
#    parser.add_option(
#        '-o', '--outputdir', dest='odir', help='output directory')
    parser.add_option('-t', '--textime', dest='tex_time',
                      action='append', help='List of time texture (outputs)')

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)

#    # Load the reduced matrix
#    matrix = aims.read(options.group_matrix)
#    reduced_matrix = numpy.transpose(numpy.asarray(matrix)[:, :, 0, 0])
#
#    # Compute linkage
#    clusterid = ward_method(
#        options.dmat, reduced_matrix.shape[0], options.odir, options.kmax)

    # load dendrogram
    clusterid = numpy.load("/neurospin/archi-public/Users/lefranc/probabilistic27seeds/subjects/group_analysis/01to20-41to60/connectivity_clustering/concat/fsgroup/G25/smooth3.0/01to20-41to60_concat_fsgroup_G25_clusterid.npy")

    # load mesh
    avg_mesh = aims.read(options.avg_mesh)
    nb_vertices = avg_mesh.vertex().size()

    # generate all individual textures time
    min_index = 0
    for i in xrange(len(options.gyri_tex)):
        # load a texture while identifying vertices of patch
        tex = aims.read(options.gyri_tex[i])
        vertices_patch = numpy.where(tex[0].arraydata() == options.patch)[0]
        max_index = min_index + len(vertices_patch)
        print max_index
        clusters = texture_time(options.kmax, clusterid, vertices_patch,
                                nb_vertices, 2, min_index, max_index)
        min_index += len(vertices_patch)
        print min_index
        aims.write(clusters, str(options.tex_time[i]))

if __name__ == "__main__":
    main()
