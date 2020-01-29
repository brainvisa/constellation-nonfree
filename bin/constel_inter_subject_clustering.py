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

Main dependencies: PyAims library
"""

#----------------------------Imports-------------------------------------------

# python system modules
from __future__ import print_function
from optparse import OptionParser
import pylab
import numpy
import sys

# soma
from soma import aims

# constellation
from constel.lib.utils.texturetools import texture_time

# scipy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster

# ---------------------------Functions-----------------------------------------


class MatplotlibFig(object):
    def __init__(self, fig):
        self._fig = fig
    def __del__(self):
        mainThreadActions().call(pylab.close, self._fig)


def validate( self ):
  import constel.lib.evaluation.indices as cv


def parseOpts(argv):
    description = """Connectivity-based parcellation of the patch.

    The joint connectivity matrix (concatenated or averaged) is computed.
    WARNING
    -------
        For concatenated approach, the number of subjects should be small
        to be easily used in this context.
        Alternatively, use "constelClusteringWard.py"

    The patch is clustered using the classical kmedoids algorithm applied on
    the euclidean distance matrix between the joint connectivity profiles.
    """

    parser = OptionParser(description)

    parser.add_option('-m', '--avgmesh', dest='mesh',
                      help='input: averaged mesh for a group of subjects')
    parser.add_option('-l', '--label', dest='patch', type='int',
                      help='region of interest number')
    parser.add_option('-t', '--gyri', dest='gyri_texture',
                      action='append',
                      help='input: list of gyri segmentation')
    parser.add_option('-a', '--kmax', dest='kmax', type='int')
    parser.add_option('-s', '--study', dest='study',
                      help='study : Average or Concatenate')
    parser.add_option('-g', '--avgmatrix', dest='group_matrix',
                      help='output: (averaged or concatenated) group matrix')
    parser.add_option('-p', '--time', dest='clustering_time',
                      action='append', help='list of textures')

    return parser, parser.parse_args(argv)


# ---------------------------Main program--------------------------------------


def main():
    parser, (options, args) = parseOpts(sys.argv)

    #########################################################################
    #                        clustering algorithm                           #
    #########################################################################
    # load reduced matrix: M(vertices_patch, basins)
    reduced_matrix = aims.read(options.group_matrix)
    reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]

    r = range(reduced_matrix.shape[0])
    numpy.random.shuffle(r)
    s = [r.index(i) for i in range(reduced_matrix.shape[0])]
    # Compute the distance matrix
    distmat = pdist(reduced_matrix[r], metric='euclidean')

    if options.study == 'avg':
        # this variant uses kmedoids from the pycluster module
        from Pycluster import kmedoids

        # generate the squareform distance matrix
        distance_matrix = squareform(distmat)

        # generate several array containing the number of the cluster to which
        # each item was assigned
        # i.e., one array by given nclusters (number_clusters)
        item_number = []
        import time
        numpy.random.seed(int(time.time() * 1000000) % 0x100000000)
        for number_clusters in range(2, options.kmax + 1):
            # implements the k-medoids clustering algorithm
            # (min of two clusters)
            # **ref**
            # http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster.pdf
            idx_medoids = numpy.random.randint(number_clusters, 
                                               size=reduced_matrix.shape[0])
            clusterid, err, nfound = kmedoids(
                distance_matrix, nclusters=number_clusters, npass=10)
            clusterid = clusterid[s]
            # rename the clusters from 1 to kmax
            for i, item in enumerate(numpy.unique(clusterid)):
                clusterid[clusterid == item] = i + 1
            item_number.append(clusterid)
    else:
        # this variant uses ward from the fastcluster module
        import fastcluster

        # compute linkage ward
        Z = fastcluster.linkage(distmat, method='ward', preserve_input=False)
        # labelization between 1 to nb
        clusterid = []
        for nb in range(1, options.kmax + 1):
            print("Trying {nb} cluster(s)".format(nb=nb))
            clusters = fcluster(Z, criterion='maxclust', t=nb)
            clusterid.append(clusters)


    #########################################################################
    #                        time texture of clusters                       #
    #########################################################################
    # load mesh
    mesh = aims.read(options.mesh)
    nb_vertices = mesh.vertex().size()

    min_index = 0
    for index, texture in enumerate(options.gyri_texture):
        # load a texture while identifying vertices of patch        
        tex = aims.read(texture)
        vertices_patch = numpy.where(tex[0].arraydata() == options.patch)[0]

        # two cases:
        # (1) concatenated approach: several textures of clusters (individual)
        # (2) averaged approach: a single texture of clusters
        if options.study == 'concat': 
            max_index = min_index + len(vertices_patch)
            clusters = texture_time(options.kmax, 
                                    clusterid, 
                                    vertices_patch, 
                                    nb_vertices, 
                                    2,
                                    min_index, 
                                    max_index)
            min_index += len(vertices_patch)
        else: #avg
            clusters = texture_time(
                options.kmax, item_number, vertices_patch, nb_vertices, 1)

        aims.write(clusters, str(options.clustering_time[index]))


if __name__ == "__main__":
    main()
