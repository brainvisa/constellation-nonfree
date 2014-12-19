#!/usr/bin/env python

# python system modules
from optparse import OptionParser
import numpy
import sys

# soma
from soma import aims

# pycluster
from Pycluster import kmedoids

# constellation
from constel.lib.texturetools import texture_time

# scipy
from scipy.spatial.distance import pdist, squareform


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
    parser.add_option('-l', '--label', dest='patch',type='int',
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


def main():
    parser, (options, args) = parseOpts(sys.argv)
    
    #########################################################################
    #                        clustering algorithm                           #
    #########################################################################
    # load reduced matrix: M(vertices_patch, basins)
    reduced_matrix = aims.read(options.group_matrix)
    reduced_matrix = numpy.asarray(reduced_matrix)[:, :, 0, 0]

    # generate the squareform distance matrix
    distance_matrix = squareform(pdist(reduced_matrix, metric='euclidean'))
    
    # generate several array containing the number of the cluster to which 
    # each item was assigned
    # i.e., one array by given nclusters (number_clusters)
    item_number = []
    for number_clusters in range(2, options.kmax + 1):
        # implements the k-medoids clustering algorithm (min of two clusters)
        # **ref**: http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster.pdf
        clusterid, err, nfound = kmedoids(
            distance_matrix, nclusters=number_clusters, npass=100)
        # rename the clusters from 1 to kmax
        for i, item in enumerate(numpy.unique(clusterid)):
            clusterid[clusterid == item] = i + 1
        item_number.append(clusterid)
        
    #########################################################################
    #                        time texture of clusters                       #
    #########################################################################
    # load mesh
    mesh = aims.read(options.mesh)
    nb_vertices = mesh.vertex().size()

    min_index = 0
    print options.patch
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
                                    item_number, 
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
