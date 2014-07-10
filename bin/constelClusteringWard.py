#!/usr/bin/env python

# python system modules
from optparse import OptionParser
import numpy as np
import system

# soma
from soma import aims

# constel
from constel.lib.clustering.clusterstools import ward_method
from constel.lib.texturetools import texture_time

def parseOpts(argv):
    description = 'Clustering Ward: one texture per subject'

    parser = OptionParser(description)

    parser.add_option('-k', '--kmax', dest='kmax', type='int',
                      help='Number of clusters')
    parser.add_option('-l', '--label', dest='patch', 
                      help='Region of interest\'s, between 1 to parcelRegionsNb' )
    parser.add_option('-x', '--matrix', dest='cmat',
                      help='input reduced connectivity matrix')
    parser.add_option('-c', '--distmatrix', dest = 'dmat', 
                      help = 'input distance matrix file')
    parser.add_option('-m', '--avgmesh', dest = 'avg_mesh', 
                      help = 'Input average mesh to the correct group')
    parser.add_option('-g', '--gyritex', dest='gyri_tex', action='append', 
                      help='List of gyri segmentation (inputs)')
    parser.add_option('-o', '--outputdir', dest='odir', 
                      help='output directory')
    parser.add_option('-t', '--textime', dest = 'tex_time', action='append',
                      help='List of time texture (outputs)')

    return parser, parser.parse_args(argv)

def main():
    parser, (options, args) = parseOpts(sys.argv)
    
    # Load the distance matrix
    matrix = aims.read(options.cmat)
    cmat = np.asarray(matrix)[:, :, 0, 0]
    n, m = cmat.shape
    dist_mat_file = np.load(options.dmat)
    
    # Compute linkage
    clusterid = ward_method(dist_mat_file, n, options.odir, options.kmax)
    
    count_vertex = 0
    avg_mesh = aims.read(options.avg_mesh)
    nb_vertices = avg_mesh.vertex().size()
    nb_tex = len(options.gyri_tex)

    # generate all individual textures time 
    num = 0
    for subject in xrange(nb_tex):
        tex = aims.read(options.gyri_tex[num])
        vertices_patch = []
        for i in xrange(tex[0].nItem()):
            if tex[0][i] == int(options.patch):
                vertices_patch.append(i)
        nb_vertex = len(vertices_patch)
        min_index = count_vertex
        max_index = count_vertex + nb_vertex
        clusters = texture_time(options.kmax, 
                                clusterid, 
                                vertices_patch, 
                                nb_vertices, 
                                2,
                                min_index, 
                                max_index)
        count_vertex += nb_vertex
        aims.write(clusters, str(options.tex_time[num]))
        num += 1

if __name__ == "__main__" : main()