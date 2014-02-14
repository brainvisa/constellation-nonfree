import constel.lib.clustering.clustertools as cct
import constel.lib.texturetools as tt
from optparse import OptionParser
from soma import aims
import numpy as np
import sys

def parseOpts(argv):
    description = 'Clustering Ward: one texture per subject'

    parser = OptionParser(description)

    parser.add_option('-k', '--kmax', dest='kmax', type='int',
        help='Number of clusters')
    parser.add_option('-l', '--label', dest='patch', 
        help='Region of interest\'s, between 1 to parcelRegionsNb' )
    parser.add_option('-c', '--concatmatrix', dest = 'cmat', 
        help = 'output concatenated matrix')
    parser.add_option('-m', '--avgmesh', dest = 'avg_mesh', 
        help = 'Input average mesh to the correct group')
    parser.add_option('-g', '--gyritex', dest='gyri_tex', action='append', 
        help='List of gyri segmentation (inputs)')
    parser.add_option('-o', '--outputdir', dest='odir', 
        help='output directory')
    parser.add_option('-t', '--textime', dest = 'tex_time',  action='append',
        help='List of time texture (outputs)')

    return parser, parser.parse_args(argv)

def main():
    parser, (options, args) = parseOpts(sys.argv)
    
    # Load the matrix
    matrix = aims.read(self.group_matrix.fullPath())
    cmat = np.asarray(matrix)[:, :, 0, 0]
    
    # Compute linkage
    clusterid = cct.ward_method(cmat, options.odir)
    count_vertex = 0
    avg_mesh = aims.read(options.avg_mesh)
    nb_vertices = avg_mesh.vertex().size()
    nb_tex = len(options.gyri_tex)

    n = 0
    for subject in xrange(nb_tex):
        tex = aims.read(options.gyri_tex[n])
        vertices_patch = []
        for i in xrange(tex[0].nItem()):
        if tex[0][i] == int(options.patch):
            vertices_patch.append(i)
        nb_vertex = len(vertices_patch)
        min_index = count_vertex
        max_index = count_vertex + nb_vertex
        clusters = tt.texture_time(k, clusterid, nb_vertices, vertices_patch, 2,
            min_index, max_index)
        count_vertex += nb_vertex
        tex_time = aims.write(clusters, str(options.tex_time[n]))
        n += 1

if __name__ == "__main__" : main()