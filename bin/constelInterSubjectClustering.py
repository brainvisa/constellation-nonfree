#!/usr/bin/env python
from optparse import OptionParser
from soma import aims
import numpy as np
import Pycluster as pc
import constel.lib.clustering.K_optimization as CK
import constel.lib.texturetools as tt
import scipy.spatial.distance as ssd
import sys
import time

def parseOpts(argv):
  description = 'Clustering of the subject : The gyrus vertices connectivity profiles of all the subjects are concatenated into a big matrix. The clustering is performed with the classical kmedoids algorithm and the Euclidean distance between profiles as dissimilarity measure...'

  parser = OptionParser(description)

  #parser.add_option( '-r', '--codeR', dest = 'code', help = 'R clustering filename' )
  parser.add_option( '-m', '--avgmesh', dest = 'avg_mesh', help = 'Input average mesh to the correct group' )
  parser.add_option( '-l', '--label', dest = 'patch', help = 'region of interest\'s, between 1 to parcelRegionsNb' )
  parser.add_option( '-t', '--avggyri', dest = 'avg_gyri_texture', action='append', help = 'Input average gyri segmentation' )
  parser.add_option( '-a', '--kmax', dest='kmax', type='int' )
  #parser.add_option( '-a', '--areamin', dest = 'areaMin_threshold', type='int', help = 'minimum area for parcells in the resulting parcellation, in mm2' )
  parser.add_option( '-s', '--study', dest = 'study', help = 'study : Average or Concatenate' )
  parser.add_option( '-g', '--avgmatrix', dest = 'avg_matrix', help = 'output average matrix' )
  #parser.add_option( '-x', '--text', dest = 'clustering_result', help = 'output clustering result' )
  #parser.add_option( '-o', '--silhouette', dest = 'clustering_silhouette', action='append', )
  parser.add_option( '-p', '--time', dest = 'clustering_time',  action='append' )
  #parser.add_option( '-q', '--kopt', dest = 'clustering_Kopt', action='append' )

  return parser, parser.parse_args(argv)

def main():
  parser, (options, args) = parseOpts(sys.argv)

  kmin = 2
  avg_mesh = aims.read( options.avg_mesh )
  subjectsNb = len(options.clustering_time)

  kmax = options.kmax
  print 'kmax = ', kmax

  reduced_matrix = aims.read(options.avg_matrix)
  if options.study == 'Concatenate':
    reduced_matrix = np.asarray(reduced_matrix)[:,:,0,0]
  if options.study == 'Average':
    reduced_matrix = np.transpose(np.asarray(reduced_matrix)[:,:,0,0])
  print "Reduced matrix of size (Nsample, Ndim): M", reduced_matrix.shape
  
  distance = ssd.squareform(ssd.pdist(reduced_matrix, 'euclidean'))
  
  nb_vertices = avg_mesh.vertex().size()
  
  NiterK = 100
  clusterID = []
  for K in range(2, kmax+1):
    clusterid, err, nfound = pc.kmedoids(distance, K, NiterK)
    labels = 1
    for i in np.unique(clusterid):
      clusterid[clusterid == i] = labels
      labels += 1
    print 'clusterid', np.unique(clusterid), 'for K =', K
    clusterID.append(clusterid)
  
  countProcessedVertex = 0
  listOfClusteringTime = []
  n = 0

  for subject in xrange(subjectsNb):
    if options.study == 'Concatenate':
      tex = aims.read( options.avg_gyri_texture[n] )
    else:
      tex = aims.read( options.avg_gyri_texture[0] )
    vertices_patch = []
    for i in xrange( tex[0].nItem() ):
      if tex[0][i] == int(options.patch):
        vertices_patch.append(i)
    subjectPatchVertex_nb = len(vertices_patch)
    all_subjects_labels_list_index_min = countProcessedVertex
    all_subjects_labels_list_index_max = countProcessedVertex + subjectPatchVertex_nb

    if options.study == 'Concatenate':
      k = kmax-1
      clusters = tt.texture_time(k, clusterID, nb_vertices, vertices_patch, 2, all_subjects_labels_list_index_min, all_subjects_labels_list_index_max)
    if options.study == 'Average':
      k = kmax-1
      clusters = tt.texture_time(k, clusterID, nb_vertices, vertices_patch, 1)

    countProcessedVertex += subjectPatchVertex_nb

    texTime = aims.write(clusters, str( options.clustering_time[n] ) )
    listOfClusteringTime.append( texTime )
    n += 1
  # One specific parcellation for each subject (a direct matching between the parcels across the subjects)
  listOfClusteringTime = options.clustering_time

if __name__ == "__main__" : main()