#!/usr/bin/env python

from brainvisa.processes import *
from soma import aims
import optparse
import numpy as np
import Pycluster as pc
import constel.lib.clustering.K_optimization as CK
import constel.lib.texturetools as tt

def parseOpts(argv):
  desc="""
  For each subject, clustering was performed on the reduced connectivity matrix.
  The classical kmedoids algorithm has been used, with the euclidean distance
  between profiles as dissimilarity measure.
  """

  parser = optparse.OptionParser(desc)
  parser.add_option( '-m', '--matrix', dest='reduced_matrix', help='reduced \
                     connectivity matrix' )
  parser.add_option( '-p', '--patch', dest='patch', type='int', help='patch' )
  parser.add_option( '-g', '--gyri-segmentation', dest='gyri_segmentation' )
  parser.add_option( '-w', '--white-mesh', dest='mesh' )
  #parser.add_option( '-a', '--area-tresh', dest='thresholded_area' )
  parser.add_option( '-a', '--kmax', dest='kmax', type='int' )
  #parser.add_option( '-k', '--cluster-kopt', dest='clustering_kopt' )
  parser.add_option( '-t', '--cluster-time', dest='clustering_time' )
  #parser.add_option( '-d', '--cluster-kmed', dest='clustering_k_medoids' )
  #parser.add_option( '-s', '--cluster-silhouette', dest='clustering_silhouette' )
  #parser.add_option( '-l', '--cluster-vertex', dest='clustering_vertex_silhouette' )
  parser.add_option( '-r', '--result-gyrus', dest='clustering_result_gyrus' )
  #parser.add_option( '-f', '--result-full', dest='clustering_result_full' )

  return parser, parser.parse_args(argv)

def main():
  parser, (options, args) = parseOpts(sys.argv)

  Rclustering_filename = os.path.join( os.path.dirname( mainPath ),\
                                       'R', 'clustering', 'kmedoids.R' )

  mesh = aims.read(options.mesh)
  gyri_seg = aims.read(options.gyri_segmentation)
  reduced_matrix = aims.read(options.reduced_matrix)

  reduced_matrix = np.transpose(np.asarray(reduced_matrix)[:,:,0,0])
  print "Reduced matrix of size (Nsample, Ndim): M", reduced_matrix.shape
  
  distance = pc.distancematrix(reduced_matrix, dist='e')
 
  nb_vertices = mesh.vertex().size()
  #mesh_patch = aims.SurfaceManip.meshExtract(mesh, gyri_seg, options.patch)[0]
  #area_patch = aims.SurfaceManip.meshArea(mesh_patch)
  #kmax = round(area_patch / (2 * int(options.thresholded_area))) + 1
  kmax = options.kmax
  print "Clustering for K = 2 :", kmax
  
  kmin = 2
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
 
#subject_PatchCl_dict,k_ordering=CK.clusteringResults(reduced_matrix,kmin,kmax,Rclustering_filename)

  vertices_patch = []
  for i in xrange( gyri_seg[0].nItem() ):
    if gyri_seg[0][i] == options.patch:
      vertices_patch.append(i)
  #vertices_patch = np.array(vertices_patch)
  #vertices_patch = vertices_patch.tolist()
  #print vertices_patch
  n=kmax-1
  clusters = tt.textureTime(n, clusterID, nb_vertices, vertices_patch)
  
  #clusters, clus_avg_width_Time_tex = CK.texturesCreation(
#subject_PatchCl_dict, nb_vertices, vertices_patch )
 
  #kmedoidsTime_tex,vertex_silhouette_width_Time_tex=CK.textureKmedoidsCreation(subject_PatchCl_dict
#,nb_vertices, vertices_patch )

  ##outputs:
  #aims.write( kmedoidsTime_tex, options.clustering_k_medoids )
  #aims.write( vertex_silhouette_width_Time_tex, options.clustering_vertex_silhouette )
  #aims.write( clus_avg_width_Time_tex, options.clustering_silhouette )
  aims.write( clusters, options.clustering_time )

  #CK.writeClusteringResultsTxt( subject_PatchCl_dict, k_ordering, options.clustering_result_gyrus )
  #clusteringKopt_tex = aims.TimeTexture_S16()
  #clusteringKopt_tex[0].reserve( clusters[0].nItem() )
  #kopt_index_in_timetexture = k_ordering[0]-kmin
  #for i in xrange( clusters[0].nItem() ):
    #clusteringKopt_tex[0].push_back(clusters[kopt_index_in_timetexture][i] )
  #aims.write( clusteringKopt_tex, options.clustering_kopt )

  #CK.addOneSubjectIntraClusteringResultsTxt( subject_PatchCl_dict, k_ordering,
#options.clustering_result_full, options.patch )

if __name__ == "__main__" : main()