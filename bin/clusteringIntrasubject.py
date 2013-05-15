#!/usr/bin/env python

from brainvisa.processes import *
from soma import aims
import optparse
import numpy

import constel.lib.clustering.K_optimization as CCK
import roca.lib.clustering.K_optimization as CK

def parseOpts(argv):
  desc="""The connectivity profiles of a given region cortex are clustered using k-medoids approach."""

  parser = optparse.OptionParser(desc)
  parser.add_option( '-m', '--matrix', dest='connectivity_matrix_reduced', help='connectivity matrix reduced' )
  parser.add_option( '-p', '--patch', dest='patch_label', type='int', help='patch_label' )
  parser.add_option( '-g', '--gyri-segmentation', dest='gyri_segmentation' )
  parser.add_option( '-w', '--white-mesh', dest='white_mesh' )
  parser.add_option( '-a', '--area-tresh', dest='areaMin_threshold' )
  parser.add_option( '-v', '--vertex-index', dest='vertex_index' )
  parser.add_option( '-k', '--cluster-kopt', dest='clustering_kopt' )
  parser.add_option( '-t', '--cluster-time', dest='clustering_time' )
  parser.add_option( '-d', '--cluster-kmed', dest='clustering_k_medoids' )
  parser.add_option( '-s', '--cluster-silhouette', dest='clustering_silhouette' )
  parser.add_option( '-l', '--cluster-vertex', dest='clustering_vertex_silhouette' )
  parser.add_option( '-r', '--result-gyrus', dest='clustering_result_gyrus' )
  parser.add_option( '-f', '--result-full', dest='clustering_result_full' )

  return parser, parser.parse_args(argv)

def main():
  parser, (options, args) = parseOpts(sys.argv)

  patchsLabeled_tex = aims.read( options.gyri_segmentation )
  Rclustering_filename = os.path.join( os.path.dirname( mainPath ), 'R', 'clustering', 'kmedoids.R' )
  subject_reducedConnMatrix = numpy.asarray( aims.read(
    options.connectivity_matrix_reduced ) )
  subject_reducedConnMatrix = subject_reducedConnMatrix.reshape(
    subject_reducedConnMatrix.shape[0], subject_reducedConnMatrix.shape[1] )

  subject_white_mesh = aims.read( options.white_mesh )
  subjectWhiteMesh_vertexNb = subject_white_mesh.vertex().size()
  current_patch_mesh = aims.SurfaceManip.meshExtract( subject_white_mesh, patchsLabeled_tex, options.patch_label )[0]
  current_patch_area = aims.SurfaceManip.meshArea( current_patch_mesh )

  kmin = 2
  kmax = round( current_patch_area / ( 2*int(options.areaMin_threshold) ) ) + 1

  subject_PatchCl_dict, k_ordering = CCK.clusteringResults( subject_reducedConnMatrix, kmin, kmax, Rclustering_filename )

  subjectPatch_vertexIndex =  numpy.loadtxt( options.vertex_index )
  subjectPatch_vertexIndex = subjectPatch_vertexIndex.tolist()
  clustersTime_tex, clus_avg_width_Time_tex = CK.texturesCreation( subject_PatchCl_dict, subjectWhiteMesh_vertexNb, subjectPatch_vertexIndex )
  kmedoidsTime_tex, vertex_silhouette_width_Time_tex = CK.textureKmedoidsCreation( subject_PatchCl_dict,subjectWhiteMesh_vertexNb, subjectPatch_vertexIndex )

  #outputs:
  aims.write( kmedoidsTime_tex, options.clustering_k_medoids )
  aims.write( vertex_silhouette_width_Time_tex, options.clustering_vertex_silhouette )
  aims.write( clus_avg_width_Time_tex, options.clustering_silhouette )
  aims.write( clustersTime_tex, options.clustering_time )

  CK.writeClusteringResultsTxt( subject_PatchCl_dict, k_ordering, options.clustering_result_gyrus )
  clusteringKopt_tex = aims.TimeTexture_S16()
  clusteringKopt_tex[0].reserve( clustersTime_tex[0].nItem() )
  kopt_index_in_timetexture = k_ordering[0]-kmin
  for i in xrange( clustersTime_tex[0].nItem() ):
    clusteringKopt_tex[0].push_back( clustersTime_tex[kopt_index_in_timetexture][i] )
  aims.write( clusteringKopt_tex, options.clustering_kopt )

  CK.addOneSubjectIntraClusteringResultsTxt( subject_PatchCl_dict, k_ordering, options.clustering_result_full, options.patch_label )

if __name__ == "__main__" : main()