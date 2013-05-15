#!/usr/bin/env python
from optparse import OptionParser
from soma import aims
import numpy as N
import roca.lib.interMeshParcellation.processes_lib as T
import constel.lib.clustering.K_optimization as CK
import roca.lib.clustering.K_optimization as CCK
import sys
import time

def parseOpts(argv):
  description = 'Clustering of the subject : The gyrus vertices connectivity profiles of all the subjects are concatenated into a big matrix. The clustering is performed with the classical kmedoids algorithm and the Euclidean distance between profiles as dissimilarity measure...'

  parser = OptionParser(description)

  parser.add_option( '-r', '--codeR', dest = 'code', help = 'R clustering filename' )
  parser.add_option( '-m', '--avgmesh', dest = 'average_mesh', help = 'Input average mesh to the correct group' )
  parser.add_option( '-l', '--label', dest = 'patch_label', help = 'region of interest\'s, between 1 to parcelRegionsNb' )
  parser.add_option( '-t', '--avggyri', dest = 'average_gryi_segmentation', action='append', help = 'Input average gyri segmentation' ) 
  parser.add_option( '-a', '--areamin', dest = 'areaMin_threshold', type='int', help = 'minimum area for parcells in the resulting parcellation, in mm2' )
  parser.add_option( '-v', '--vertex', dest = 'vertex_index', action='append', help = 'index of vertex' )
  parser.add_option( '-s', '--study', dest = 'study', help = 'study : Average or Concatenate' )
  parser.add_option( '-g', '--avgmatrix', dest = 'avg_matrix', help = 'output average matrix' )
  parser.add_option( '-x', '--text', dest = 'clustering_result', help = 'output clustering result' )
  parser.add_option( '-o', '--silhouette', dest = 'clustering_silhouette', action='append', )
  parser.add_option( '-p', '--time', dest = 'clustering_time',  action='append' )
  parser.add_option( '-q', '--kopt', dest = 'clustering_Kopt', action='append' )

  return parser, parser.parse_args(argv)

def main():
  parser, (options, args) = parseOpts(sys.argv)

  Rclustering_filename = options.code
  kmax = None  # the number of clusters is chosen to maximize the average silhouette width of the clustering.
  kmin = 2
  n = 0
  if options.study == 'Concatenate':
    for oneMatrix in self.individual_reduced_matrix:
      patchsLabeled_tex = aims.read( options.average_gryi_segmentation[n] )
      n += 1
  else:
    patchsLabeled_tex = aims.read( options.average_gryi_segmentation[n] )

  average_mesh = aims.read( options.average_mesh )
  current_patch_mesh = aims.SurfaceManip.meshExtract( average_mesh, patchsLabeled_tex, int(options.patch_label) )[0]
  current_patch_area = aims.SurfaceManip.meshArea( current_patch_mesh )

  if kmax == None:
    #areaMin_threshold = 200
    kmax = round( current_patch_area / ( 2*options.areaMin_threshold ) ) + 1
    print 'kmax = ', kmax
    print 'current_patch_area: ', current_patch_area
    #print 'areaMin_threshold = ', areaMin_threshold
    print 'areaMin_threshold = ', options.areaMin_threshold

  if options.study == 'Concatenate':
    matrix = T.readMatrixImaAsNumpyArray( options.avg_matrix, True )
    all_subjects_PatchCl_dict, k_ordering_def = CK.clusteringsResults(matrix, kmin, kmax, Rclustering_filename )
  elif options.study == 'Average':
    matrix = T.readMatrixImaAsNumpyArray( options.avg_matrix, True )
    all_subjects_PatchCl_dict, k_ordering_def = CK.clusteringsResults(matrix, kmin, kmax, Rclustering_filename )
  k_opt = k_ordering_def[0]
  print 'k_opt = ', k_opt

  countProcessedVertex = 0
  WhiteMeshVertexNb = average_mesh.vertex().size()
  print 'WhiteMeshVertexNb : ', WhiteMeshVertexNb
  listofClusteringSilhouette = []
  listOfClusteringTime = []
  listOfClusteringKopt = []
  n = 0

  for subject in options.vertex_index:
    vertex_filename = options.vertex_index[n]
    split_ima_filename = vertex_filename.split("\\")
    len_split = len(split_ima_filename)
    if len_split >1:
      vertex_filename = split_ima_filename[0]
      for i in xrange(1,len_split):
        vertex_filename += split_ima_filename[i]
    subjectPatch_vertexIndex = N.loadtxt(vertex_filename)
    subjectPatch_vertexIndex = subjectPatch_vertexIndex.tolist()
    subjectPatchVertex_nb = len(subjectPatch_vertexIndex)
    all_subjects_labels_list_index_min = countProcessedVertex
    all_subjects_labels_list_index_max = countProcessedVertex + subjectPatchVertex_nb

    if options.study == 'Concatenate':
      clustersTime_tex, clus_avg_width_Time_tex = CCK.texturesCreationMultiSubjects(all_subjects_PatchCl_dict, WhiteMeshVertexNb, subjectPatch_vertexIndex, all_subjects_labels_list_index_min, all_subjects_labels_list_index_max)
    elif options.study == 'Average':
      clustersTime_tex, clus_avg_width_Time_tex = CK.texturesCreation(all_subjects_PatchCl_dict, WhiteMeshVertexNb, subjectPatch_vertexIndex)

    countProcessedVertex += subjectPatchVertex_nb
    kopt_tex = aims.TimeTexture_S16()
    k_opt_Time = k_opt - kmin
    print 'k_opt_Time = ', k_opt_Time
    kopt_tex[0].assign( clustersTime_tex[k_opt_Time] )

    texTime = aims.write(clustersTime_tex, str( options.clustering_time[n] ) )
    texSilhouette = aims.write(clus_avg_width_Time_tex, str( options.clustering_silhouette[n] ) )
    texKopt = aims.write(kopt_tex, str( options.clustering_Kopt[n] ) )
    listOfClusteringTime.append( texTime )
    listofClusteringSilhouette.append( texSilhouette )
    listOfClusteringKopt.append( texKopt )
    n += 1
  # One specific parcellation for each subject (a direct matching between the parcels across the subjects)
  listOfClusteringTime = options.clustering_time
  listofClusteringSilhouette = options.clustering_silhouette
  listOfClusteringKopt = options.clustering_Kopt
  script_filename = options.clustering_result
  script_file = open(script_filename, "a")
  script_file.write( "---\n" + time.strftime( '%d/%m/%y %H:%M',time.localtime() ) + "\n---\n" )
  script_file.write("Clustering of patch " + str(options.patch_label) + " of area " + str(current_patch_area) + " mm2 \n")
  script_file.write("kmax :"+ str(kmax) + "\n")
  script_file.write("k : \n")
  for k in k_ordering_def:
    script_file.write(str(k) + " clusters, avg.width : " + str(all_subjects_PatchCl_dict[k]['avg.width']) + "\n")
  script_file.write("\n")
  script_file.close()

if __name__ == "__main__" : main()