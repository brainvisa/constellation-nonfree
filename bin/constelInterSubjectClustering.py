#!/usr/bin/env python
from optparse import OptionParser
from soma import aims
import numpy as np
import Pycluster as pc
import constel.lib.clustering.K_optimization as CK
import constel.lib.texturetools as tt
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

  #Rclustering = options.code
  #kmax = None  # the number of clusters is chosen to maximize the average silhouette width of the clustering.
  kmin = 2
  #n = 0
  #lkmax = []
  avg_mesh = aims.read( options.avg_mesh )
  subjectsNb = len(options.clustering_time)

  #if options.study == 'Concatenate':
    #for oneMatrix in xrange(subjectsNb):
      #labels = aims.read( options.avg_gyri_texture[n] )
      ##patch_mesh = aims.SurfaceManip.meshExtract( avg_mesh, labels, int(options.patch) )[0]
      ##patch_area = aims.SurfaceManip.meshArea( patch_mesh )
      ##kmax = round( patch_area / ( 2*options.areaMin_threshold ) ) + 1
      #lkmax.append( kmax )
      #n += 1
  #else:
    #labels = aims.read( options.avg_gyri_texture[n] )
    ##patch_mesh = aims.SurfaceManip.meshExtract( avg_mesh, labels, int(options.patch) )[0]
    ##patch_area = aims.SurfaceManip.meshArea( patch_mesh )
    ##kmax = round( patch_area / ( 2*options.areaMin_threshold ) ) + 1
    #lkmax.append( kmax )

  #kmax = max( lkmax )
  kmax = options.kmax
  #print 'list of k: ', lkmax
  print 'kmax = ', kmax
  #print 'Patch area = ', patch_area
  #print 'Area min threshold = ', options.areaMin_threshold

  reduced_matrix = aims.read(options.avg_matrix)
  if options.study == 'Concatenate':
    reduced_matrix = np.asarray(reduced_matrix)[:,:,0,0]
  if options.study == 'Average':
    reduced_matrix = np.transpose(np.asarray(reduced_matrix)[:,:,0,0])
  print "Reduced matrix of size (Nsample, Ndim): M", reduced_matrix.shape
  
  distance = pc.distancematrix(reduced_matrix, dist='e')
  
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
  
  
  #cl_dict, korder = CK.clusteringResults( reduced_matrix, kmin, kmax, Rclustering )
  #print 'k order: ', korder
  #print 'print dico: ', cl_dict
  #k_opt = korder[0]
  #print 'k optimale = ', k_opt

  countProcessedVertex = 0
  #listofClusteringSilhouette = []
  listOfClusteringTime = []
  #listOfClusteringKopt = []
  n = 0

  for subject in xrange(subjectsNb):
    #vertex_filename = options.vertex_index[n]
    #split_ima_filename = vertex_filename.split("\\")
    #len_split = len(split_ima_filename)
    #if len_split >1:
      #vertex_filename = split_ima_filename[0]
      #for i in xrange(1,len_split):
        #vertex_filename += split_ima_filename[i]
    if options.study == 'Concatenate':
      tex = aims.read( options.avg_gyri_texture[n] )
    else:
      tex = aims.read( options.avg_gyri_texture[0] )
    vertices_patch = []
    for i in xrange( tex[0].nItem() ):
      if tex[0][i] == int(options.patch):
        vertices_patch.append(i)
    #print 'vertex_index:', vertices_patch
    #print options.patch
    #vertices_patch = numpy.loadtxt(vertex_filename)
    #vertices_patch = vertices_patch.tolist()
    subjectPatchVertex_nb = len(vertices_patch)
    all_subjects_labels_list_index_min = countProcessedVertex
    all_subjects_labels_list_index_max = countProcessedVertex + subjectPatchVertex_nb

    if options.study == 'Concatenate':
      #clustersTime_tex, clus_avg_width_Time_tex = CK.texturesCreationMultiSubjects(cl_dict, nb_vertices, vertices_patch, all_subjects_labels_list_index_min, all_subjects_labels_list_index_max)
      k = kmax-1
      clusters = tt.textureTime(k, clusterID, nb_vertices, vertices_patch, 2, all_subjects_labels_list_index_min, all_subjects_labels_list_index_max)
    if options.study == 'Average':
      #clustersTime_tex, clus_avg_width_Time_tex = CK.texturesCreation(cl_dict, nb_vertices,
#vertices_patch)
      k = kmax-1
      clusters = tt.textureTime(k, clusterID, nb_vertices, vertices_patch)

    countProcessedVertex += subjectPatchVertex_nb
    #kopt_tex = aims.TimeTexture_S16()
    #k_opt_Time = k_opt - kmin
    #kopt_tex[0].assign( clustersTime_tex[k_opt_Time] )

    texTime = aims.write(clusters, str( options.clustering_time[n] ) )
    #texSilhouette = aims.write(clus_avg_width_Time_tex, str( options.clustering_silhouette[n] ) )
    #texKopt = aims.write(kopt_tex, str( options.clustering_Kopt[n] ) )
    listOfClusteringTime.append( texTime )
    #listofClusteringSilhouette.append( texSilhouette )
    #listOfClusteringKopt.append( texKopt )
    n += 1
  # One specific parcellation for each subject (a direct matching between the parcels across the subjects)
  listOfClusteringTime = options.clustering_time
  #listofClusteringSilhouette = options.clustering_silhouette
  #listOfClusteringKopt = options.clustering_Kopt
  #script_filename = options.clustering_result
  #script_file = open(script_filename, "a")
  #script_file.write( "---\n" + time.strftime( '%d/%m/%y %H:%M',time.localtime() ) + "\n---\n" )
  #script_file.write("Clustering of patch " + str(options.patch) + " of area " +
#str(patch_area) + " mm2 \n")
  #script_file.write("kmax :"+ str(kmax) + "\n")
  #script_file.write("k : \n")
  #for k in korder:
    #script_file.write(str(k) + " clusters, avg.width : " + str(cl_dict[k]['avg.width']) + "\n")
  #script_file.write("\n")
  #script_file.close()

if __name__ == "__main__" : main()