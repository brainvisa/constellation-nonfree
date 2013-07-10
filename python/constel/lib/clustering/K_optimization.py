
import constel.lib.clustering.kmedoids as rkmedoids
import numpy as N
import numpy.random as NR
from soma import aims
import decimal as pydeci
import time
import constel



def clusteringResults(X, k_min, k_max, Rclustering_filename, diss=False,
    init_kmedoids=None):
  """
  inputs:
            X: data = n observations to cluster (shape = (n,p) )
            diss: dissimilarity matrix: if true : X is the dissimilarity matrix otherwise, X is a data matrix: each row corresponds to an observation, and each column corresponds to a variable.
            init_medoids: numpy array of length k of integer indices (in 1:n) specifying initial medoids instead of using the build algorithm
  """
  #k_min = 2
  #Rclustering_filename = "kmedoids.R"
  #X = NR.randn(10,2)
  #A = N.concatenate([N.ones((7,2)),N.zeros((3,2))])
  #X = X+3*A;
  #k_max = 5

  k = k_min
  cl_dict = {}
  avg_width_ar = N.zeros((k_max - k_min + 1), dtype = N.float)

  while(k <= k_max):
    print "X.shape:", X.shape
    clusteringMethod = rkmedoids.KMEDOIDS_R(Rclustering_filename)
    cl = clusteringMethod.predict(X, k, 0, diss, init_kmedoids)
    #print "save clustering results:"
    cl_dict[k] = {}
    cl_dict[k]['avg.width'] = cl['silinfo']['avg.width']
    #print "save clus.avg.widths..."
    cl_dict[k]['clus.avg.widths'] = cl['silinfo']['clus.avg.widths']
    #print "save clustering...", N.unique((N.asarray(cl["clustering"])).copy()).tolist()
    cl_dict[k]['labels_list'] = (N.asarray(cl["clustering"])).copy().tolist()
    #print "save medoids..."
    cl_dict[k]['medoids'] = N.asarray(cl["id.med"]) - 1
    avg_width_ar[k - k_min] = cl['silinfo']['avg.width']
    cl_dict[k]['widths']=cl['silinfo']['widths'][:,2]
    k += 1
  k_ordering = N.argsort(avg_width_ar)
  k_ordering = k_ordering + k_min
  i = k_ordering.size - 1
  k_ordering_def = N.zeros((k_ordering.size), dtype = N.int)
  for i in xrange(k_ordering.size):
    k_ordering_def[i] = k_ordering[k_ordering.size - i - 1 ]
  print "k in order of validity:", k_ordering_def
  return cl_dict, k_ordering_def


def texturesCreation(clusterings_dict, meshVertex_nb, seedVertex_index):
  """
  inputs:
            clusterings_dict: result of clusteringResults(X, k_min, k_max, Rclustering_filename) method
            meshVertex_nb: vertex number for output textures
            seedVertex_index: vertex index of the seed region (observation in the matrix X

  outputs:
            clustersTime_tex: for each time step: results of the k th clustering
            clus_avg_width_Time_tex: for each time step, clus.avg.width associated to each cluster
  """
  stepCount = 0
  clustersTime_tex = aims.TimeTexture_S16()
  clus_avg_width_Time_tex = aims.TimeTexture_FLOAT()
  for k in clusterings_dict.keys():
    clustersTime_tex[stepCount].resize( meshVertex_nb, 0 )
    clustersTime_tex[stepCount].arraydata()[ seedVertex_index ] \
      = clusterings_dict[k]['labels_list']
    current_clust_avg_width_tex = constel.oneTargetDensityTargetsRegroupTexture(N.asarray(clusterings_dict[k]['clus.avg.widths']), clustersTime_tex, stepCount)
    clus_avg_width_Time_tex[stepCount].assign( current_clust_avg_width_tex[0] )
    stepCount += 1

  return clustersTime_tex, clus_avg_width_Time_tex


def texturesCreationMultiSubjects( clusterings_dict, meshVertex_nb, 
                                   seedVertex_index, labels_list_min_index, 
                                   labels_list_max_index ):
  """
  inputs:
            clusterings_dict: result of clusteringResults(X, k_min, k_max, Rclustering_filename) method
            meshVertex_nb: vertex number for output textures
            seedVertex_index: vertex index of the seed region (observation in the matrix X
            labels_list_min_index
            labels_list_max_index

  outputs:
            clustersTime_tex: for each time step: results of the k th clustering
            clus_avg_width_Time_tex: for each time step, clus.avg.width associated to each cluster
  """
  stepCount = 0
  clustersTime_tex = aims.TimeTexture_S16()
  clus_avg_width_Time_tex = aims.TimeTexture_FLOAT()
  for k in clusterings_dict.keys():
    clustersTime_tex[stepCount].resize( meshVertex_nb, 0 )
    clustersTime_tex[stepCount].arraydata()[ seedVertex_index ] \
      = clusterings_dict[k]['labels_list'][
        labels_list_min_index:labels_list_max_index]
    current_clust_avg_width_tex \
      = constel.oneTargetDensityTargetsRegroupTexture(
        N.asarray(clusterings_dict[k]['clus.avg.widths']), clustersTime_tex, 
        stepCount )
    clus_avg_width_Time_tex[stepCount].assign( current_clust_avg_width_tex[0] )
    stepCount += 1

  return clustersTime_tex, clus_avg_width_Time_tex


def textureKmedoidsCreation(clusterings_dict, meshVertex_nb, seedVertex_index):
  """
  Creation of an aims Time texture_S16
  for each time step: labeled texture of kmedoids of the associated clustering
  
  imputs:
            clusterings_dict: result of clusteringResults(X, k_min, k_max, Rclustering_filename) method
            meshVertex_nb: vertex number for output textures
            seedVertex_index: vertex index of the seed region (observation in the matrix X
            
  outputs:
            kmedoidsTime_tex: for each time step: k mdedoids of the k th clustering
            vertex_silhouette_width_Time_tex: for each time step, width associated to each vertex
  """
  stepCount = 0
  kmedoidsTime_tex = aims.TimeTexture_S16()
  vertex_silhouette_width_Time_tex = aims.TimeTexture_FLOAT()
  for k in clusterings_dict.keys():
    kmedoidsTime_tex[stepCount].resize( meshVertex_nb, 0 )
    vertex_silhouette_width_Time_tex[stepCount].resize( meshVertex_nb, 0 )
    vertsilarr = vertex_silhouette_width_Time_tex[stepCount].arraydata()
    vertsilarr[ seedVertex_index ] = clusterings_dict[k]['widths'].astype(
      vertsilarr.dtype )
    kmedoids_array = clusterings_dict[k]['medoids']
    for kmed in kmedoids_array:
      kmedIndexInWhiteMesh = seedVertex_index[kmed]
      kmedoidsTime_tex[stepCount][kmedIndexInWhiteMesh]=k
    stepCount += 1

  return kmedoidsTime_tex, vertex_silhouette_width_Time_tex


def writeClusteringResultsTxt(clusterings_dict, k_ordering, txt_filename):
  """
  Write clustering results information in a .txt filename
  For intra subject clustering results, for one subject
  inputs:
          clusterings_dict: result of clusteringResults(X, k_min, k_max, Rclustering_filename) method
          k_ordering: result of clusteringResults(X, k_min, k_max, Rclustering_filename) method
          txt_filename: input txt filename
  """
  file = open(txt_filename, 'w')
  k_nb = len(k_ordering)
  file.write( "---\n" + time.strftime( '%d/%m/%y %H:%M',time.localtime() ) + "\n---\n" )
  file.write("k ordering: ")
  for i in xrange(k_nb):
    k = k_ordering[i]
    if clusterings_dict.has_key(k):
      if i != k_nb-1:
        file.write(str(k) + ", ")
      else:
        file.write(str(k) + "\n")

  for i in xrange(k_nb):
    k = k_ordering[i]
    if clusterings_dict.has_key(k):
      file.write("k = " + str(k) + ": " + str(clusterings_dict[k]['avg.width']) + "\n")
  file.close()


def addOneSubjectIntraClusteringResultsTxt(clusterings_dict, k_ordering, 
                                           txt_filename, gyrus ):
  """
  Write all clustering results information in a .txt filename
  add intra subject clustering results, to txt_filename (sum up for a set of subjects)
  inputs:
          clusterings_dict: result of clusteringResults(X, k_min, k_max, Rclustering_filename) method
          k_ordering: result of clusteringResults(X, k_min, k_max, Rclustering_filename) method
          txt_filename: input txt filename
  """
  precision = pydeci.Decimal("0.001")
  file = open(txt_filename, 'a')
  k_nb = len(k_ordering)
  file.write( "---\n" + time.strftime( '%d/%m/%y %H:%M',time.localtime() ) + "\n---\n" )
  file.write( "Gyrus " + str(gyrus) + "\n" )
  for i in xrange(k_nb):
    k = k_ordering[i]
    if clusterings_dict.has_key(k):
      avg_width = str(clusterings_dict[k]['avg.width'])
      avg_width_good_prec = str(pydeci.Decimal(avg_width).quantize(precision))
      if i != k_nb-1:
        file.write("k=" + str(k) + ": " + avg_width_good_prec + ", ")
      else:
        file.write("k=" + str(k) + ": " + avg_width_good_prec + "\n")
  file.close()


