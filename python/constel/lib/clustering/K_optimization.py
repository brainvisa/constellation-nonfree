
import constel.lib.clustering.kmedoids as rkmedoids
import numpy as N
import numpy.random as NR
from soma import aims


def clusteringResults(X, k_min, k_max, Rclustering_filename, diss=False,
    init_kmedoids=None):
  """
  inputs:
            X: data = n observations to cluster (shape = (n,p) )
            diss: dissimilarity matrix: if true : X is the dissimilarity matrix otherwise, X is a data matrix: each row corresponds to an observation, and each column corresponds to a variable.
            init_medoids: numpy array of length k of integer indices (in 1:n) specifying initial medoids instead of using the build algorithm
  """
  #k_min = 2
  #Rclustering_filename = "/home/pr216633/p4/roca-stable/R/roca/clustering/kmedoids.R"
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
    print "'avg.width':" , cl['silinfo']['avg.width']
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


