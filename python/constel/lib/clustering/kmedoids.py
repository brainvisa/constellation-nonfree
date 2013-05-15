#!/usr/bin/env python
import numpy as np
#import rpy
#try:
  #import rpy2.rpy_classic as rpy
#except:
  #import rpy2.rpy_classic as rpy

from rpy2 import rinterface
import rpy2.robjects as robjects
#from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
if hasattr( rpy2.robjects.numpy2ri, 'activate' ):
  # new versions of rpy2 need this.
  rpy2.robjects.numpy2ri.activate()

class KMEDOIDS_R(object):
  def __init__(self, clustering_file_name):
    #s = 'source(' + clustering_file_name + ')'
    #rpy2.robject.r( s )
    robjects.r.source(clustering_file_name) #to source rpy.r.clustering function
  #def fit(self,X,Y,k):
    #X = np.asarray(X)
    #Y = np.asarray(Y)
    #(n,p)=X.shape
    #centroids = numpy.zeros((k,p),float)
    ##Calcul of centroids:
    #for i in xrange(k):
      #centroids[i,:]=

  def predict(self,X,k,display=0, diss=False, init_kmedoids=None):
    """
    X: data matrix (size:(n,p)) or dissimilarity matrix (size (n,n))
    diss: dissimilarity matrixmode: if true : X is the dissimilarity matrix otherwise, X is a data matrix: each row corresponds to an observation, and each column corresponds to a variable.
    init_medoids: numpy array of length k of integer indices (in 1:n) specifying initial medoids instead of using the build algorithm
    """
    X = np.asarray(X)
    print "X shape:", X.shape
    print "k in R:", k
    print "start binding avec R : Computing kmedoids algo"
    
    if diss:
      if init_kmedoids:
        cl = robjects.r.kmedoidsondist(X,k,init_kmedoids )[0]
      else:
        cl = robjects.r.kmedoidsondist2(X,k)
      #cl = robjects.r.pam(x, k, diss=TRUE, medoids=init_medoids)
    else:
      #cl = robjects.r.coucoukmedoids(X,k,display)[0] #type(parcellData) = dict
      cl = robjects.r.coucoukmedoids(X,k,display)
    print type(cl)
    #print cl
    #print avg_width
    dict = {}
    dict['silinfo']={}
    dict["id.med"]=np.asarray(np.array(tuple(cl)[1]))
     #- 1
    #print tuple(tuple(cl)[6])
    #print tuple(tuple(cl)[6])[0]
    dict["clustering"]=np.array(tuple(cl)[2])
    dict['silinfo']['widths']=np.array(tuple(tuple(cl)[6])[0])
    dict['silinfo']['clus.avg.widths']=np.array(tuple(tuple(cl)[6])[1])
    dict['silinfo']['avg.width']=np.array(tuple(tuple(cl)[6])[2])[0]
    #cl_medoids_id = np.asarray(cl["id.med"]) - 1 #index array in R are between 1 and array.size
    #cl_silinfo = cl['silinfo']['avg.width']
    #return cl["medoids"],cl_medoids_id, np.asarray(cl["clustering"]),cl_silinfo['avg.width']
    return dict
    
  def predictOnlyMedoids(self,X,k):
    X = np.asarray(X)
    cl = robjects.r.kmedoidsWithDistMatrix(X,k)
    return np.asarray(np.array(tuple(cl)[1]))
  
    

    
    