#!/usr/bin/env python
from soma import aims
import numpy as np
import Pycluster as pc
from scipy.spatial.distance import pdist, squareform
import plot as p

def silhouette_sample(X, K):
  n = X.shape[0]
  dist = pc.distancematrix(X, dist='e')
  NiterK = 100
  clusterid, error, nfound = pc.kmedoids(dist, K, NiterK)
  labels = 0
  for i in np.unique(clusterid):
    clusterid[clusterid == i] = labels
    labels += 1
  cl = len(np.unique(clusterid))
  distance = squareform(pdist(X))
  kIndices = [np.flatnonzero(clusterid == k) for k in range(cl)]
 
  a = np.zeros(n)
  b = np.zeros(n)
  for i in range(n):
    a[i] = np.mean( [distance[i][ind] for ind in kIndices[clusterid[i]] if ind!=i] )
    b[i] = np.min( [np.mean(distance[i][ind]) for k,ind in enumerate(kIndices) if clusterid[i] != k] )
  s = (b - a)/np.maximum(a, b)
  #p.silhouette_plot(X, s, clusterid)
  return s
  
def silhouette_score(X, K):
  s = silhouette_sample(X, K)
  score = np.mean(s)
  return score