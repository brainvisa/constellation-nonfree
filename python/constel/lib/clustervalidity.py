#!/usr/bin/env python
from soma import aims
import numpy as np
import Pycluster as pc

def silhouette_sample(X, K):
  n = X.shape[0]
  dist = pc.distancematrix(X, dist='e')
  NiterK = 100
  clusterid, error, nfound = pc.kmedoids(dist, K, NiterK)
  cl = len(np.unique(clusterid))
  
  # compute pairwise distance matrix
  fact = 2
  distance = [np.array([])]
  for i in range(1, X.shape[0]):
    distance.append((dist[i]**fact).copy())
  
  # indices belonging to each cluster
  kIndices = [np.flatnonzero(clusterid == k) for k in range(cl)]
  print kIndices
  
  # compute a,b,s for each instance
  a = np.zeros(n)
  b = np.zeros(n)
  for i in range(n):
    # instances in same cluster other than instance itself
    a[i] = np.mean( [distance[i][ind] for ind in kIndices[clusterid[i]] if ind!=i] )
    b[i] = np.min( [np.mean(distance[i][ind]) for k,ind in enumerate(kIndices) if clusterid[i] != k] )
  s = (b - a)/np.maximum(a, b)
  return s
  
def silhouette_score(X, labels):
  pass