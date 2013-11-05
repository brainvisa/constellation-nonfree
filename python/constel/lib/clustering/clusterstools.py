#!/usr/bin/env python
import numpy as np

def getCenters(clust, K):
  i = 0
  centers = np.array([])
  while (centers.size < K):
    c = clust[i]
    if ((np.where(centers == c)[0]).size == 0):
      centers = np.append(centers, array([int(c)]))
    i += 1
  return centers

def wcDist(distance, clusterid, K, fact):
  centers = getCenters(clusterid, K)
  W = zeros(K)
  addW = 0.0
  for l in range(K):
    sum = 0.0
    centerid = np.where(clusterid == centers[l])[0]
    for i in centerid:
      for j in centerid:
        if (i < j):
          sum += (distance[j][i] ** fact)
        elif (j < i):
          sum += (distance[i][j] ** fact)
    W[l] = sum / float(2.0 * centerid.size)
    addW += W[l]
  return addW