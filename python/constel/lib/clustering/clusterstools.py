#!/usr/bin/env python

import numpy as np
import scipy.spatial.distance as ssd
import Pycluster as pc
import math

def getCenters(clust, K):
  i = 0
  centers = np.array([])
  while (centers.size < K):
    c = clust[i]
    if ((np.where(centers == c)[0]).size == 0):
      centers = np.append(centers, np.array([int(c)]))
    i += 1
  return centers

def wcDist(distance, clusterid, K, fact):
  centers = getCenters(clusterid, K)
  W = np.zeros(K)
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

def centroid(dist):
     s = dist.sum(axis = 0)
     c = where(s == s.min())[0][0]
     return c
     
def distToCenters(distance, centers, clusterid, K):
     i = 0
     distCent = zeros(K)
     for i in range(K):
          c = centers[i]
          group = where(clusterid == c)[0]
          dgroup = distance[c, group]
          distCent[i] = dgroup.mean()
     return distCent

def sse(tab):
  """Sum square error function"""
  v = 0.5*np.sum(np.square(tab))
  return v

def entropy( labels ):
  '''Computes entropy of label distribution'''
  n_labels = len( labels )  
  if n_labels <= 1:
    return 0    
  counts = np.bincount( labels )
  counts = np.array( [ x for x in counts if x != 0 ] ).astype(float)
  counts_sum = np.sum( counts )
  probs = counts / n_labels
  n_classes = np.count_nonzero( probs )  
  if n_classes <= 1:
    return 0  
  entropy = -np.sum( probs * ( np.log( counts ) - np.log( counts_sum ) ) )
  return entropy

import scipy
import scipy.cluster.vq
import scipy.spatial.distance
dst = scipy.spatial.distance.euclidean

def gap(data, refs=None, nrefs=20, ks=range(1,11)):
  """
  Compute the Gap statistic for an nxm dataset in data.
  
  Either give a precomputed set of reference distributions in refs as an (n,m,k) scipy array,
  or state the number k of reference distributions in nrefs for automatic generation with a
  uniformed distribution within the bounding box of data.
  
  Give the list of k-values for which you want to compute the statistic in ks.
  """
  shape = data.shape
  if refs==None:
    tops = data.max(axis=0)
    bots = data.min(axis=0)
    dists = scipy.matrix(scipy.diag(tops-bots))
    
    rands = scipy.random.random_sample(size=(shape[0],shape[1],nrefs))
    for i in range(nrefs):
      rands[:,:,i] = rands[:,:,i]*dists+bots
  else:
    rands = refs
  
  gaps = scipy.zeros((len(ks),))
  for (i,k) in enumerate(ks):
    print "ks", ks, "k", k, "i", i
    #clusterid, error, nfound = pc.kmedoids(data, k, 100)
    #labels, error, nfound = pc.kcluster(clusterid, k)
    #kml = labels
    #centers = getCenters(clusterid, k)
    #kmc = centers
    (kmc,kml) = scipy.cluster.vq.kmeans2(data, k)
    print "labels -->", kmc, "centers -->", kml
    disp = sum([dst(data[m,:],kmc[kml[m],:]) for m in range(shape[0])])
    
    refdisps = scipy.zeros((rands.shape[2],))
    for j in range(rands.shape[2]):
      (kmc,kml) = scipy.cluster.vq.kmeans2(rands[:,:,j], k)
      refdisps[j] = sum([dst(rands[m,:,j],kmc[kml[m],:]) for m in range(shape[0])])
    gaps[i] = scipy.log(scipy.mean(refdisps))-scipy.log(disp)
    
  return gaps

def calculate_dispersion(data, clusters, compactness_dispersion):
  """
  W in the Tibshirani paper
  """
  W = np.zeros( len(clusters) + 1)
  Nsamples = data.shape[0]
  distmat = ssd.squareform( ssd.pdist( data, 'euclidean' ) )
  for k in xrange( 1, len( clusters ) ):
    if ( k == 1 ):
      cluster_mean = np.mean( data )
      if compactness_dispersion:
        D = sum( sum( dissmat ) ) / ( 2 * Nsamples )
      else:
        D = sum( sum( data - np.tile( cluster_mean, ( Nsamples, 1) )**2 ) ) / (2 * Nsamples )
      W[k] = math.log( D )
    else:
      clusterid, error, nfound = pc.kmedoids(distmat, k, 100)
      clusterID = []
      labels = 1
      for i in np.unique(clusterid):
        clusterid[clusterid == i] = labels
        labels += 1
      clusterID.append(clusterid)
      clusterid = clusterID
      for j in xrange( k ):
        if compactness_dispersion:
          D[j] = sum( sum( distmat[ np.where( clusterid[j] == k ) ] ) ) / (2 * len( np.where( clusterid == k ) ) )
        else:
          print "j", j
          cluster_mean = np.mean( data[ np.where( clusterid[j] == k ) ] )
          NsamplesC, varDim = ( data[ np.where( clusterid[j] == k ) ] ).shape
          d = sum( ( data[ np.where( clusterid[j] == k) ] - np.tile( cluster_mean, ( NsamplesC, 1 ) ) ) ** 2 )
          print "d", d
          D = sum( d ) / ( 2.0 * NsamplesC )
      W[k] = math.log( sum(D) )
  return W