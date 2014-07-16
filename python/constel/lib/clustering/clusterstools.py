#!/usr/bin/env python

# python system modules
from tempfile import mkdtemp
import os.path as path
import numpy as np
import math
import os

# scipy
from scipy.cluster.hierarchy import fcluster

# Pycluster
import Pycluster as pc

# soma
from soma import aims

#fastcluster
import fastcluster


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
        print "labels -->", kml, "centers -->", kmc
        disp = sum([dst(data[m,:],kmc[kml[m],:]) for m in range(shape[0])])
    
      refdisps = scipy.zeros((rands.shape[2],))
      for j in range(rands.shape[2]):
          (kmc,kml) = scipy.cluster.vq.kmeans2(rands[:,:,j], k)
          refdisps[j] = sum([dst(rands[m,:,j],kmc[kml[m],:]) for m in range(shape[0])])
      gaps[i] = scipy.log(scipy.mean(refdisps))-scipy.log(disp)
    
    return gaps
  
def ward_method(dist_mat_file, n, output_dir, n_clusters):
    '''Ward's method is a criterion applied in hierarchical cluster analysis. 
    Ward's minimum variance criterion minimizes the total within-cluster 
    variance.
   
    Parameters:
        dist_mat_file: binary file of the distance matrix
        n: shape[0] of the reduced matrix
        output_dir_Z: directory to save the dendrogram
        n_clusters: k max for the clustering
        
    Returns:
        clusterid: for each vertex, a label is attributed
                   between 1 to k max
   
    '''
   
    # distance matrix file loading
    f = open(dist_mat_file, "r")
    
    # initialization distance matrix as vector
    filename = path.join(mkdtemp(), 'newfile.dat')
    dist_mat = np.memmap(filename, dtype='float32', mode='w+', shape=(n*(n-1)/2,))
    #dist_mat = np.zeros((n*(n-1)/2,))
    
    # number of iteration to generate the distance matrix
    n_iter = 100000 # not perfect
    nb_iteration = dist_mat.shape[0] / n_iter
    nb_add = dist_mat.shape[0] % n_iter
    list_iteration = [nb_iteration] * n_iter
    idx = 0
    for j in xrange(n_iter):
        if j < nb_add:
            list_iteration[j] += 1
        x = f.read(list_iteration[j] * 8)
        dist_mat.ravel()[idx:idx + list_iteration[j]] = struct.unpack('d' * list_iteration[j], x)
        idx += list_iteration[j]
    f.close()

    # compute linkage ward  
    Z = fastcluster.linkage(dist_mat, method='ward', preserve_input=False)
    
    # save the dendrogram
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_dir_Z = output_dir + '/dendrogram.npy'
    np.save(output_dir_Z, Z)
    
    # labelization between 1 to nb
    clusterid = []
    for nb in range(1, n_clusters + 1):
        print "Trying {nb} cluster(s)".format(nb=nb)
        clusters = fcluster(Z, criterion='maxclust', t=nb)
        clusterid.append(clusters)
    return clusterid