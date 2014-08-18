#!/usr/bin/env python

# 
import numpy as np
import itertools

#soma
from soma import aims
   
    
def resize_matrix(matrix, rows_length=100.0, cols_length=80.0):
    """Resize a matrix.
    
    The visual effect may be easier.    
    
    Parameters
    ----------
        matrix: (aims.Volume)
        rows_length: fixed value
        cols_length: fixed value
        
    Return
    ------
        matrix: (aims.Volume)
    """
    if matrix.getSizeX() == matrix.getSizeY():
      matrix.header()['voxel_size'] = aims.vector_FLOAT(
          [rows_length / matrix.getSizeX(), rows_length / matrix.getSizeY(), 1,1])
    else:
      matrix.header()['voxel_size'] = aims.vector_FLOAT(
          [rows_length / matrix.getSizeX(), cols_length / matrix.getSizeY(), 1,1])
    
    return matrix

def permutationResampling(feat):
  '''Resampling by permutation of the features'''
  Nsamples = feat.shape[0]
  Ndim = feat.shape[1]  
  perm = np.random.permutation(feat[:, 0].reshape((Nsamples, 1)))
  for f in range(1, Ndim):
    perm = np.hstack((perm, np.random.permutation(feat[:, f].reshape((Nsamples, 1))) ))
  return perm
  
def caseResampling(feat):
  '''
  newFeat: new features matrix
  bootVect: list of indices resampling
  survive: list of items that survived the resampling
  '''
  N = feat.shape[0]
  bootVect = np.random.random_integers( N, size = N )
  bootVect = bootVect - 1
  newFeat = feat[bootVect[0], :]
  for i in range(1, bootVect.size):
    newFeat = np.vstack( ( newFeat, feat[bootVect[i], :] ) )
  survive = np.array([])
  for i in bootVect:
    if ( ( np.where(survive == i)[0] ).size == 0 ):
      survive = np.append(survive, np.array([int(i)]))
  return newFeat, bootVect, survive
  
def GenerateUniform(whiteFeat, Nsample, Ndim):
  '''
  Generate samples from uniform distrib with same shape than original samples.
  '''
  n0 = ( whiteFeat[:, 0].max() - whiteFeat[:, 0].min()) * np.random.random_sample( ( Nsample, 1 ) ) + whiteFeat[:, 0].min()
  n1 = ( whiteFeat[:, 1].max() - whiteFeat[:, 1].min()) * np.random.random_sample( ( Nsample, 1 ) ) + whiteFeat[:, 1].min()
  s = np.hstack( (n0, n1) )
  for i in range(2, Ndim):
    ni = (whiteFeat[:,i].max() - whiteFeat[:, i].min()) * np.random.random_sample( ( Nsample, 1 ) ) + whiteFeat[:, i].min()
    s = np.hstack( (s, ni) )
  print 'Generated uniform matrix of shape: ', s.shape
  return s
  
def orderDataMatrix(mat, labels):
  """
  inputs:
            mat: matrix of observations, shape = (n,p), n observations, p attributes
            labels: numpy array of shape (n) corresponding to the labels of observations
            
  output:   order_mat: ordered matrix, of shape (n,p)
            labels: ordered labels
  """
  (n,p) = mat.shape
  obs_nb = labels.size
  if n != obs_nb:
    raise exceptions.ValueError("matrix dimensions and labels size are not compatible")
  labels_argsort = labels.argsort()
  order_mat = np.zeros((n,p),dtype = np.float32)
  for i in xrange(n):
    order_mat[i,:] = mat[labels_argsort[i],:]
  sortLabels = labels.copy()
  sortLabels.sort()
  return order_mat, sortLabels

def euclidianDistance(v1,v2):
  """
  input:
          v1, v2: two vectors, shape (1,p) (numpy arrays)
  output:
          dist(v1,v2)
  """
  dist = ((v1-v2)**2).sum()
  return np.sqrt(dist)
   
def euclidianDistanceMatrix(matrix):
  (n,p) = matrix.shape
  euclidian_dist_matrix = np.zeros((n,n), dtype = np.float)
  for i in xrange(n):
    v1 = matrix[i]
    dist_value = 0
    euclidian_dist_matrix[i][i] = 0
    for j in xrange(0,i):
      v2 = matrix[j]
      dist_value = euclidianDistance(v1,v2)
      euclidian_dist_matrix[i][j] = dist_value
      euclidian_dist_matrix[j][i] = dist_value
  return euclidian_dist_matrix

def compute_mclusters_by_nbasins_matrix(reducedmatrix, parcels, timestep = 0, mode = "meanOfProfiles"):
  """
  Compute the mean connectivity profile of each parcels and create the associated matrix:
  inputs:
    reducedmatrix: size (patchVertex_nb, targets_nb)
    vertices_patch: size (patchVertex_nb): index of the patch vertex in the parcels
    parcels: labeled aims Time Texture_S16, parcellation of the patch
    timestep: if the input parcels has several time steps, indicates the chosen step. (each step corresponding to a parcellation of the patch, time_step = 0, parcellation into 2 clusters)
    mode: "mean": normalized by the number of parcels vertices
       or "sum": not normalized
  outputs:
    matrix: size (parcels_nb, targets_nb)
  """
  (n,p) = reducedmatrix.shape
  print "reducedmatrix.shape", reducedmatrix.shape
  labels = np.zeros((n), dtype=int)
  
  vertices_patch = []

  for i in xrange( parcels[0].nItem() ):
    if parcels[0][i] != 0:
      vertices_patch.append(i)
      
  vertices_patch = np.array( vertices_patch ) 
  vertices_patch =  vertices_patch.tolist()
  for i in xrange(n):#iteration on patch Vertex
    index_i = np.uint32(vertices_patch[i])
    labels[i]= parcels[timestep][index_i]
  
  labels_unique = np.unique(labels).tolist()
  
  matrix = np.zeros((len(labels_unique),p),dtype=np.float32)
  labelCount = 0
  for i in xrange(len(labels_unique)):
    label = labels_unique[i]
    label_connMatrixToTargets_array = reducedmatrix[np.where(labels==label)]
    if mode == "meanOfProfiles":
      matrix[labelCount]=label_connMatrixToTargets_array.mean(axis = 0)
    elif mode == "sumOfProfiles":
      matrix[labelCount]=label_connMatrixToTargets_array.sum(axis = 0)
    else:
      raise exceptions.ValueError("the mode must be \"mean\" or \"sum\" ")
    labelCount+=1
  
  return matrix
  
def contingency_matrix( labels1, labels2 ):
  classes = list( set( labels1 ) )
  n = len( classes )
  contingency_matrix = np.array( [ zip(labels1, labels2 ).count(x) for x in itertools.product( classes, repeat = 2 ) ] ).reshape( n, n )
  #contingency_matrix = np.bincount( n * ( labels1 - 1 ) + ( labels2 - 1 ), minlength = n * n).reshape( n, n )
  return contingency_matrix
  
def partialWhiten(features):
  '''
  Normalize the data (sigma = 1) per category.
  '''
  dist = features[:, 0:6]
  direc = features[:, 6:12]
  sdist = np.std( dist )
  sidrec = np.std( direc )
  direc = direc / sidrec
  dist = dist / sdist
  white = hstack( ( dist, direc ) )
  return white