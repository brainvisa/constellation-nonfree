#!/usr/bin/env python
from soma import aims
import numpy as np

def writeConnMatrixAsIma(mat, filename, lines_length = 100.0, cols_length = 80.0):
  print "Writing Matrix As .ima"
  (n,p)=mat.shape
  mat_ima = aims.Volume_FLOAT(n,p,1,1)
  mat_ima_array = mat_ima.arraydata()
  mat_ima_array[0,0,:,:] = mat.transpose()
  if n==p:
    mat_ima.header()['voxel_size'] = aims.vector_FLOAT([lines_length/n, lines_length/p, 1,1])
  else:
    mat_ima.header()['voxel_size'] = aims.vector_FLOAT([lines_length/n, cols_length/p, 1,1])
    
  aims.write(mat_ima, filename)
  print "done."

def permutationResampling(feat):
  '''Resampling by permutation of the features'''
  Nsamples = feat.shape[0]
  Ndim = feat.shape[1]
    
  perm = np.random.permutation(feat[:, 0].reshape((Nsamples, 1)))
  for f in range(1, Ndim):
    perm = np.hstack((perm, np.random.permutation(feat[:, f].reshape((Nsamples, 1))) ))
  return perm
  
def orderDataMatrix(mat, labels):
  """
  inputs:
            mat: matrix of observations, shape = (n,p), n observations, p attributes
            labels: numpy array of shape (n) corresponding to the labels of observations
            
  output:   order_mat: ordered matrix, of shape (n,p)
            labels: ordered labels
  """
  (n,p)=mat.shape
  obs_nb = labels.size
  if n!=obs_nb:
    raise exceptions.ValueError("matrix dimensions and labels size are not compatible")
  labels_argsort = labels.argsort()
  order_mat = np.zeros((n,p),dtype = np.float32)
  for i in xrange(n):
    order_mat[i,:]=mat[labels_argsort[i],:]
  sortLabels= labels.copy()
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
  (n,p)=matrix.shape
  euclidian_dist_matrix = np.zeros((n,n), dtype = np.float)
  for i in xrange(n):
    v1 = matrix[i]
    dist_value = 0
    euclidian_dist_matrix[i][i]=0
    for j in xrange(0,i):
      v2 = matrix[j]
      dist_value = euclidianDistance(v1,v2)
      euclidian_dist_matrix[i][j]= dist_value
      euclidian_dist_matrix[j][i]= dist_value
  return euclidian_dist_matrix
  
def connMatrixParcelsToTargets(reducedmatrix, parcels, timestep = 0, mode = "meanOfProfiles"):
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
    if mode =="meanOfProfiles":
      matrix[labelCount]=label_connMatrixToTargets_array.mean(axis = 0)
    elif mode == "sumOfProfiles":
      matrix[labelCount]=label_connMatrixToTargets_array.sum(axis = 0)
    else:
      raise exceptions.ValueError("the mode must be \"mean\" or \"sum\" ")
    labelCount+=1
  
  return matrix