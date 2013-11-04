
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