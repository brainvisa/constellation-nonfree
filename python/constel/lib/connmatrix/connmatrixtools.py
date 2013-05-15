
from soma import aims


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

