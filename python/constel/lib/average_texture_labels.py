from soma import aims
import numpy
import sys

def usage():
  print "Create average texture labels"
  print "usage : python average_texture_labels.py output.gii(.tex) subject1.gii(.tex) ... subjectN.gii(.tex)"

def average_texture_labels( output, inputs ):
  # read textures
  tex = []
  for fname in inputs:
    tex.append( aims.read( fname ) )
  # make a 2D array from a series of textures
  ar = numpy.vstack( [ t[0].arraydata() for t in tex ] )
  # count occurrences
  N = numpy.max(ar)
  def bin_resized( x ):
    y = numpy.bincount(x)
    y.resize( N+1 ) # labels: 1 to 72
    return y
  cnt = numpy.apply_along_axis( bin_resized, 0, ar )
  # get max of occurrences in each vertex
  maj = numpy.argmax( cnt, axis=0 )
  # make an aims texture from result (numpy array)
  otex = aims.TimeTexture( 'S16' )
  otex[0].assign( maj )

  aims.write( otex, output )
  
if __name__ == "__main__":
    if len(sys.argv)<=3:
        usage()
        sys.exit(1)

    li_labels = []
    output = sys.argv[1]
    inputs = sys.argv[2:]
    average_texture_labels(output, inputs)