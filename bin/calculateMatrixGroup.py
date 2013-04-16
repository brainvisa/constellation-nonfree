#!/usr/bin/env python

import optparse
from soma import aims
import numpy
import sys
import os

import roca.lib.interMeshParcellation.processes_lib as T
import roca.lib.textureTools_stats as TTS

def parseOpts(argv):
  desc = """Calculate Matrix Group: average or concatenate."""
  
  parser = optparse.OptionParser(desc)

  parser.add_option( '-m', '--matrices', dest='individual_matrix', action='append', help='All individual matrices reduced to the correct group.' )
  parser.add_option( '-o', '--output', dest='matrix', help='matrix group' )
  parser.add_option( '-s', '--study', dest= 'study', help='It\'s a String. Choice "Average" or "Concatenate".' )

  return parser, parser.parse_args(argv)

def main():
  parser, (options, args) = parseOpts(sys.argv)

  files = options.individual_matrix
  
  count = 0
  for matrixf in files:
    subject_reducedConnMatrix = T.readMatrixImaAsNumpyArray( str(matrixf), True )
    print 'averaging', matrixf, ', min/max:', numpy.min( subject_reducedConnMatrix ), numpy.max( subject_reducedConnMatrix )
    if options.study == 'Concatenate':
      if count == 0:
        subjectsReducedConnMatrix = subject_reducedConnMatrix
        #subjectsReducedConnMatrix = TTS.writeConnMatrixAsIma(subjectsReducedConnMatrix, options.matrix )
      else:
        #subjectsReducedConnMatrix = T.readMatrixImaAsNumpyArray(options.matrix, True)
        subjectsReducedConnMatrix = numpy.concatenate( ( subjectsReducedConnMatrix, subject_reducedConnMatrix ) )
        #subjectsReducedConnMatrix = TTS.writeConnMatrixAsIma(subjectsReducedConnMatrix, options.matrix )
      count += 1
      TTS.writeConnMatrixAsIma(subjectsReducedConnMatrix, options.matrix )
    if options.study == 'Average':
      if count == 0:
        avgsubjectsReducedConnMatrix = subject_reducedConnMatrix
      else:
        avgsubjectsReducedConnMatrix = avgsubjectsReducedConnMatrix + subject_reducedConnMatrix
      count += 1
  if options.study == 'Average':
    avgsubjectsReducedConnMatrix = (1./count)*avgsubjectsReducedConnMatrix
    avgReducedConnMatrix = TTS.writeConnMatrixAsIma( avgsubjectsReducedConnMatrix, options.matrix )
  print 'matrix : ', options.matrix 
  print 'count : ', count
  print 'sum min/max:', numpy.min( avgsubjectsReducedConnMatrix ), numpy.max( avgsubjectsReducedConnMatrix )
  print 'average min/max:', numpy.min( avgsubjectsReducedConnMatrix ), numpy.max( avgsubjectsReducedConnMatrix )

if __name__ == "__main__" : main()