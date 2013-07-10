#!/usr/bin/env python
from soma import aims
import optparse
import numpy
import sys

def validation():
  try:
    import constel
  except:
    raise ValidationError( 'constellation module is not here.' )

import constel.lib.texturetools as TT

def parseOpts( argv ):
  description = 'Threshold an aimsTimeTexture (with one time step only ! Warning !)'
  parser = optparse.OptionParser( description )
  parser.add_option('-i', '--itex', dest = 'tex',
    metavar = 'FILE',
    help = 'input gyri texture' )
  parser.add_option('-m', '--mesh', dest = 'mesh',
    metavar = 'FILE',
    help = 'input mesh' )
  parser.add_option('-o', '--otex', dest = 'otex',
    metavar = 'FILE', action = 'store', default = None,
    help = 'output gyri texture with one connected component per gyrus (default : don\'t store resulting texture)' )
  return parser, parser.parse_args( argv )

def main():
  parser, ( options, args ) = parseOpts( sys.argv )
  tex = aims.read( options.tex )
  mesh = aims.read( options.mesh )
  cont = True
  count = 0
  while cont and count < 10:
    print "clean up number ", str( count+1 ), " :"
    tex = TT.cleanGyriTexture( mesh, tex )
    wrong_labels = TT.findWrongLabels( mesh, tex )
    cont = False
    count += 1
    if len( wrong_labels ) > 0:
      cont = True
  aims.write( tex, options.otex )

if __name__ == "__main__" : main()