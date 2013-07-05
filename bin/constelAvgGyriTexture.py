#!/usr/bin/env python
from soma import aims
import optparse
import sys

def validation():
  try:
    import constel
  except:
    raise ValidationError( 'constellation module is not here.' )

import constel.lib.texturetools as TT

def parseOpts( argv ):
  description = 'Create average texture labels. usage : python average_texture_labels.py output.gii(.tex) subject1.gii(.tex) ... subjectN.gii(.tex)'
  parser = optparse.OptionParser( description )
  parser.add_option('-i', '--itex', dest = 'itex',
    metavar = 'FILE', action='append',
    help = 'inputs texture (list of)' )
  parser.add_option('-o', '--otex', dest = 'otex',
    metavar = 'FILE',
    help = 'output texture' )
  return parser, parser.parse_args( argv )

def main():
  parser, ( options, args ) = parseOpts( sys.argv )
  TT.average_texture_labels( options.otex, options.itex )

if __name__ == "__main__" : main()