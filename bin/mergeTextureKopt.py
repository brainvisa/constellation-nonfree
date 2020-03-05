#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from optparse import OptionParser
from soma import aims
import numpy as np
import sys
import six

import constel.lib.utils.texturetools as TT

def usage():
  print("Fusion multi texture for all gyri (by subject)")
  print("usage: python mergeTextureKopt.py output.tex file1.tex ... fileN.tex")

def parseOpts(argv):
  description = 'Fusion of a texture in another.'
  parser = OptionParser(description)
  parser.add_option('-i', '--texIn', dest='texIn',
    metavar = 'FILE',action='store', default = None,
    help='input texture name 1')
  parser.add_option('-o', '--texOut', dest='texOut',
    metavar = 'FILE', help='input tex name 2')

  return parser, parser.parse_args(argv)

def main():
  parser, (options, args) = parseOpts(sys.argv)
  l = 0
  j = 0
  disp = 10 #for id_scale, we want to keep the same numbering (gyrus by gryrus)
  for i in options.texIn:
    if l == 0:
      tex_init = aims.read( options.texIn[0] )
    else:
      tex_init = aims.read( tex_final )
    tex_destination = aims.read( options.texIn[j+1] )
    tex_prov = TT.immerseTexInAnother( tex_init, tex_destination, disp ) # textureTools
    tex_final = aims.write( tex_prov, options.texOut )
    l = 1
    j += 1
    disp += 10

if __name__ == "__main__" : main()

#def mergeTextureKopt( output, inputs ):
  #l = 0
  #j = 0
  #disp = 10 #for id_scale, we want to keep the same numbering (gyrus by gryrus)
  #for i in inputs:
    #if l == 0:
      #tex_init = aims.read( inputs[0] )
    #else:
      #tex_init = aims.read( tex_final )
    #tex_destination = aims.read( inputs[j+1] )
    #tex_prov = TT.immerseTexInAnother( tex_init, tex_destination, disp ) # textureTools
    #tex_final = aims.write( tex_prov, output )
    #l = 1
    #j += 1
    #disp += 10

#if __name__ == "__main__":
  #if len(sys.argv)<=2:
        #usage()
        #sys.exit(1)
    #output = sys.argv[1]
    #inputs = sys.argv[2:]
    #mergeTextureKopt( output, inputs )










































  
  #ntex = aims.TimeTexture_S16()
  #t_init = aims.read( inputs[0] )
  #vertex_nb = t_init[0].nItem()
  #ntex.reserve(vertex_nb)
  #dispersion = 0 #for id_scale, we want to keep the same numbering (gyrus by gryrus)
  #count = 0
  
  #for t in inputs:
    #tex = aims.read( t )
    #tex_array = tex[0].arraydata()
    #print(tex_array)
    #for v in six.moves.xrange(vertex_nb):
      #ntex[0].push_back(0)
    #final_valid_vertex_sum = np.ones((0))
    #final_valid_vertex_sum = ntex[0].arraydata()
    #valid_vertex = tex !=0
    #final_valid_vertex_sum[ valid_vertex ] += 1
    #if count > 0:
      
      
      #texture_fusion = tex_array
      #print(texture_fusion)
    #else:
      #print('je rentre dans c=autre')
      #texture_fusion = np.concatenate((texture_fusion, tex_array))
      #print(texture_fusion, len(texture_fusion))
    #count += 1
  #valid_vertex = np.where( tex_array != dispersion )
  #ntex = aims.TimeTexture( 'S16' )
  #ntex.assign(texture_fusion)[0]
  #aims.write( ntex, '/tmp/tex.gii' )
    #vertex_nb = tex[0].nItem()
    #ntex[0].reserve( vertex_nb )
    #for i in six.moves.xrange( vertex_nb ):
      #ntex[0].push_back(0)
    #tex_sum = ntex[0].arraydata()
    #valid_vertex = tex_array != 0
    #ntex.assign( np.where( valid_vertex )[0] )




   #print(np.sum( tex_array ))
   #tex_array_filter = [ x for x in tex_array if x!=0  ] 