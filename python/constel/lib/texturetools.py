#!/usr/bin/env python
import soma.aims.texturetools as MTT
from soma import aims
import numpy

def removeLabelsFromTexture( tex, labels_list ):
  """
  inputs:
          tex: labeled texture between 1 and labelsNb, background = 0, aimsTimeTexture_S16
          labels_list: list of labels to remove        
  output:
          out_tex: labeled texture without region of labels_list, and renumbered between 1 and NewLabelsNb (= LabelsNb - len(labels_list))
  """
  out_tex = aims.TimeTexture_S16()
  vertex_nb = int( tex.nItem() )
  tex_ar = tex[0].arraydata()
  out_tex[0].assign( tex_ar )
  labels_nb = tex_ar.max()
  removedLabels_nb = len( labels_list )
  outtex_ar = out_tex[0].arraydata()
  for label in labels_list:
    outtex_ar[outtex_ar==label]=0
  outtex_kept_labels = numpy.unique( outtex_ar )
  outtex_kept_labels_list = outtex_kept_labels.tolist()
  if outtex_kept_labels_list.count(0) != 0:
    outtex_kept_labels_list.remove(0)
  for i in xrange( len(outtex_kept_labels_list) ):
    current_label = outtex_kept_labels_list[i]
    print "current label:", current_label, " new:", str(i+1)
    outtex_ar[ tex_ar==current_label ] = i+1
  return out_tex


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


def findWrongLabels( mesh, gyriTex ):
  """
  inputs:
         mesh
         gyriTex: gyri texture
  outputs:
         wrong_labels: list of wrong labels
         
         [ cctex: connectedComponentTex: aimsTimeTexture_S16, time step = LabelsNb, for each time step (label in the tex), texture of the connected components corresponding to this label (background = -1, and connected components = values between 1 and ccNb)
         areas_measures = python dictionary, areas_measures[label] = [16.5, 6.0] (numpy array) if label (in tex) has two connected Components 1 and 2 with area = 16.5 and 6.0 respectively, areas are in square mm ]
  """
  meshNeighborsVector = aims.SurfaceManip.surfaceNeighbours( mesh )
  cctex, areas_measures = MTT.connectedComponents( mesh, gyriTex, areas_mode = 1 )
  wrong_labels = []
  for label in areas_measures.keys():
    if areas_measures[ label ].size != 1:
      wrong_labels.append( label )
  return wrong_labels


def changeLabelOfWrongLabelVertices( cclabel, label, gyriTex, meshNeighborsVector, cctexLabel_ar ):
  """
  inputs:
          cclabel: label of connected component in cctexLabel_ar
          label: label of associated vertices in gyri texture
          gyriTex: gyri texture
          meshNeighborsVector : aims.SurfaceManip.surfaceNeighbours(mesh) (mesh = white mesh associated to gyriTex)
          cctexLabel_ar : cctex[label-1].arraydata()
  outputs:
          gyriTex: new gyriTex texture
          winner_label
  """
  indexes=numpy.where( cctexLabel_ar == cclabel )[0]
  neighbourLabels = []
  print 'Nb of wrong indexes: ', indexes.size
  for i in indexes:
    for n in meshNeighborsVector[i]:
      n_label = gyriTex[0][n]
      if n_label != label:
        neighbourLabels.append( n_label )
  neighLabels_list = numpy.unique( neighbourLabels )
  maxCount = 0
  winnner_label = -1
  for neighLabel in neighLabels_list:
    neighLabelCount = neighbourLabels.count( neighLabel )
    if neighLabelCount > maxCount:
      print 'neighbourLabels.count( neighLabel ): ', neighbourLabels.count( neighLabel )
      winnner_label = neighLabel
      maxCount = neighLabelCount
  for i in indexes:
    gyriTex[0][i] = winnner_label
  return gyriTex, winnner_label


def cleanGyriTexture( mesh, gyriTex ):
  """
  inputs:
         mesh
         gyriTex: gyri texture
  outputs:
         gyriTex: new gyri texture

         [ cctex: connectedComponentTex: aimsTimeTexture_S16, time step = LabelsNb, for each time step (label in the tex), texture of the connected components corresponding to this label (background = -1, and connected components = values between 1 and ccNb)
         areas_measures = python dictionary, areas_measures[label] = [16.5, 6.0] (numpy array) if label (in tex) has two connected Components 1 and 2 with area = 16.5 and 6.0 respectively, areas are in square mm ]
  """
  meshNeighborsVector = aims.SurfaceManip.surfaceNeighbours( mesh )
  print 'Mesh neighbors vector: ', meshNeighborsVector
  cctex, areas_measures = MTT.connectedComponents( mesh, gyriTex, areas_mode = 1 )
  wrong_labels = []
  for label in areas_measures.keys():
    if areas_measures[ label ].size != 1:
      wrong_labels.append( label )
  for label in wrong_labels:
    cctexLabel_ar = cctex[ label-1 ].arraydata()
    areas_measures_cc = areas_measures[ label ]
    cc_nb = areas_measures_cc.size
    for l in range( 1, cc_nb ):
      gyriTex, win = changeLabelOfWrongLabelVertices( l+1, label, gyriTex, meshNeighborsVector, cctexLabel_ar )
      print "Changed in ", win
  return gyriTex

def textureTime(n_time, list_number, nb_vertices, vertices_patch, mode, minid=None, maxid=None):
  """
  inputs:
        n_time: 
        list_number: 
        nb_vertices: vertex number for output textures
        vertices_patch: vertex index of the seed region
        mode: 1 = average, 2 = concatenate
  output:
        tex: for each time step, results of the k-th clustering
  """
  count = 0
  tex = aims.TimeTexture_S16()
  for k in range(n_time):
    tex[count].resize(nb_vertices, 0)
    if mode == 1:
      tex[count].arraydata()[vertices_patch] = list_number[k].astype(numpy.int16)
    if mode == 2:
      tex[count].arraydata()[vertices_patch] = list_number[k][minid:maxid].astype(numpy.int16)
    count += 1
  
  return tex