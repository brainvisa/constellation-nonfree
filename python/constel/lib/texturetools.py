#!/usr/bin/env python
from soma import aims
import numpy


def removeLabelsFromTexture(tex, labels_list):
  """
  inputs:
          tex: labeled texture between 1 and labelsNb, background = 0, aimsTimeTexture_S16
          labels_list: list of labels to remove
          
  output:
          out_tex:
                labeled texture without region of labels_list, and renumbered between 1 and NewLabelsNb (= LabelsNb - len(labels_list))
  """
  out_tex = aims.TimeTexture_S16()
  vertex_nb = int(tex.nItem())
  initTex(out_tex, vertex_nb, 0, 0)
  tex_ar = tex[0].arraydata()
  fillTex(out_tex, tex_ar, 0)
  labels_nb = tex_ar.max()
  removedLabels_nb = len(labels_list)
  outtex_ar = out_tex[0].arraydata()
  for label in labels_list:
    outtex_ar[outtex_ar==label]=0
  outtex_kept_labels = numpy.unique(outtex_ar)
  outtex_kept_labels_list = outtex_kept_labels.tolist()
  if outtex_kept_labels_list.count(0) != 0:
    outtex_kept_labels_list.remove(0)
  for i in xrange(len(outtex_kept_labels_list)):
    current_label = outtex_kept_labels_list[i]
    print "current label:", current_label, " new:", str(i+1)
    outtex_ar[tex_ar==current_label] = i+1
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

