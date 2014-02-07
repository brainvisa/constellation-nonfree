#!/usr/bin/env python
import numpy as np
import colorsys
from soma import aims
import random

def generateIntPairsNames(elements_nb):
  """
  input:
          elements_nb: nb of elements (elements generated are int labels between 0 and elements_nb)
  output:
          pairsNames_list: list of all the combination (of two elements) possible : for example : = [0_0, 0_1, 0_2, 1_1, 1_2, 2_2] if elements_nb == 2
          names_dict: for all the labels between 1 and elements_nb, list of the names they belong to. For example, names_dict[1] = [0_1, 1_1, 1_2]
  """
  pairsNames_list = []
  names_dict = {}
  for i in xrange(elements_nb+1):
    label1 = i
    new_name = str(label1)+"_"+str(label1)
    pairsNames_list.append(new_name)
    names_dict[label1] = []
    names_dict[label1].append(new_name)
    for j in xrange(min(elements_nb, i)):
      label2 = j
      new_name = str(min(label1, label2)) + "_" + str(max(label1, label2))
      pairsNames_list.append(new_name)
      names_dict[label1].append(new_name)
      names_dict[label2].append(new_name)
  #print pairsNames_list
  return pairsNames_list, names_dict
    
def generateIntPairsNamesfromNamesList(names_list):
  pairsNames_list = []
  names_dict = {}
  for label1 in names_list:
    new_name = str(label1)+"_"+str(label1)
    pairsNames_list.append(new_name)
    names_dict[label1] = []
    names_dict[label1].append(new_name)
    for label2 in names_list:
      if names_dict.has_key(label2) == False:
        names_dict[label2] = []
      if label2 < label1:
        new_name = str(min(label1, label2)) + "_" + str(max(label1, label2))
        pairsNames_list.append(new_name)
        names_dict[label1].append(new_name)
        names_dict[label2].append(new_name)
  #print pairsNames_list
  return pairsNames_list, names_dict
      
def addStringToFilename(filename, added_str):
  print "filename:", filename
  split_list=filename.split('.')
  out_filename = filename.replace(split_list[-2], split_list[-2]+added_str)
  return out_filename

def reductionHisto(input_histo, input_lengths, output_lengths):
  """
  inputs:
    input histo: array of size (lengths.size, vertices_nb)
    input_lengths = N.asarray(range(20,201))#absciss of input histogram, size = 181, type: list
  outputs:
    output_lengths = N.asarray(range(2,11))*20# absciss of output histogram = [40, 60, 80, 100, 120, 140, 160, 180, 200], size = 10, type: list
  example:
input_histo = N.ones((10,4))
for i in range(10):
  for j in xrange(4):
    input_histo[i,j]=i

input_lengths=N.asarray(range(1,11))
output_lengths=N.asarray(range(1,5))*2
output_histo = reductionHisto(input_histo, input_lengths, output_lengths)
    """
  inHisto_dims = input_histo.shape
  output_histo = np.zeros((output_lengths.size,inHisto_dims[1]))
  valid = False
  for i in xrange(output_lengths.size):
    length = output_lengths[i]
    if i==0:
      valid = input_lengths<=length
    else:
      validMaxlength = input_lengths<=length
      validMinLength = input_lengths >output_lengths[i-1]
      valid = validMaxlength*validMinLength
    coords_in_input_histo = np.where(valid)
    output_histo[i,:]=input_histo[valid].sum(axis = 0)
  return output_histo

def sameNbElements(listA, listB, NonZeroRealNb=True):
  '''Return lists with the same number of elements'''
  for index in range(len(listA)):
    if listA[index] == 0 and listB[index] != 0:
      listB[index] = 0
    elif listA[index] != 0 and listB[index] == 0:
      listA[index] = 0
  if NonZeroRealNb:
    listA = [x for x in listA if x != 0]
    listB = [x for x in listB if x != 0]
  if len(listA) != len(listB):
    print 'ERROR: The A size is different from the B.'
  return listA, listB


def makeColormapGradients( cmap ):
  # SV variants in HSV space
  cols = [ [ 0.2,0.33 ], [ 0.2,.5], [ 0.2,0.66],
    [ 0.6,0.33], [0.6,0.44], [0.6, 0.55], [ 0.6, .66],
    [ 1., 0.33], [1., 0.5], [1., .66] ]
  outcmap = aims.Volume_RGB( cmap.getSizeX() * len( cols ) )
  index = 0
  for i in xrange( cmap.getSizeX() ):
    color = cmap.value( i )
    hsv = colorsys.rgb_to_hsv( color[0]/255., color[1]/255., color[2]/255. )
    col0 = hsv[0]
    coltable = [ [ col0, s, v ] for s,v in cols ]
    rgb = [ colorsys.hsv_to_rgb( *x ) for x in coltable ]
    for i, x in enumerate( rgb ):
      outcmap.setValue( aims.AimsRGB(
        int(x[0]*256), int(x[1]*256), int(x[2]*256) ), index )
      index += 1
  return outcmap


def makeDoubleRandomHueColormap( n = 36 ):
  hue = list( np.arange( n ) / float(n) )
  random.shuffle( hue )
  hue2 = hue + hue
  rgb = [ colorsys.hsv_to_rgb( x, 0.7, 0.7 ) for x in hue2 ]
  vol = aims.Volume_RGB( n * 2 )
  for i, x in enumerate( rgb ):
    vol.setValue( aims.AimsRGB( int(x[0]*255.9), int(x[1]*255.9),
      int(x[2]*255.9) ), i )
  return vol

def crop(img, x1, x2, y1, y2):
    """Return the cropped image at the x1, x2, y1, y2 coordinates
    """
    if x2 == -1:
        x2 = img.shape[1]-1
    if y2 == -1:
        y2=img.shape[0]-1
    mask = np.zeros(img.shape)
    mask[y1:y2 + 1, x1:x2 + 1] = 1
    m = mask > 0
    return img[m].reshape((y2 + 1 - y1, x2 + 1 - x1))