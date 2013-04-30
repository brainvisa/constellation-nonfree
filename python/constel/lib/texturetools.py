#!/usr/bin/env python
from soma import aims
import numpy as np


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
  outtex_kept_labels = np.unique(outtex_ar)
  outtex_kept_labels_list = outtex_kept_labels.tolist()
  if outtex_kept_labels_list.count(0) != 0:
    outtex_kept_labels_list.remove(0)
  for i in xrange(len(outtex_kept_labels_list)):
    current_label = outtex_kept_labels_list[i]
    print "current label:", current_label, " new:", str(i+1)
    outtex_ar[tex_ar==current_label] = i+1
  return out_tex


