#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def silhouette_sample_plot(X, s, clusterid):
  order = np.lexsort((-s, clusterid))
  K = len(np.unique(clusterid))
  indices = [np.flatnonzero(clusterid[order] == k) for k in range(K)]
  ytick = [(max(ind) + min(ind)) / 2 for ind in indices]
  ytickLabels = ["%d" % x for x in range(K)]
  cmap = cm.jet(np.linspace(0, 1, K) ).tolist()
  clr = [cmap[i] for i in clusterid[order]]
  
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.barh(range(X.shape[0]), s[order], height = 1.0,edgecolor = 'none', color = clr)
  ax.set_ylim(ax.get_ylim()[::-1])
  plt.yticks(ytick, ytickLabels)
  plt.xlabel('Silhouette Value')
  plt.ylabel('Cluster')

  return plt.show()

def silhouette_score_plot(s, k):
  fig = plt.figure()
  plt.plot(k, s, 'k')
  #plt.axis([0, 12, 0, 1.5])
  plt.ylabel('Silhouette Value')
  plt.xlabel('Cluster')
  plt.show()
  return fig