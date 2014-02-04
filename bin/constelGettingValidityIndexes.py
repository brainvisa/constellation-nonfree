#!/usr/bin/env python

from optparse import OptionParser
import Pycluster as pc
from soma import aims
import numpy as np

import constel.lib.clustervalidity as cv

def parseOpts(argv):
  description = 'Cluster criterion: David-Bouldin Index, Dunn Index and Calinski-Harabasz Index.'

  parser = OptionParser(description)

  parser.add_option( '-m', '--matrix', dest = 'matrix', help = 'Reduced connectivity matrix (Ndim, Nsample)' )
  parser.add_option( '-n', '--nbiter', dest = 'nbiter', help = 'Number of iterations of same kmedoid' )
  parser.add_option( '-k', '--kmax', dest = 'kmax', help = 'K max for the clustering' )
  parser.add_option( '-f', '--file', dest = 'indexFile', help = 'Validity indexes file' )
  
  return parser, parser.parse_args(argv)
  
def main():
  parser, (options, args) = parseOpts(sys.argv)
  
  matrixFile = str(options.matrix)
  matrix = aims.read( options.matrix )
  feat = np.asarray(matrix)[:, :, 0, 0].transpose()
  Nsample = feat.shape[0]
  Ndim = feat.shape[1]
  print 'Number of samples: ', Nsample, 'and Dimension: ', Ndim
  
  nbIter = options.nbiter
  kmax = options.kmax
  
  distance1 = pc.distancematrix(feat, dist='e')
  
  print 'Euclidean distance at power 2'
  distMat1 = [array([])]
  for i in range(1, feat.shape[0]):
    distMat1.append((distance1[i] ** 2).copy())
  
  print 'Converting distance matrix'
  distMat = np.zeros((Nsample, Nsample))
  for i in range(1, Nsample):
    for j in range (0, i):
      distMat[i,j] = distMat1[i][j]
      distMat[j,i] = distMat1[i][j]
  Nsample = distMat.shape[0]  
  print 'There are', Nsample, 'samples'
  
  sDunn = zeros(Kmax + 1)
  sDB = zeros(Kmax + 1)
  sCH = zeros(Kmax + 1)
  
  for k in range(Kmax):
    print 'Runnning K-medoids with K = ', K
    uniclusterid, unierror, uninfound = pc.kmedoids(distMat, k + 1, nbIter)
    print 'getting validity indexes'
    sDB[k + 1], sDunn[k + 1], sCH[k + 1] = cv.computeClusterValidity(distMat, uniclusterid, k + 1)
    print 'sDB = ',sDB[k + 1]
    print 'sDunn = ', sDunn[k + 1]
    print 'sCH = ', sCH[k + 1]
  
  sDB[1] = sDB.min()
  sDB[0] = sDB.min()
  sDunn[1] = sDunn.min()
  sDunn[0] = sDunn.min()
  sCH[1] = sCH.min()
  sCH[0] = sCH.min()
  
  print 'Normalizing scores'
  sDB = (sDB - sDB.min()) / float(sDB.max() - sDB.min())
  sDunn = (sDunn - sDunn.min()) / float(sDunn.max() - sDunn.min())
  sCH = (sCH - sCH.min()) / float(sCH.max() - sCH.min())
  
  #scores = np.vstack((sDB, sCH))
  scores = vstack((sDunn, sDB))
  scores = vstack((scores, sCH))
        
  SFA = mean(scores, axis = 0)
  SFG = st.mstats.gmean(scores, axis = 0)
  SFMed = median(scores, axis = 0)
  
  fileR = open(options.indexFile, 'w')
  fileR.write( 'Reduced connectivity matrix is ' + matrixFile + '\n')
  linesDB = ''
  linesDunn = ''
  linesCH = ''
  lineSFA = ''
  lineSFG = ''
  lineSFMed = ''
  
  linesDB = linesDB + str(sDB) + ' '
  linesDunn = linesDunn + str(sDunn) + ' '
  linesCH = linesCH + str(sCH) + ' '
  lineSFA =lineSFA + str(SFA) + ' '
  lineSFG = lineSFG + str(SFG) + '' 
  lineSFMed = lineSFMed + str(SFMed) + ''
  
  fileR.write( 'Clustering validity index --> sDB: ' + linesDB + '\n')
  fileR.write( 'Clustering validity index --> sDunn: ' + linesDunn + '\n')
  fileR.write( 'Clustering validity index --> sCH: ' + linesCH + '\n')
  fileR.write( 'Fusion-based scores --> SFA: ' + lineSFA + '\n')
  fileR.write( 'Fusion-based scores --> SFG: ' + lineSFG + '\n')
  fileR.write( 'Fusion-based scores --> SMed: ' + lineSFMed + '\n')

  fileR.close()
