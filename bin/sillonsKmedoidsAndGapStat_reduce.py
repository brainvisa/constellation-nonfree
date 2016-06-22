#!/usr/bin/env python

import constel.lib.utils.matrixtools as clcm
import constel.lib.clustering.clusterstools as clcc
import scipy.spatial.distance as ssd
import scipy.cluster.vq as scv
import Pycluster as pc
from soma import aims
import numpy as np
import optparse
import glob
import sys

def parseOpts(argv):
    desc = """Step 1: cluster the observed data, varying the total number of clusters from k = 1:K, giving within-dispersion measures Wk.
              Step 2: Compute Gap(k)"""
    parser = optparse.OptionParser(desc)
    parser.add_option('-i', '--input', dest='input', metavar='FILE', action='append',
                      help='input permutation files (list)')
    parser.add_option('-m', '--matrix', dest='matrix', metavar='FILE',
                      help='Input connectivity matrix')
    parser.add_option('-k', '--kmax', dest='kmax', type='int',
                      help='Max number of clusters K')
    parser.add_option('-n', '--niter', dest='niter', type='int',
                      help='Number of iterations of same kmedoid algorithm')
    parser.add_option('-o', '--output', dest='output', metavar='FILE',
                      help='Output gap file')
    parser.add_option('-p', '--permutations', dest='permut', type='int',
                      help='number of permutations')                   
    parser.add_option('-e', '--distance', dest='dist', type='str',
                      help='Distance : sqeuclidean, euclidean, cityblock')
    parser.add_option('-f', '--wflag', dest='wflag', type='int',
                      help='Whitening of features or not ? 2 per feature, 1 per category, 0 no')                  
    return parser, parser.parse_args(argv)

def main():
    parser, (options, args) = parseOpts(sys.argv)

    print 'Opening and reading: ', options.matrix
    matrix = aims.read(options.matrix)
    mat = np.asarray(matrix)[:, :, 0, 0].transpose()
    n_sample, n_dim = mat.shape
    print 'Number of samples: ', n_sample, 'Dimension: ', n_dim
    
    if (options.wflag == 2):
        print 'Per feature whitening'
        features = scv.whiten(mat)
    elif (options.wflag == 1):
        print 'Per category (dist and dir) whitening'
        features = clcm.partialWhiten(mat)
    else:
        print 'No whitening'
        features = mat
        
    if (options.dist == 'sqeuclidean'):
        fact = 2
    else:
        fact = 1
    
    distmat = ssd.squareform(ssd.pdist(features, options.dist))
    W = np.zeros(options.kmax + 1)
    print 'Computing Wk for data'
    for K in range(1, options.kmax + 1):
        clusterid, error, nfound = pc.kmedoids(distmat, K, options.niter)
        wc = clcc.wcDist(distmat, clusterid, K, fact)
        print 'K = ', K, ', error = ', error, ', wcDist = ', wc
        #W[K] = np.log(wc)
        W[K] = wc
    print 'Wk = ', W
    
    inputs = glob.glob('*.npy')
    Wkbparts = []
    for infile in inputs:
      Wkbparts.append(np.load(infile))
    Wkb = np.vstack(Wkbparts)
    del Wkbparts
              
    print 'Computing Gap'
    meanWkb = np.mean(Wkb, axis = 0)
    sdWkb = np.std(Wkb, axis = 0)
    gap = np.zeros(options.kmax + 1)
    sW = np.zeros(options.kmax + 1)
    for K in range(1, options.kmax + 1):
        gap[K] = (meanWkb[K] - W[K])
        sW[K] = sdWkb[K] * np.sqrt(1 + 1.0 / float(options.permut))
    
    print 'MeanWkb =', meanWkb
    print 'sdWkb =', sdWkb 
    print 'Gap =', gap
    print 'STD =', sW  
    print 'Writing results to ', options.output
    
    fileR = open(options.output, 'w')
    fileR.write('--> Reading:' + options.matrix + '\n')
    lineW = ''
    lineMW = ''
    lineG = ''
    lineS = ''
    for K in range(1, options.kmax + 1):
        lineW = lineW + str(W[K]) + ' '
        lineMW = lineMW + str(meanWkb[K]) + ' '
        lineG = lineG + str(gap[K]) + ' '
        lineS = lineS + str(sW[K]) + ' '
    fileR.write('--> MeanWkb: ' + lineW + '\n')
    fileR.write('--> sdWkb: ' + lineMW + '\n')
    fileR.write('--> Gap: ' + lineG + '\n')
    fileR.write('--> STD: ' + lineS + '\n')

    fileR.close()

if __name__ == "__main__":
     main()
