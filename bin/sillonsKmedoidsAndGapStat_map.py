#!/usr/bin/env python

import constel.lib.connmatrix.connmatrixtools as clcm
import constel.lib.clustering.clusterstools as clcc
import scipy.spatial.distance as ssd
import scipy.cluster.vq as scv
import Pycluster as pc
from soma import aims
import numpy as np
import optparse
import math
import sys

def parseOpts(argv):
    desc = """blabla."""
    parser = optparse.OptionParser(desc)
    parser.add_option('-m', '--matrix', dest='matrix', metavar='FILE',
                       help='input connectivity matrix')
    parser.add_option('-k', '--kmax', dest='kmax', type='int',
                       help='K max')
    parser.add_option('-n', '--niter', dest='niter', type='int',
                       help='number of iterations of same kmedoid')
    parser.add_option('-p', '--permutations', dest='permut', type='int',
                       help='number of permutations')
    parser.add_option('-o', '--output', dest='output', metavar='FILE',
                       help='output permutations file')
    parser.add_option('-r', '--re', dest='typer', type='str',
                      help='Type of resampling: b=bootstrap p=permutations m=montecarlo ')                 
    parser.add_option('-e', '--distance', dest='dist', type='str',
                      help='Distance : sqeuclidean, euclidean, cityblock')
    parser.add_option('-f', '--wflag', dest='wflag', type='int',
                      help='Whitening of features or not ? 2 per feature, 1 per category, 0 no')                    
    return parser, parser.parse_args(argv)
     
def main():
    parser, (options, args) = parseOpts(sys.argv)
   
    print 'Opening and reading: ', options.matrix
    matrix = aims.read(options.matrix)
    matrixar = np.asarray(matrix)[:, :, 0, 0].transpose()
    n_sample, n_dim = matrixar.shape
    print 'Number of samples: ', n_sample, 'Dimension: ', n_dim

    if (options.wflag == 2):
        print 'Per feature whitening'
        features = scv.whiten(matrixar)
    elif (options.wflag == 1):
        print 'Per category (dist and dir) whitening'
        features = clcm.partialWhiten(matrixar)
    else:
        print 'No whitening'
        features = matrixar
    
    if (options.dist == 'sqeuclidean' ):
        fact = 2
    else:
        fact = 1
    
    if (options.typer == 'm'):
        print 'Generating uniform distribution with Monte-carlo'
        uniform = clcm.GenerateUniform(features, n_sample, n_dim)
    elif (options.typer == 'b'):
        print 'Generating uniform distribution that will be bootstrapped'
        uniform = clcm.GenerateUniform(features, n_sample, n_dim)
    elif (options.typer == 'p'):
        print 'Reference distribution will be sampled by permutations of the original one'
        uniform = features.copy()
    
    WB = np.zeros((options.permut, options.kmax + 1))
    print 'Re-sampling ', options.permut, ' times to estimate uniform Wk'
    for iteration in range(options.permut):
        if ((iteration%1) == 0):
            print '  --> Iteration ', iteration
        if (options.typer == 'm'):
            featuni = clcm.GenerateUniform(features, n_sample, n_dim)
        elif (options.typer == 'b'):
            featuni, sampuni, suruni = clcm.caseResampling(uniform)
        elif (options.typer == 'p'):
            featuni = clcm.permutationResampling(uniform)
            print '      Permutation done'
        
        rdist = ssd.squareform(ssd.pdist(featuni, options.dist))
        print '      Distance matrix done'
        for K in range(1, options.kmax + 1):
            uniclusterid, unierror, uninfound = pc.kmedoids(rdist, K, options.niter)
            uniwc = clcc.wcDist(rdist, uniclusterid, K, fact)
            WB[iteration, K] = math.log(uniwc)
        print '      Clustering assessed'
              
    print 'Resampling done'
    np.save(options.output, WB)

if __name__ == "__main__":
     main()