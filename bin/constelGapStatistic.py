#!/usr/bin/env python

# python system modules
import numpy as np
import optparse
import sys

# pycluster
import Pycluster as pc

# constellation
import constel.lib.connmatrix.connmatrixtools as clcm
from constel.lib.clustering.clusterstools import wcDist

# scipy
import scipy.spatial.distance as ssd
import scipy.cluster.vq as scv

# soma
from soma import aims

def parseOpts(argv):
    desc = """(1) Cluster the observed data. 
              (2) Compute the standard deviation. 
              (3) Choose the number of clusters."""
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
                      help='Type of resampling:' 
                           'b=bootstrap' 
                           'p=permutations' 
                           'm=montecarlo ')                 
    parser.add_option('-e', '--distance', dest='dist', type='str',
                      help='Distance:' 
                           'sqeuclidean' 
                           'euclidean' 
                           'cityblock')
    parser.add_option('-f', '--wflag', dest='wflag', type='int',
                      help='Whitening of features or not ?' 
                           '2 per feature' 
                           '1 per category' 
                           '0 no')                    
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
    
    if (options.dist == 'sqeuclidean' ):
        fact = 2
    else:
        fact = 1
    
    if (options.typer == 'm'):
        print 'Generating uniform distribution with Monte-carlo'
        uniform = clcm.generate_uniform_matrix(features)
    elif (options.typer == 'b'):
        print 'Generating uniform distribution that will be bootstrapped'
        uniform = clcm.generate_uniform_matrix(features)
    elif (options.typer == 'p'):
        print 'Reference distribution will be sampled by permutations of the original one'
        uniform = features.copy()
    
    WB = np.zeros((options.permut, options.kmax + 1))
    print 'Re-sampling ', options.permut, ' times to estimate uniform Wk'
    for iteration in range(options.permut):
        if ((iteration%1) == 0):
            print '  --> Iteration ', iteration
        if (options.typer == 'm'):
            featuni = clcm.generate_uniform_matrix(features)
        elif (options.typer == 'b'):
            featuni, sampuni, suruni = clcm.bootstrap_resampling(uniform)
        elif (options.typer == 'p'):
            featuni = clcm.permutation_resampling(uniform)
            print '      Permutation done'
        
        rdist = ssd.squareform(ssd.pdist(featuni, options.dist))
        print '      Distance matrix done'
        for K in range(1, options.kmax + 1):
            uniclusterid, unierror, uninfound = pc.kmedoids(rdist, K, options.niter)
            uniwc = wcDist(rdist, uniclusterid, K, fact)
            #WB[iteration, K] = np.log(uniwc)
            WB[iteration, K] = uniwc
        print '      Clustering assessed'
              
    print 'Resampling done'
    np.save(options.output, WB)

if __name__ == "__main__":
     main()