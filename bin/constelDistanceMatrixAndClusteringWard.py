#!/usr/bin/env python

# python system modules
import numpy as np
import optparse
import pickle
import struct
import sys

# scipy
import scipy.cluster._hierarchy_wrap as hier_wrap

# soma_workflow
import soma_workflow.client as swf

def parseOpts(argv):
    desc = """Workflow to generate the clustering Ward"""
    parser = optparse.OptionParser(desc)
    
    parser.add_option('-m', '--matrix', dest='matrix',
                      help='input reduced connectivity matrix')
    parser.add_option('-o', '--output', dest='output', metavar='FILE',
                      help='workflow file')
                      
    return parser, parser.parse_args(argv)             

def main():
    parser, (options, args) = parseOpts(sys.argv)
    

    matrix = np.load(options.matrix)
    f = open("/neurospin/tmp/sandrine/blop", "r")
    n, m = matrix.shape
    a = np.zeros((n*(n-1)/2,))
    k = 100000
    nb_iteration = a.shape[0] / k
    nb_add = a.shape[0] % k
    list_iteration = [nb_iteration] * k
    idx = 0
    for j in xrange(k):
        if j<nb_add:
            list_iteration[j] += 1
        x = f.read(list_iteration[j]*8)
        a.ravel()[idx:idx + list_iteration[j]] = struct.unpack('d'*list_iteration[j], x)

        idx += list_iteration[j]
    f.close()

    # clustering: linkage ward
    Z = np.zeros((n-1,4))

    hier_wrap.linkage_euclid_wrap(a, Z, matrix, m, n, 5)   
    
    workflow = swf.Workflow(jobs=1, name='clustering_ward')
    pickle.dump(workflow, open(options.output, 'w'))
    
if __name__ == "__main__":
    main()