#!/usr/bin/env python

# python system modules
import numpy as np
import optparse
import pickle
import sys
import six

# soma_workflow
import soma_workflow.client as swf
import nibabel as nib
# soma
from soma.path import find_in_path

def parseOpts(argv):
    desc = """Workflow to generate the distance matrix in parallel"""
    parser = optparse.OptionParser(desc)
    parser.add_option('-m', '--matrix', dest='matrix',
                      help='input concatenated matrix')                
    parser.add_option('-j', '--jobs', dest='jobs', type='int',
                      help='number of jobs')
    parser.add_option('-i', '--input', dest='indir', metavar='FILE',
                      help='file to construct the distance matrix')
    parser.add_option('-o', '--output', dest='output', metavar='FILE',
                      help='workflow file')
    parser.add_option('-d', '--method', dest='method', default='scipy',
                      help='two methods to generate the distance matrix:'
                           '(1) scipy --> risk of Memory Error,' 
                           '(2) square --> partition of distance matrix,'
                           '[default: %default]'
                           '(3) vector --> to use kmedoid algorithm')
    return parser, parser.parse_args(argv)             

def main():
    parser, (options, args) = parseOpts(sys.argv)
    
    # matrix loading
    img = nib.load()
    matrix = img.get_data()
    #matrix = np.load(options.matrix)
    m, n = matrix.shape
    
    # create the file to receive the elements of the distance matrix
    if options.method == 'square':
        # square distance matrix
        f = open(options.indir, 'w')
        f.seek((m * m * 8) - 1)
        f.write('x')
        f.close()
    else:
        # distance matrix as vector
        f = open(options.indir, 'w')
        f.seek(m * (m - 1) / 2 * 8 - 1)
        f.write('x')
        f.close()

    # determines the numbers of iterations by job
    jobs = []
    nb_iteration = m / options.jobs
    nb_add = m % options.jobs
    list_iteration = [nb_iteration] * options.jobs

    # hacks...
    starting_iteration = 0
    k = starting_iteration
    x = m - 1
    y = 1

    # run command
    for j in six.moves.xrange(options.jobs):
        # limit on the numbers of job
        if j < nb_add:
            list_iteration[j] += 1
        cmd = ['python', find_in_path('constelSquareDistanceMatrix.py'),
               '-m', options.matrix,
               '-s', str(starting_iteration),
               '-t', str(list_iteration[j]),
               '-k', str(k),
               '-i', options.indir,
               '-d', options.method]
        # parameters for each job    
        job = swf.Job(command=cmd, name='job_%d' % j)
        jobs.append(job)
        
        # definition of starting iteration, corresponding to matrix line
        starting_iteration += list_iteration[j]
        
        # indice of k in triangular distance matrix as vector (hack...)
        for i in six.moves.xrange(list_iteration[j]):
            k += x
            x = (m - 1) - y
            y += 1

    # create workflow
    workflow = swf.Workflow(jobs=jobs, name='square_distance_matrix')
    pickle.dump(workflow, open(options.output, 'w'))
    
if __name__ == "__main__":
    main()
