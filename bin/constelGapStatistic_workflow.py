#!/usr/bin/env python

# python system modules
import optparse
import pickle
import sys
import os
import six

#soma_workflow
import soma_workflow.client as swf

# soma
from soma.path import find_in_path

def parseOpts(argv):
    desc = """Workflow to generate gap statistic in parallel"""
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
                      help='workflow file')
    parser.add_option('-r', '--re', dest='typer', type='str',
                      help='Type of resampling:' 
                           'b --> bootstrap' 
                           'p --> permutations' 
                           'm --> montecarlo')                 
    parser.add_option('-e', '--distance', dest='dist', type='str',
                      help='Distance:' 
                           '- sqeuclidean' 
                           '- euclidean' 
                           '- cityblock')
    parser.add_option('-f', '--wflag', dest='wflag', type='int',
                      help='Whitening of features or not ?' 
                           '2 per feature' 
                           '1 per category' 
                           '0 no')
    parser.add_option('-d', '--directory', dest='directory', metavar='FILE',
                      help= 'working directory for intermediary files')
    parser.add_option('-j', '--jobs', dest='jobs', type='int',
                      help='number of jobs')
    return parser, parser.parse_args(argv)             

def main():
    parser, (options, args) = parseOpts(sys.argv)
    
    # determines the numbers of iterations by job
    jobs = []
    npermut = options.permut / options.jobs
    nadd = options.permut % options.jobs
    nperm = [npermut] * options.jobs
    
    # run command
    for j in six.moves.xrange(options.jobs):
        # limit on the numbers of job
        if j < nadd:
            nperm[j] += 1
        outfile = os.path.join(options.directory, 'gap_job_%d.npy' % j)
        cmd = ['python', find_in_path('sillonsKmedoidsAndGapStat_map.py'),
               '-m', options.matrix, 
               '-k', str(options.kmax),
               '-n', str(options.niter),
               '-o', outfile,
               '-p', str(nperm[j]),
               '-r', str(options.typer),
               '-e', str(options.dist),
               '-f', str(options.wflag)]
        # parameters for each job    
        job = swf.Job(command=cmd, name='job_%d' % j)
        jobs.append(job)

    # create workflow   
    workflow = swf.Workflow(jobs=jobs, name='GAP')    
    pickle.dump(workflow, open(options.output, 'w'))    
    
if __name__ == "__main__":
    main()