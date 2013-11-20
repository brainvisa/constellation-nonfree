#!/usr/bin/env python
from optparse import OptionParser
import soma_workflow.client as swf
import os, sys
import pickle
import math

parser = OptionParser( 'blabla' )
parser.add_option( '-m', '--matrix', dest = 'matrix', metavar = 'FILE',
                    help = 'input connectivity matrix')
parser.add_option( '-k', '--kmax', dest = 'kmax', type='int',
                    help = 'K max')
parser.add_option( '-n', '--niter', dest = 'niter', type='int',
                    help = 'number of iterations of same kmedoid')
parser.add_option( '-p', '--permutations', dest = 'permut', type='int',
                    help = 'number of permutations')
parser.add_option( '-o', '--output', dest = 'output', metavar = 'FILE',
                    help = 'workflow file')
parser.add_option( '-d', '--directory', dest = 'directory', metavar = 'FILE',
                    help = 'working directory for intermediary files')
parser.add_option( '-j', '--jobs', dest = 'jobs', type='int',
                    help = 'number of jobs')
options, args = parser.parse_args(sys.argv)

jobs = []
npermut = int( math.ceil( float( options.permut ) / options.jobs ) )
nperm = [npermut] * (options.jobs-1)
nperm.append( options.permut - (options.jobs-1) * npermut )
for j in xrange( options.jobs ):
  outfile = os.path.join( options.directory, 'gap_job_%d.npy' % j )
  cmd = [ '/volatile/sandrine/svn/brainvisa/connectomist/constellation-private/trunk/bin/sillonsKmedoidsAndGapStat_map.py', '-m', options.matrix, '-k', str(options.kmax),
    '-n', str(options.niter), '-o', outfile, '-p', str(nperm[j]) ]
  job = swf.Job( command=cmd, name='job_%d' % j )
  jobs.append( job )
  
workflow = swf.Workflow( jobs=jobs, name='GAP' )

pickle.dump( workflow, open( options.output, 'w' ) )

  
