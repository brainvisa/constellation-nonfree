#!/usr/bin/env python
# Filename: sillonsKMedoids.py

from numpy import *
from scipy import cluster as cl
import scipy
from soma import aims
import anatomist.direct.api as anatomist
import Pycluster as pc
import sys
import math
from optparse import OptionParser

#featFile='/Users/olivier/dataLocal/SillonsVlad/S.F.inf./Left/spamGlobalLocalMontreal/SFinfFeatures.txt'
#featFile='/Users/olivier/dataLocal/SillonsVlad/S.F.inf./Left/spamGlobalLocalMontreal/uniformFeatures.txt'


distK=[] #classe K, les distances.
indexK=[] # classe K, les index
centerK=[]
subjK=[] # classe K, les sujets
subjCentK=[] # classe K, le sujet le plus proche du centre

errE_K=[] #error as a function of K, this is plotted to find the optimal K
errRmean_K=[] # the same divided by the number of clusters
Nresamp=2000 # How many resmapling during validation

plot=0 # plot de l'erreur globale vs nbre de clusters


####################################################################
# normalize the data (sigma=1) per category
####################################################################

def partialWhiten(features):
     dist=features[:,0:6]
     direc=features[:,6:12]
     sdist=std(dist)
     sidrec=std(direc)
     direc=direc/sidrec
     dist=dist/sdist
     white=hstack((dist, direc))
     return white

####################################################################
# generate samples from uniform distrib with same shape than
# original samples
####################################################################

def GenerateUniform(whiteFeat, Nsample, Ndim):     
     n0=(whiteFeat[:,0].max()-whiteFeat[:,0].min())*random.random_sample((Nsample,1))+whiteFeat[:,0].min()
     n1=(whiteFeat[:,1].max()-whiteFeat[:,1].min())*random.random_sample((Nsample,1))+whiteFeat[:,1].min()
     s=hstack((n0,n1))
     
     for i in range(2, Ndim):
          ni=(whiteFeat[:,i].max()-whiteFeat[:,i].min())*random.random_sample((Nsample,1))+whiteFeat[:,i].min()
          s=hstack((s, ni))
          
     print 'Generated uniform matrix of shape: ', s.shape
     return s
     
     
####################################################################
# case resampling, i.e. bootstrap with replacement
####################################################################
  
def caseResampling(feat):
     N=feat.shape[0]
     bootVect=random.random_integers(N, size=N)
     bootVect=bootVect-1
     newFeat=feat[bootVect[0],:]
     for i in range(1,bootVect.size):
          newFeat=vstack((newFeat, feat[bootVect[i],:]))
     survive=array([])
     for i in bootVect:
          if ((where(survive==i)[0]).size==0):
               survive=append(survive, array([int(i)]))
     return newFeat, bootVect, survive
     # newFeat est la nouvelle matrice de features
     # bootVect est la liste des indices du reechantillonage, 
     # survive est la liste des elements ayant survecu au reechantillonage

####################################################################
# resampling by permutation of the features 
####################################################################
  
def permutationResampling(feat):
     Nsamples=feat.shape[0]
     Nfeat=feat.shape[1]
     
     permF=random.permutation(feat[:,0].reshape((Nsamples,1)))
     for f in range (1, Nfeat):
          permF=hstack((permF, random.permutation(feat[:,f].reshape((Nsamples,1))) ))

     return permF     
     
####################################################################
# within cluster sum of distance 
####################################################################

def getCenters(clust,K):
     i=0
     cent=array([])
     #print clust.size
     while (cent.size < K):
          c=clust[i]
          if ((where(cent==c)[0]).size==0):
               cent=append(cent,  array([int( c )]))
          i+=1
     return cent

def wcDist(distance, clusterid, K, fact):
     centers=getCenters(clusterid, K)
     total=0
     W=zeros(K)
     Total=0.0
     for l in range(K):
          sum=0.0
          indit=where(clusterid==centers[l])[0] #le tuple des indexs
          for i in indit:
               for j in indit:
                    if (i<j):
                         sum+=(distance[j][i]**fact)
                    elif (j<i):
                         sum+=(distance[i][j]**fact)
          W[l]=sum/float(2.0*indit.size)
          Total+=W[l]
     return Total

####################################################################
# Main function
# usage : python silllonsKMedoidsAndGapStat.py 2 se p /Users/olivier/dataLocal/SillonsVlad/S.F.inf./Left/spamGlobalLocalMontreal/SFinfFeatures.txt expeDir/
####################################################################

def main(arguments):
    parser = OptionParser( 'blabla' )
    parser.add_option( '-i', '--input', dest = 'input', metavar = 'FILE', action='append',
                        help = 'input permutation files (list)')
    parser.add_option( '-m', '--matrix', dest = 'matrix', metavar = 'FILE',
                        help = 'input connectivity matrix')
    parser.add_option( '-k', '--kmax', dest = 'kmax', type='int',
                        help = 'K max')
    parser.add_option( '-n', '--niter', dest = 'niter', type='int',
                        help = 'number of iterations of same kmedoid')
    parser.add_option( '-o', '--output', dest = 'output', metavar = 'FILE',
                        help = 'output gap file')
    options, args = parser.parse_args(arguments)
    
    wflag=0 # whitening of features or not ? (2 per feature, 1 per category, 0 no)
    d='e' # distance : 'se'=squared euclidean, 'e'=euclidean, 'b'=city-block
    expeDir=arguments[5] 
    
    gapFile=options.output

    featFile=options.matrix
    
    NiterK=options.niter # nb iteration for k-medoids algorithm
    #print 'WARNING VERY SLOW - NiterK has been set to 10000 !!! Go change it if that is not what you want'

    Kmax=options.kmax # max nb of clusters
    vects=[]
    print 'Opening and reading ', featFile
    
    vects=aims.read(featFile)
    
    feat=asarray(vects)[:,:,0,0].transpose()
    
    Nsample=feat.shape[0]
    Ndim=feat.shape[1]
    
    print 'Number of samples: ', Nsample
    print 'Dimension: ', Ndim
    
    
    if d=='se':
        fact=2
        d='e'
    else:
        fact=1

    #whiteFeat=cl.vq.whiten(feat)
    if (wflag==2):
        print 'Per feature whitening'
        whiteFeat=scipy.cluster.vq.whiten(feat)
    elif (wflag==1):
        print 'Per category (dist and dir) whitening'
        whiteFeat=partialWhiten(feat)
    else:
        print 'No whitening'
        whiteFeat=feat
#     centers, disto=cl.vq.kmeans(whiteFeat, K, 100)
#     code, dist = cl.vq.vq(whiteFeat, centers)

    distance1=pc.distancematrix(whiteFeat, dist=d)
    
    # addition I should have made before : squaring the euclidean distances.
    
    distance=[array([])]
    
    if (d=='e'):
        print 'Euclidean distance at power ', fact
        for i in range(1, whiteFeat.shape[0]):
              distance.append((distance1[i]**fact).copy())
    elif (d=='b'):
        print 'City-block distance'
        distance=distance1

    W=zeros(Kmax+1)
    print 'computing Wk for data'
    for K in range(1, Kmax+1):
        clusterid, error, nfound = pc.kmedoids(distance, K, NiterK)
        wc=wcDist(distance1, clusterid, K, 2)
        print 'K=', K, ', error=', error, ', wcDist=', wc
        W[K]=math.log(wc)
    print 'OK'
    print 'Wk=', W
    
    if (re=='m'):
        print 'Generating uniform distribution with Monte-carlo'
        uniform=GenerateUniform(whiteFeat, Nsample, Ndim)
    elif (re=='b'):
        print 'Generating uniform distribution that will be bootstrapped'
        uniform=GenerateUniform(whiteFeat, Nsample, Ndim)
    elif (re=='p'):
        print 'Reference distribution will be sampled by permutations of the original one'
        uniform=whiteFeat.copy()
    
    # WB=zeros((B, Kmax+1))
    
    WBparts = []
    for infile in options.input:
      WBparts.append( load( infile ) )
    WB = vstack( WBparts )
    del WBparts
              
    print 'Computing Gap'
    meanWB=mean(WB, axis=0)
    sdWB=std(WB, axis=0)
    gap=zeros(Kmax+1)
    sW=zeros(Kmax+1)
    for K in range(1, Kmax+1):
        gap[K]=(meanWB[K] - W[K])
        sW[K]=sdWB[K]*math.sqrt(1+1.0/float(B))
    
    print 'MeanWB=', meanWB
    print 'sdWB=', sdWB 
    print 'Gap=', gap
    print 'STD=', sW
    
    print 'Writing results to ', gapFile

    fileR=open(gapFile, 'w')
    fileR.write(featFile + '\n')
    lineW=''
    lineMW=''
    lineG=''
    lineS=''
    for K in range(1,Kmax+1):
        lineW=lineW + str(W[K]) + ' '
        lineMW=lineMW + str(meanWB[K]) + ' '
        lineG=lineG + str(gap[K]) + ' '
        lineS=lineS + str(sW[K]) + ' '
    fileR.write(lineW + '\n')
    fileR.write(lineMW + '\n')
    fileR.write(lineG + '\n')
    fileR.write(lineS + '\n')

    fileR.close()

if __name__ == "__main__":
     main(sys.argv)