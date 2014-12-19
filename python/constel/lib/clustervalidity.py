#!/usr/bin/env python

import numpy
import Pycluster as pc
from scipy.spatial.distance import pdist, squareform
from constel.lib.clustering import clusterstools


def silhouette_sample(X, K):
    n = X.shape[0]
    dist = pc.distancematrix(X, dist='e')
    NiterK = 100
    clusterid, error, nfound = pc.kmedoids(dist, K, NiterK)
    labels = 0
    for i in numpy.unique(clusterid):
        clusterid[clusterid == i] = labels
        labels += 1
    cl = len(numpy.unique(clusterid))
    distance = squareform(pdist(X))
    kIndices = [numpy.flatnonzero(clusterid == k) for k in range(cl)]

    a = numpy.zeros(n)
    b = numpy.zeros(n)
    for i in range(n):
        a[i] = numpy.mean(
            [distance[i][ind] for ind in kIndices[clusterid[i]] if ind != i])
        b[i] = numpy.min(
            [numpy.mean(distance[i][ind]) for k,
             ind in enumerate(kIndices) if clusterid[i] != k])
    s = (b - a) / numpy.maximum(a, b)
    print s
    # p.silhouette_plot(X, s, clusterid)
    return s


def silhouette_score(X, K):
    s = silhouette_sample(X, K)
    score = numpy.mean(s)
    return score


def computeClusterValidity(distance, clusterid, K):
    centers = clusterstools.getCenters(clusterid, K)
    distC = clusterstools.distToCenters(distance, centers, clusterid, K)
    sDB = 0
    sCH = 0
    centAll = clusterstools.centroid(distance)

    Nsample = distance.shape[0]

    print '     step1'
    for i in range(K):
        Di = numpy.zeros(K)
        for j in range(K):
            if (i != j):
                Di[j] = (distC[i] + distC[j]) / distance[
                    centers[i], centers[j]]
        sDB += Di.max()
    sDB = sDB / float(K)

    print '     step2'
    if (K > 2):
        dinter = numpy.array([])
        dintra = numpy.array([])
        for i in range(Nsample):
            for j in range(i + 1, Nsample):
                if (clusterid[i] != clusterid[j]):
                    dinter = numpy.append(dinter, distance[i, j])
                else:
                    dintra = numpy.append(dintra, distance[i, j])
        print dinter.shape, dintra.shape
        sDunn = dinter.min() / dintra.max()
    else:
        sDunn = 0

    print '     step3'
    if (K > 1):
        B = 0
        W = 0
        for i in range(K):
            B += (distance[centers[i], centAll] * numpy.where(
                clusterid == centers[i])[0].size)
            for j in range(Nsample):
                if (clusterid[j] == centers[i]):
                    W += distance[j, centers[i]]
        sCH = B * (Nsample - K) / (W * 10.0 * (K - 1))
    else:
        sCH = 0
    print '     OK'
    return (1.0 - sDB), sDunn, sCH