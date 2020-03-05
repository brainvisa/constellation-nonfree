#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

"""
This script does the following:
*

Main dependencies:

Author: Sandrine Lefranc, 2015
"""


#----------------------------Imports-------------------------------------------


# python system modules
from __future__ import print_function
from __future__ import absolute_import
import os
import sys
import numpy
import argparse
import textwrap

# scipy module
from scipy import stats

# soma module
from soma import aims

# pycluster module
import Pycluster as pc

# constel module
import constel.lib.evaluation.indices as cv
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
from six.moves import range


#----------------------------Functions-----------------------------------------


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Cluster criterion: David-Bouldin Index, Dunn Index and
            Calinski-Harabasz Index.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        'matrix', type=str,
        help='Reduced connectivity matrix (Ndim, Nsample)')
    parser.add_argument(
        "cortical_region", type=str,
        help="The study cortical region.")
    parser.add_argument(
        'nbiter', type=int,
        help='Number of iterations of same kmedoid')
    parser.add_argument(
        'kmax', type=int,
        help='K max for the clustering')
    parser.add_argument(
        "ofile", type=str,
        help='Validity indexes file')

    # parsing arguments
    return parser, parser.parse_args(argv)


def create_page(title, measures, name_y, cortical_region, ybound=[0., 1.],
                ignore_k2=False):
    """
    """
    # Create a figure instance (ie. a new page)
    fig = pyplot.figure()

    # Add a centered title to the figure
    fig.suptitle(title,
                 fontsize=14,
                 fontweight="bold",
                 color="blue")

    ax1 = fig.add_axes([0.1, 0.35, 0.35, 0.5],
                       xlabel="Clusters (K)",
                       ylabel=name_y,
                       title=cortical_region,
                       color_cycle='g',
                       autoscale_on=False,
                       ybound=ybound,
                       xbound=[2, len(measures) + 1])
    ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.2])
    ax3 = fig.add_axes([0.55, 0.55, 0.55, 0.2])

    tvals = []
    tkeys = []
    tkeys.append("K")
    tvals.append(name_y + " (%)")
    dict_clusters = {}
    for idx, k in enumerate(range(len(measures))):
        tkeys.append(idx + 2)
        tvals.append(round(measures[idx], 4) * 100)
        dict_clusters[idx + 2] = measures[idx]

    table = []
    table.append(tkeys)
    table.append(tvals)

    # give the larger ASW
    k_opt = max(dict_clusters.values())

    # search the id of the optimal number of clusters
    kopt = list(dict_clusters.keys())[list(dict_clusters.values()).index(k_opt)]

    ax1.plot(list(dict_clusters.keys()), list(dict_clusters.values()), "k",
             color="red",
             linewidth=3)

    # put a grid on the curve
    ax1.grid(True)

    ax2.table(cellText=table, loc="center")
    ax2.axis("off")

    if ignore_k2 and kopt == 2:
        del dict_clusters[kopt]
        k_opt = max(dict_clusters.values())
        kopt = list(dict_clusters.keys())[list(dict_clusters.values()).index(k_opt)]
        ax1.annotate("kopt",
                     xy=(kopt, k_opt),
                     xytext=(kopt + 0.5, k_opt + 0.05),
                     arrowprops=dict(facecolor='black', shrink=0.05),)
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt) + ".\n (You"
                 " have decided to ignore Kopt=2 clusters.)",
                 color="red")
    else:
        ax1.plot(list(dict_clusters.keys()), list(dict_clusters.values()), "k",
                 color="red",
                 linewidth=3)
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt) + ".",
                 color="red")

    ax3.axis("off")


def main():
    """Execute the command to create and write ...
    """
    # Load the arguments of parser (delete script name: sys.arg[0])
    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    # load the matrix file
    matrix = aims.read(args.matrix)

    # transpose the matrix
    feat = numpy.asarray(matrix)[:, :, 0, 0].T

    Nsample = feat.shape[0]
    #Ndim = feat.shape[1]
    #print('Number of samples: ', Nsample, 'and Dimension: ', Ndim)

    nbIter = args.nbiter
    kmax = args.kmax

    distance1 = pc.distancematrix(feat, dist='e')

    print('Euclidean distance at power 2')
    distMat1 = [numpy.array([])]
    for i in range(1, feat.shape[0]):
        distMat1.append((distance1[i] ** 2).copy())

    print('Converting distance matrix')
    distMat = numpy.zeros((Nsample, Nsample))
    for i in range(1, Nsample):
        for j in range(0, i):
            distMat[i, j] = distMat1[i][j]
            distMat[j, i] = distMat1[i][j]
    Nsample = distMat.shape[0]
    print('There are', Nsample, 'samples')

    sDunn = numpy.zeros(kmax + 1)
    sDB = numpy.zeros(kmax + 1)
    sCH = numpy.zeros(kmax + 1)

    for k in range(kmax):
        print('Runnning K-medoids with K = ', k)
        uniclusterid, unierror, uninfound = pc.kmedoids(distMat, k + 1, nbIter)
        print('getting validity indexes')
        sDB[k + 1], sDunn[k + 1], sCH[k + 1] = cv.compute_cluster_validity(
            distMat, uniclusterid, k + 1)
        print('sDB = ', sDB[k + 1])
        print('sDunn = ', sDunn[k + 1])
        print('sCH = ', sCH[k + 1])

    sDB[1] = sDB.min()
    sDB[0] = sDB.min()
    sDunn[1] = sDunn.min()
    sDunn[0] = sDunn.min()
    sCH[1] = sCH.min()
    sCH[0] = sCH.min()

    print('Normalizing scores')
    sDB = (sDB - sDB.min()) / float(sDB.max() - sDB.min())
    sDunn = (sDunn - sDunn.min()) / float(sDunn.max() - sDunn.min())
    sCH = (sCH - sCH.min()) / float(sCH.max() - sCH.min())

    # scores = numpy.vstack((sDB, sCH))
    scores = numpy.vstack((sDunn, sDB))
    scores = numpy.vstack((scores, sCH))

    SFA = numpy.mean(scores, axis=0)
    SFG = stats.mstats.gmean(scores, axis=0)
    SFMed = numpy.median(scores, axis=0)

    all_measures = []
    for element in (sDB, sDunn, sCH, SFA, SFG, SFMed):
        all_measures.append(element)

    scores = ["Davies Bouldin Index", "Dunn Index", "Calinski Harabasz Index",
              "SFA", "SFG", "SFMed"]

    # the PDF document
    pp = PdfPages(args.ofile)

    for idx, measures in enumerate(all_measures):
        title_page = scores[idx] + ": \n" + os.path.basename(args.matrix)
        create_page(title_page,
                    measures,
                    scores[idx],
                    args.cortical_region)
        # Done with the page
        pp.savefig()

    # Write the PDF document to the disk
    pp.close()


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
