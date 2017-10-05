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

Main dependencies: PyAims library

Author: Sandrine Lefranc, 2015
"""


#----------------------------Imports-------------------------------------------


# system module
import os
import sys
import json
import argparse
import textwrap
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

# soma module
from soma import aims

# constel module
try:
    from constel.lib.clustering.clusterstools import entropy
    import constel.lib.evaluation.indices as measure
    from constel.lib.utils.misctools import sameNbElements
except:
    pass


#----------------------------Functions-----------------------------------------


def mylist(string):
    return json.loads(string)


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Calculate the differents scores from time textures in order to find
            the optimal number of clusters.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        "cortical_parcellation_1", type=str,
        help="Cortical parcellation resulting from clustering algorithm.")
    parser.add_argument(
        "cortical_parcellation_2", type=str,
        help="Cortical parcellation resulting from clustering algorithm.")
    parser.add_argument(
        "timestep_max", type=int,
        help="The max number of clusters.")
    parser.add_argument(
        "outputdir", type=str,
        help="The results will write in this PDF file.")
    parser.add_argument(
        "title", type=str,
        help="The title of the pdf file.")
    parser.add_argument(
        "cortical_region", type=str,
        help="The study cortical region.")
    parser.add_argument(
        "-s", "--scaley", type=mylist, dest="ybound",
        help="Do the scale given on the axe Y.")
    parser.add_argument(
        "-r", "--removek2", action="store_true", dest="ignore_k2",
        help="Ignore K=2 in the research of the optimal number of clusters.")

    # parsing arguments
    return parser, parser.parse_args(argv)


def mkdir_path(path):
    if not os.access(path, os.F_OK):
        os.makedirs(path)


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
    kopt = dict_clusters.keys()[dict_clusters.values().index(k_opt)]

    ax1.plot(dict_clusters.keys(), dict_clusters.values(), "k",
             color="red",
             linewidth=3)

    # put a grid on the curve
    ax1.grid(True)

    ax2.table(cellText=table, loc="center")
    ax2.axis("off")

    if ignore_k2 and kopt == 2:
        del dict_clusters[kopt]
        k_opt = max(dict_clusters.values())
        kopt = dict_clusters.keys()[dict_clusters.values().index(k_opt)]
        ax1.annotate("kopt",
                     xy=(kopt, k_opt),
                     xytext=(kopt + 0.5, k_opt + 0.05),
                     arrowprops=dict(facecolor='black', shrink=0.05),)
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt) + ".\n (You"
                 " have decided to ignore Kopt=2 clusters.)",
                 color="red")
    else:
        ax1.plot(dict_clusters.keys(), dict_clusters.values(), "k",
                 color="red",
                 linewidth=3)
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt) + ".",
                 color="red")

    ax3.axis("off")


#----------------------------Main program--------------------------------------


def main():
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    # directory for all validation results organized by gyrus
    odir = os.path.join(args.outputdir, "validation")
    mkdir_path(odir)

    # directory for all measures
    ofile = os.path.join(odir, args.title)

    # read clustering texture (aims)
    c1 = aims.read(args.cortical_parcellation_1)
    c2 = aims.read(args.cortical_parcellation_2)

    # calculate measure between two clustering (1 and 2)
    randindex_values = []
    cramer_values = []
    mutualinformation_values = []
    homogeneity = []
    completeness = []
    v_measure = []
    all_measures = []
    for i in range(2, args.timestep_max + 1):
        # list of labels (1 and 2)
        labels1 = c1[i].arraydata()
        labels2 = c2[i].arraydata()

        # keep the same number of elements
        l1, l2 = sameNbElements(labels1, labels2)

        # measure value for K
        ri = measure.rand_index(l1, l2)
        cv = measure.cramer_v(l1, l2)
        mi = measure.mutual_information(l1, l2)

        # calculate entropy for each clustering
        entropy1 = entropy(l1)
        entropy2 = entropy(l2)

        # calculate homogeneity and completeness
        homog = mi / (entropy1)
        compl = mi / (entropy2)

        # calculate V measure: equivalent to the mutual information
        if homog + compl == 0.0:
            v = 0.0
        else:
            v = (2.0 * homog * compl / (homog + compl))

        # list of measure values for K = 2:kmax
        randindex_values.append(ri)
        cramer_values.append(cv)
        mutualinformation_values.append(mi)
        homogeneity.append(homog)
        completeness.append(compl)
        v_measure.append(v)

    all_measures.append(randindex_values)
    all_measures.append(cramer_values)
    all_measures.append(mutualinformation_values)
    all_measures.append(homogeneity)
    all_measures.append(completeness)
    all_measures.append(v_measure)

    scores = ["Rand Index", "Cramer_V", "Mutual Information",
              "Homogeneity", "Completeness", "V_measure"]

    # the PDF document
    pp = PdfPages(ofile + ".pdf")

    for idx, measures in enumerate(all_measures):
        title_page = scores[idx] + ": " + args.title
        create_page(title_page,
                    measures,
                    scores[idx],
                    args.cortical_region,
                    args.ybound,
                    args.ignore_k2)
        # Done with the page
        pp.savefig()

    # Write the PDF document to the disk
    pp.close()


if __name__ == "__main__":
    main()
