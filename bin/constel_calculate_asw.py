#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################


#----------------------------Imports-------------------------------------------


# system module
import os
import sys
import csv
import json
import numpy
import argparse
import textwrap
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages


# soma module
from soma import aims

# constel module
try:
    from constel.lib.evaluation.indices import silhouette_index
except:
    pass


#----------------------------Functions-----------------------------------------

# Command parameters
doc = """
Average Silhouette Width
~~~~~~~~~~~~~~~~~~~~~~~~

Calculate the average silhouette width in order to find the optimal number of
clusters.
"""

def mylist(string):
    return json.loads(string)


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(doc))

    # adding arguments
    parser.add_argument(
        "matrices", type=mylist,
        help="List of the individual reduced matrix.")
    parser.add_argument(
        "kmax", type=int,
        help="The max number of clusters.")
    parser.add_argument(
        "ofile", type=str,
        help="The results will write in this PDF file.")
    parser.add_argument(
        "-s", "--scaley", type=mylist, dest="ybound",
        help="Do the scale given on the axe Y.")
    parser.add_argument(
        "-r", "--removek2", action="store_true", dest="ignore_k2",
        help="Ignore K=2 in the research of the optimal number of clusters.")
    parser.add_argument(
        "-c", "--scores", action="store_true", dest="save_scores",
        help="Write the ASW values (.npy) on the disk.")

    # parsing arguments
    return parser, parser.parse_args(argv)


def create_page(name_matrix, matrix, kmax, ybound=[0., 1.], ignore_k2=False):
    """
    """
    # Create a figure instance (ie. a new page)
    fig = pyplot.figure()

    # Add a centered title to the figure
    fig.suptitle("Matrix name: " + os.path.basename(str(name_matrix)),
                 fontsize=14,
                 color="blue",
                 style="oblique",
                 fontweight="bold",
                 family="sans-serif")

    ax1 = fig.add_axes([0.1, 0.35, 0.35, 0.5],
                       xlabel="Clusters (K)",
                       ylabel="Average Silhouette Width (ASW)",
                       color_cycle='g',
                       autoscale_on=False,
                       ybound=ybound,
                       xbound=[2, kmax])
    ax2 = fig.add_axes([0.1, 0.1, 0.7, 0.2])
    ax3 = fig.add_axes([0.55, 0.55, 0.55, 0.2])

    tvals = []
    tkeys = []
    tkeys.append("K")
    tvals.append("ASW (%)")
    dict_clusters = {}
    for idx, k in enumerate(range(2, kmax + 1)):
        asw_score, sample_silhouette, clusterid = silhouette_index(matrix.T, k)
        tkeys.append(k)
        tvals.append(round(asw_score, 4) * 100)
        dict_clusters[k] = asw_score
    print dict_clusters

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
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt) + ".\n (You "
                 " have decided to ignore Kopt=2 clusters.)",
                 color="red")
    else:
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt) + ".",
                 color="red")
    ax3.axis("off")

    return dict_clusters


#----------------------------Main program--------------------------------------


def main():
    """
    """
    if len(sys.argv[:]) == 6:
        # load the arguments of parser (delete script name: sys.arg[0])
        arguments = (json.dumps(eval(sys.argv[1])),
                     sys.argv[2],
                     sys.argv[3],
                     sys.argv[4],
                     sys.argv[5])
    elif len(sys.argv[:]) == 5:
        # load the arguments of parser (delete script name: sys.arg[0])
        arguments = (json.dumps(eval(sys.argv[1])),
                     sys.argv[2],
                     sys.argv[3],
                     sys.argv[4])
    else:
        # not autoscale Y
        arguments = (json.dumps(eval(sys.argv[1])),
                     sys.argv[2],
                     sys.argv[3])
    parser, args = parse_args(arguments)

    # the PDF document
    pp = PdfPages(args.ofile)

    for matrix in args.matrices:
        red_mat = aims.read(matrix)
        rmat = numpy.asarray(red_mat)[:, :, 0, 0]
        if rmat.shape[0] > rmat.shape[1]:
            rmat = rmat.T
        dict_clusters = create_page(matrix,
                                    rmat,
                                    args.kmax,
                                    args.ybound,
                                    args.ignore_k2)
        # Done with the page
        pp.savefig()
        if args.save_scores:
            matrix_name = '.'.join(matrix.split('.')[:-1])
            asw_name = matrix_name + "_asw.csv"
            with open(asw_name, 'wb') as csv_file:
                writer = csv.writer(csv_file)
                for key, value in dict_clusters.items():
                    writer.writerow([key, value])

    # Write the PDF document to the disk
    pp.close()


if __name__ == "__main__":
    main()
