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
import numpy
import argparse
import textwrap
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

# System modules
import matplotlib.cm as cm

# soma module
from soma import aims

# constel module
try:
    from constel.lib.evaluation.indices import silhouette_index
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
            Calculate the average silhouette width in order to find the optimal
            number of clusters.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        "matrix", type=str,
        help="Individual reduced matrix.")
    parser.add_argument(
        "kmax", type=int,
        help="The max number of clusters.")
    parser.add_argument(
        "ofile", type=str,
        help="The results will write in this PDF file.")

    # parsing arguments
    return parser, parser.parse_args(argv)


def create_page(matrix, cluster_number):
    """
    """
    # Create a figure instance (ie. a new page)
    fig, (ax1, ax2, ax3) = pyplot.subplots(1, 3)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(matrix) + (cluster_number + 1) * 10])

    silhouette_avg, sample_silhouette_values, clusterid = silhouette_index(
        matrix, cluster_number)
    
    y_lower = 10
    for i in range(cluster_number):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[clusterid == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.spectral(float(i) / cluster_number)
        ax1.fill_betweenx(numpy.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.spectral(clusterid.astype(float) / cluster_number)
    ax2.scatter(matrix[:, 0], matrix[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors)

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")

    #ax3.plot(dict_clusters.keys(), dict_clusters.values(), "k",
     #        color="red",
      #       linewidth=3)

    ax3.set_title("")
    ax3.set_xlabel("")
    ax3.set_ylabel("")

    pyplot.suptitle(
        ("Silhouette analysis for Kmedoids clustering on sample data "
         "with cluster_number = %d" % cluster_number),
        fontsize=14, fontweight='bold')

#----------------------------Main program--------------------------------------


def main():
    """
    """
    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    # the PDF document
    pp = PdfPages(args.ofile)

    for k in range(args.kmax):
        matrix = aims.read(args.matrix)
        mat = numpy.array(matrix)[:, :, 0, 0]
        if mat.shape[0] < mat.shape[1]:
            mat = numpy.transpose(mat)
        create_page(mat, k + 2)
        # Done with the page
        pp.savefig()

    # Write the PDF document to the disk
    pp.close()


if __name__ == "__main__":
    main()
