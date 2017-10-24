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

# matplotlib. Set it to a non-graphical PDF backend
import matplotlib
matplotlib.use("pdf")
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

#def mylist(string):
    #return json.loads(string)


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(doc))

    # adding arguments
    parser.add_argument(
        "matrices", type=str, nargs="+",
        help="List of the individual reduced matrix.")
    parser.add_argument(
        "kmax", type=int,
        help="The max number of clusters.")
    parser.add_argument(
        "ofile", type=str,
        help="The results will write in this PDF file.")
    parser.add_argument(
        "-s", "--scaley", type=float, action="append", dest="ybound",
        default=[0., 1.],
        help="Do the scale given on the axe Y.")
    parser.add_argument(
        "-r", "--removek2", action="store_true", dest="ignore_k2",
        help="Ignore K=2 in the research of the optimal number of clusters. "
        "Equivalent to --kmintype min --kmin 3.")
    parser.add_argument(
        "-c", "--scores", action="store_true", dest="save_scores",
        help="Write the ASW values (.npy) on the disk.")
    parser.add_argument(
        "-m", "--kmintype", dest="kmin_type", default="",
        help="Selection method for minimum K. Possible values are "
        "\"none\" (default), \"min\": hard-coded min "
        "(see --kmin parameter), "
        "\"avg_pts\": corresponds to an average number of points "
        "(kmin = region pts / max_avg_pts, see --max_avg_pts)")
    parser.add_argument(
        "--kmin", type=int, default=0,
        help="minimum number of clusters taken into account to search for the "
        "optimum. Used when option --kmintype is \"min\".")
    parser.add_argument(
        "--max_avg_pts", type=float,
        help="Max. average number of points for min number of "
        "clusters. Used when opton --kmintype is \"avg_pts\". "
        "Ideally it should rather be a surface but we would require "
        "to input the mesh as additional argument for that to get "
        "the surface, and this should be enough.")

    # parsing arguments
    args = parser.parse_args(argv)
    if args.ignore_k2:
        if args.kmin_type not in ("", "min") or args.kmin not in (0, 3):
            raise ValueError(
                "Cannot specify --removek2 option with other --kmintype "
                "options")
        args.kmin_type = "min"
        args.kmin = 3

    print(args)

    return parser, args


def create_page(name_matrix, matrix, kmax, ybound=[0., 1.],
                kmin_type="none", kmin=0, max_avg_pts=0):
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

    kmin_v = 2
    if kmin_type == "min" and kmin >= 2:
        kmin_v = kmin
    elif kmin_type == "max_avg_pts" and max_avg_pts > 0:
        kmin_v = matrix.shape[0] / max_avg_pts
    print('kmin:', kmin_v)

    if kopt < kmin_v:
        del dict_clusters[kopt]
        k_opt = max([v for k, v in six.iteritems(dict_clusters)
                     if v >= kmin_v])
        kopt = dict_clusters.keys()[
            dict_clusters.values().index(k_opt)]
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt)
                 + ".\n (You have decided to ignore Kopt < %d "
                 "clusters.)" % kmin_v,
                 color="red")
    else:
        ax3.text(0., 0.,
                 "The optimal number of clusters is " + str(kopt)
                 + ".",
                 color="red")
    ax3.axis("off")

    return dict_clusters


#----------------------------Main program--------------------------------------


def main():
    """
    """
    # load the arguments of parser (delete script name: sys.arg[0])
    #if len(sys.argv) < 2:
      #parse_args(("-h", ))
      #sys.exit(0)

    arguments = sys.argv[1:]
    #if sys.argv[1].startswith("["):
        #arg0 = json.dumps(eval(sys.argv[1]))
    #else:
        #arg0 = sys.argv[1]
    #arguments = (arg0, ) + tuple(sys.argv[2:])
    #if len(sys.argv[:]) == 6:
        #arguments = (json.dumps(eval(sys.argv[1])),
                     #sys.argv[2],
                     #sys.argv[3],
                     #sys.argv[4],
                     #sys.argv[5])
    #elif len(sys.argv[:]) == 5:
        ## load the arguments of parser (delete script name: sys.arg[0])
        #arguments = (json.dumps(eval(sys.argv[1])),
                     #sys.argv[2],
                     #sys.argv[3],
                     #sys.argv[4])
    #else:
        ## not autoscale Y
        #arguments = (json.dumps(eval(sys.argv[1])),
                     #sys.argv[2],
                     #sys.argv[3])
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
                                    args.kmin_type,
                                    args.kmin,
                                    args.max_avg_pts)
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
