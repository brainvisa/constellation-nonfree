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
"""


#----------------------------Imports-------------------------------------------


# python system module
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import json
import glob
import numpy
import argparse
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# aims module
from soma import aims

# constel module
from constel.lib.utils.filetools import select_ROI_number


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
            Calculate a connectome from a nomenclature and create its snapshot.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument("cortical_parcellation",
                        type=str,
                        help="")
    parser.add_argument("freesurfer_cortical_parcellation",
                        type=str,
                        help="")
    parser.add_argument("gyri_dirname",
                        type=str,
                        help="The directory where find all matrices.")
    parser.add_argument("nomenclature",
                        type=str,
                        help="The nomenclature file where each parcel name is "
                             "associated at a label.")
    parser.add_argument("freesurfer_nomenclature",
                        type=str,
                        help="The FreeSurfer nomenclature file where each "
                             "gyrus name is associated at a label.")
    parser.add_argument("reduced_matrix",
                        type=str,
                        help="The output reduced matrix wtih the size "
                             "(X regions, regions).")
    parser.add_argument("-f", "--freesurfer",
                        action="store_true",
                        dest="freesurfer_regions",
                        help="Use the FreeSufer nomenclature.")
    parser.add_argument("-a", "--aal",
                        action="store_true",
                        dest="aal_regions",
                        help="Use the AAL nomenclature.")
    parser.add_argument("-c", "--constellation",
                        action="store_true",
                        dest="constellation_regions",
                        help="Use the Constellation nomenclature.")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    """Execute the command to calculate a connectome and its snapshot.
    """
    # Load the arguments of parser (ignore script name: sys.arg[0])
    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    #--------------------------------------------------------------------------
    # STEP 1: Load all labels in the cortical parcellation
    parcellation = aims.read(args.cortical_parcellation)
    parcels = parcellation[0].arraydata()
    labels = list(numpy.unique(parcels))
    
    freesurfer_parcellation = aims.read(args.freesurfer_cortical_parcellation)
    freesurfer_parcels = freesurfer_parcellation[0].arraydata()

    # Ignore the unknow regions (right and left)
    # FreeSurfer
    if args.freesurfer_regions:
        labels.remove(-1)
        labels.remove(36)
    # AAL
    if args.aal_regions:
        labels.remove(0)
    # Constellation
    if args.constellation_regions:
        labels.remove(11)
        labels.remove(407)

    #--------------------------------------------------------------------------
    # STEP 2: Load the nomenclature (labels and names)
    names = numpy.loadtxt(args.nomenclature, dtype=str)

    # List the names of the FreeSurfer or AAL or Constellation nomenclature
    list_names = []
    for label in labels:
        for name in names:
            if int(name[0]) == label: 
                list_names.append(name[1])
    reorganize_names = list_names
    reorganize_labels = labels

    # Only with the AAL nomenclature
    if args.aal_regions:
        # Reorganize the names left/rigth
        reorganize_names = []
        for region in list_names:
            if region[-1:] == "L":
                reorganize_names.append(region)
        for region in list_names:
            if region[-1:] == "R":
                reorganize_names.append(region)
        # Reorganize the labels from names
        reorganize_labels = []
        for rname in reorganize_names:
            for name in names:
                if rname == name[1]:
                    reorganize_labels.append(int(name[0]))

    #--------------------------------------------------------------------------
    # STEP 3: Compute the reduced individual connectome (X regions, X regions)
    if args.freesurfer_regions:
        matrices = []
        for region_name in list_names:
            if region_name == "rh.unknown" or region_name == "lh.unknown":
                pass
            else:
                matrices.append(
                    glob.glob(args.gyri_dirname + "/" + region_name +
                              "/matrix/*_complete_matrix_smooth3.0*.imas")[0])
    else:
        matrices = glob.glob(args.gyri_dirname +
                             "/*/matrix/*_complete_matrix_smooth3.0*.imas")
        if not matrices:
            matrices = glob.glob(os.path.join(args.gyri_dirname, '*.imas'))

    # initialyze a matrix M(vertices, labels)
    clusters_matrix = numpy.zeros((len(parcels), len(labels)))

    all_rows = []
    for num, matrix in enumerate(matrices):
        print('add matrix', num, ' /', len(matrices), ':', matrix)
        mat = aims.read(matrix)
        m = numpy.asarray(mat.denseMatrix())[:, :, 0, 0]
        if m.shape[0] < m.shape[1]:
            m = numpy.transpose(m)
        if args.freesurfer_regions:
            columns = numpy.sum(m, axis=0)
        label_row = []
        # Add the input matrix lines corresponding to the label/region
        # M(cortical vertices, input region vertices)
        for label in reorganize_labels:
            idx = numpy.where(parcels == label)[0]
            if args.freesurfer_regions:
                label_row.append(numpy.sum(columns[idx]))
            else:
                sumcol = numpy.sum(m[idx], axis=0)
                label_row.append(sumcol)
        if args.freesurfer_regions:
            all_rows.append(numpy.asarray(label_row))
        # AAL or Ccnstellation parcellation
        # Sort all the vertices such as cortical surface vertices are sorted
        else:
            cortical_region = os.path.basename(
                os.path.dirname(os.path.dirname(matrix)))
            label_number = select_ROI_number(
                args.freesurfer_nomenclature, cortical_region)
            idx_vertices = numpy.where(
                freesurfer_parcels == int(label_number))[0]
            for l, lr in enumerate(numpy.asarray(label_row).T):
                if l < len(idx_vertices):
                    clusters_matrix[idx_vertices[l]] = lr
    # AAL or Ccnstellation parcellation
    # Add the lines of intermediate matrix corresponding to the label/region
    # M(corticcal vertices, Regions)
    if not args.freesurfer_regions:
        new_label_row = []
        for label in reorganize_labels:
            idx = numpy.where(parcels == label)[0]
            new_label_row.append(numpy.sum(clusters_matrix[idx], axis=0))
        all_rows = new_label_row
    
    # Create the reduced matrix (Regions, Regions)
    connectome = numpy.asarray(all_rows)

    # Write the reduced matrix on the disk.
    numpy.save(args.reduced_matrix, connectome)

    #--------------------------------------------------------------------------
    # STEP 4: Create the figure with matplotlib
    fig, ax = plt.subplots()
    heatmap = ax.imshow(connectome,
                        cmap=plt.cm.nipy_spectral,
                        interpolation='none')

    ax.xaxis.tick_top()

    ax.set_xticks(numpy.arange(0.1, connectome.shape[0]))
    ax.set_yticks(numpy.arange(0.1, connectome.shape[1]))

    ax.tick_params(which="both", axis="both", width=0, length=0)

    if args.freesurfer_regions:
        labels_size = 9
    elif args.aal_regions:
        labels_size = 5
    else:
        labels_size = 2

    ax.set_xticklabels(reorganize_names, size=labels_size, rotation=90)
    ax.set_yticklabels(reorganize_names, size=labels_size)

    ax.set_aspect("equal")

    colorbar_title = "number of connections"
    colorbar = fig.colorbar(heatmap, ticks=[0, 1000000])
    colorbar.ax.set_yticklabels(['0', '1000000'])
    colorbar.set_label(colorbar_title, rotation=270, labelpad=20)

    fig.tight_layout()
    snapshot = os.path.splitext(args.reduced_matrix)[0]

    # Save to PNG file
    fig.savefig(snapshot, dpi=200)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
