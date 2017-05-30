###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################


# ---------------------------Imports-------------------------------------------

# System module
import os
import numpy
import argparse

# Soma library
from soma import aims

# ---------------------------Command-line--------------------------------------

description = """Compute the connectogram from the FSL outputs."""

parser = argparse.ArgumentParser(description=description)

parser.add_argument(
    'vertex_labels',
    type=str,
    help='List of labels for both hemispheres.')
parser.add_argument(
    'fdt_matrix',
    type=str,
    help='Vertex to vertex connectivity.')
parser.add_argument(
    'coords',
    type=str,
    help='Coords to vertex connectivity.')
parser.add_argument(
    'label',
    type=int,
    help='Label.')
parser.add_argument(
    'outdir',
    type=str,
    help='Ouput directory.')

args = parser.parse_args()

# ---------------------------Main program--------------------------------------

# define the variables
vlabels = args.vertex_labels
fdt_matrix = args.fdt_matrix
outdir = args.outdir
label = args.label
coords = args.coords

if not os.path.isdir(outdir):
    os.makedirs(outdir)

# load the labels texture
aims_vlabels = aims.read(vlabels)
np_vlabels = aims_vlabels[0].arraydata()

#
idx = numpy.where(np_vlabels == label)[0]
inv_idx = numpy.zeros((len(np_vlabels)), dtype=int) - 1
inv_idx[idx] = range(len(idx))

# Initialize the matrix dimension
connectome = numpy.zeros((len(idx), len(np_vlabels)))

coords_list = []
nv_hemi = len(np_vlabels) / 2
with open(coords, 'r') as f:
    for line in f:
        elements = line.strip().split()
        hemi, vertex = int(elements[-2]), int(elements[-1])
        coords_list.append(nv_hemi * hemi + vertex)

with open(fdt_matrix, 'r') as f:
    for line in f:
        elements = line.strip().split()
        v1_idx, v2_idx, nb_connections \
            = (coords_list[int(elements[0]) - 1],
               coords_list[int(elements[1]) - 1],
               float(elements[2]))
        if inv_idx[v1_idx] >= 0:
            connectome[inv_idx[v1_idx], v2_idx] += nb_connections
        if inv_idx[v2_idx] >= 0:
            connectome[inv_idx[v2_idx], v1_idx] += nb_connections

# Write the matrix as Volume on the disk
vol = aims.Volume(connectome.astype(float))
smat = aims.SparseOrDenseMatrix()
smat.setMatrix(vol)
output_connectome = os.path.join(outdir,
                                 'connectome_label' + str(label) + '.imas')
aims.write(smat, output_connectome)
