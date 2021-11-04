###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################
import argparse
import textwrap
import sys
import numpy as np

# brainvisa import
from soma import aims

# constel import
from constel.lib.utils.colormaptools import get_colormap_colors


# ---------------------------Command-line--------------------------------------

doc = """Color a parcellation with a few colors"""


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(doc))

    parser.add_argument(
        'texture',
        type=str,
        help='Texture to color')
    parser.add_argument(
        'mesh',
        type=str,
        help='Mesh associated with the texture')
    parser.add_argument(
        'nb_colors',
        type=str,
        help='Number of colors to use during the dsatur algorithm'
    )
    parser.add_argument(
        'palette',
        type=str,
        help='Output palette')

    args = parser.parse_args(argv)
    return parser, args

# ----------------------------Main program-------------------------------------


def main():
    """
    """

    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    # load the files (mesh and gyri segmentation)
    mesh = aims.read(args.mesh)
    texture = aims.read(args.texture)

    RGBA_colors = get_colormap_colors(mesh, texture, args.nb_colors)

    # write as volume
    arr = np.asarray(RGBA_colors, dtype='uint8').reshape(
        int(len(RGBA_colors)/4), 1, 1, 1, 4)
    vol = aims.Volume_RGBA(arr.shape[0])
    vol['v'] = arr
    aims.write(vol, args.palette)


if __name__ == "__main__":
    main()
