###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################


# ----------------------------Imports-------------------------------------------

import argparse
import textwrap
import sys

from soma import aims
import numpy as np
import json

from constel.lib.utils.texturetools import cut_regions_from_array

# ---------------------------Command-line--------------------------------------

doc = """Cut a region from a constellation texture"""


def mylist(string):
    return json.loads(string)


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(doc))

    parser.add_argument(
        'original_texture',
        type=str,
        help='Texture to clean')
    parser.add_argument(
        "labels",
        type=mylist,
        help='Labels to cut in texture')
    parser.add_argument(
        'cut_texture',
        type=str,
        help='Path to the destination cut texture')
    parser.add_argument(
        'cut_regions',
        type=str,
        help='Path to the destination cut region')

    args = parser.parse_args(argv)

    return parser, args

# ----------------------------Main program-------------------------------------


def main():
    """
    """

    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    # import texture as numpy array
    texture = aims.read(args.original_texture)
    array = np.asarray(texture)[0]
    labels = []
    for label in args.labels:
        labels.append(int(label))

    cut_regions_array, cut_array = cut_regions_from_array(args.labels, array)

    # create and assign arrays to textures
    cut_texture = aims.TimeTexture_S16()
    cut_regions = aims.TimeTexture_S16()
    cut_texture[0].assign(cut_array)
    cut_regions[0].assign(cut_regions_array)

    # write the texture on disk
    aims.write(cut_texture, args.cut_texture)
    aims.write(cut_regions, args.cut_regions)


if __name__ == "__main__":
    main()
