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

from constel.lib.utils.texturetools import cut_region_from_array

# ---------------------------Command-line--------------------------------------

doc = """Cut a region from a constellation texture"""


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
        'label',
        type=int,
        help='Number of the label to cut in texture')
    parser.add_argument(
        'cut_texture',
        type=str,
        help='Path to the destination cut texture')
    parser.add_argument(
        'cut_region',
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
    label = args.label

    cut_region_array, cut_array = cut_region_from_array(label, array)

    # create and assign arrays to textures
    cut_texture = aims.TimeTexture_S16()
    cut_region = aims.TimeTexture_S16()
    cut_texture[0].assign(cut_array)
    cut_region[0].assign(cut_region_array)

    # write the texture on disk
    aims.write(cut_texture, args.cut_texture)
    aims.write(cut_region, args.cut_region)


if __name__ == "__main__":
    main()
