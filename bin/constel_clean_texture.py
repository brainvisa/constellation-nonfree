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

from constel.lib.utils.texturetools import clean_array

# ---------------------------Command-line--------------------------------------

doc = """Clean labels of a constellation texture"""


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
        'cleaned_texture',
        type=str,
        help='Path to the destination cleaned texture')

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

    # clean the numpy array
    cleaned_array = clean_array(array)

    # create and assign the array to a texture
    cleaned_texture = aims.TimeTexture_S16()
    cleaned_texture[0].assign(cleaned_array)

    # write the texture on disk
    aims.write(cleaned_texture, args.cleaned_texture)


if __name__ == "__main__":
    main()
