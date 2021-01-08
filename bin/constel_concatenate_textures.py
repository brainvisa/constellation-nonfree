###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################


#----------------------------Imports-------------------------------------------

import argparse
import textwrap
import sys

from soma import aims
from constel.lib.utils.texturetools import concatenate_texture_new

# ---------------------------Command-line--------------------------------------

doc = """Concatenate textures from a list of paths"""


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(doc))

    parser.add_argument(
        'texture_list',
        type=str,
        help='List of textures to concatenate')
    parser.add_argument(
        'concatenated_texture',
        type=str,
        help='Path to the destination texture')

    args = parser.parse_args(argv)

    return parser, args

# ----------------------------Main program-------------------------------------


def main():
    """
    """

    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    textures = eval(args.texture_list)
    final_rseg = concatenate_texture_new(textures, 10)

    aims.write(final_rseg, args.concatenated_texture)


if __name__ == "__main__":
    main()
