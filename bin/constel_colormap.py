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
import json

# brainvisa import
from soma import aims

# constel import
from constel.lib.utils.colormaptools import get_colormap_colors


# ---------------------------Command-line--------------------------------------

doc = """Color a parcellation with a few colors"""


def mylist(string):
    return json.loads(string)


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
        'colors_dict',
        type=str,
        help='Optional dictionary of colors'
    )
    parser.add_argument(
        'nb_colors',
        type=str,
        help='Number of colors to use during the dsatur algorithm'
    )
    parser.add_argument(
        'default_labels',
        type=mylist,
        help='Labels for which the palette will be set to default color'
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

    #  default list of colors : Default, Red, Green, Blue, Yellow, Purple,
    #  Cyan, Orange, Grey
    default_colors = {0: [255, 255, 255, 255],
                      1: [203, 0, 0, 255],
                      2: [64, 173, 38, 255],
                      3: [0, 0, 142, 255],
                      4: [255, 217, 0, 255],
                      5: [186, 85, 211, 255],
                      6: [0, 255, 255, 255],
                      7: [255, 165, 0, 255],
                      8: [100, 100, 100, 255],
                      9: [102, 0, 51, 255],
                      10: [204, 255, 153, 255],
                      11: [0, 153, 153, 255],
                      12: [255, 153, 153, 255],
                      13: [255, 255, 204, 255]}

    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    # load the files (mesh and gyri segmentation)
    mesh = aims.read(args.mesh)
    texture = aims.read(args.texture)

    default_labels = []
    for label in args.default_labels:
        default_labels.append(int(label))

    try:
        aims_colors_dict = aims.read(args.colors_dict)
        colors_dict = {}
        for key in aims_colors_dict.keys():
            colors_dict[int(key)] = [int(x) for x in aims_colors_dict[key]]
    except:
        colors_dict = default_colors

    mode = len(colors_dict[list(colors_dict.keys())[0]])

    if 0 not in colors_dict.keys():
        if mode == 3:
            colors_dict[0] = [255, 255, 255]
        elif mode == 4:
            colors_dict[0] = [255, 255, 255, 255]

    # Support RGB
    RGBA_colors = get_colormap_colors(mesh,
                                      texture,
                                      colors_dict,
                                      args.nb_colors,
                                      default_labels)
    # write as volume
    arr = np.asarray(RGBA_colors, dtype='uint8').reshape(
        int(len(RGBA_colors)/4), 1, 1, 1, mode)
    if mode == 3:
        vol = aims.Volume('RGB', arr.shape[0])
    elif mode == 4:
        vol = aims.Volume('RGBA', arr.shape[0])
    vol['v'] = arr
    aims.write(vol, args.palette)


if __name__ == "__main__":
    main()
