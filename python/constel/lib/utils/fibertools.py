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
Fiber utilities.
"""

# ---------------------------Imports-------------------------------------------


# python module
import logging
import os

# soma module
from soma import aims

# define logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create handler (console)
steam_handler = logging.StreamHandler()
steam_handler.setLevel(logging.DEBUG)
logger.addHandler(steam_handler)


# ---------------------------Functions-----------------------------------------


def load_fiber_tracts(directory, formats):
    """Load all fiber tracts from a given directory.

    Parameters
    ----------
    directory: str (mandatory)
        pattern of fiber tracts files to download
    formats: str (mandatory)
        fiber tracts formats (.bundles or .trk)

    Return
    ------
    fnames: list
        list of bundles files
    """
    bundles = []
    if formats == "bundles":
        for root, dirs, files in os.walk(directory):
            for f in files:
                end = os.path.basename(os.path.splitext(f)[0])[-1]
                if f.endswith(".bundles") and end.isdigit():
                    bundles.append(os.path.join(root, f))
    elif formats == "trk":
        for root, dirs, files in os.walk(directory):
            for f in files:
                end = os.path.basename(os.path.splitext(f)[0])[-1]
                if f.endswith(".trk") and end.isdigit():
                    bundles.append(os.path.join(root, f))
    else:
        logger.error("Format unsupported. Use .bundles or .trk formats.")

    return bundles


def write_transfo2file(transformation_matrix, filename):
    """Write the transfo (dw tot t1) in the fibertracts file.

    Parameters
    ----------
    transformation_matrix: str (mandatory)
        filename of the transformation matrix
    filename: str (mandatory)
        filename of the fiber tracts

    Return
    ------
    fibertracts: soma.aims.Graph
    """
    fibertracts = aims.read(filename)
    transfo = aims.read(transformation_matrix)
    fibertracts["AffineTransformation3d"] = [transfo.toMatrix()]
    return fibertracts


def is_empty(bundle_filename):
    """Check if the bundle file is empty.

    Parameter
    ---------
    bundle_filename: (string)
        name of bundles with directory

    Return
    ------
    empty: (boolean)
        True if empty
        False otherwise
    """
    f = aims.Finder()
    f.check(bundle_filename)
    curves_count = f.header()['curves_count']
    empty = (curves_count <= 0)
    return empty
