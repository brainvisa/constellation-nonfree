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
Misc utilities.
"""

# ---------------------------Imports-------------------------------------------

# System module
from __future__ import print_function
from __future__ import absolute_import
import logging
from six.moves import range

# Define logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create handler (console)
steam_handler = logging.StreamHandler()
steam_handler.setLevel(logging.INFO)
logger.addHandler(steam_handler)


# ---------------------------Functions-----------------------------------------


def sameNbElements(listA, listB, NonZeroRealNb=True):
    """Return lists with the same number of elements.

    Parameters
    ----------
        listA
        listB
        NonZeroRealNb

    Returns
    -------
        listA
        listB
    """
    for index in range(len(listA)):
        if listA[index] == 0 and listB[index] != 0:
            listB[index] = 0
        elif listA[index] != 0 and listB[index] == 0:
            listA[index] = 0
    if NonZeroRealNb:
        listA = [x for x in listA if x != 0]
        listB = [x for x in listB if x != 0]
    if len(listA) != len(listB):
        print('ERROR: The A size is different from the B.')
    return listA, listB


def check_no_empty_list(listof):
    """Check no empty list.

    Parameters
    ----------
        listof: (mandatory)
    """
    if not listof:
        logger.error("List is empty.")
    else:
        logger.info("List has elements.")
