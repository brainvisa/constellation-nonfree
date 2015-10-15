###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

"""
This script does the following:
*
*

Main dependencies:

Author: sandrine.lefranc@cea.fr
"""

#----------------------------Imports-------------------------------------------


#----------------------------Functions-----------------------------------------


def read_file(filename):
    """Read the first column of a txt file.

    Parameter
    ---------
    filename : str (mandatory)
        the file gives the nomenclature of the ROIs (name and number)
        extension of this file: .txt

    Return
    ------
    nomenclature : numpy.array
        the first column of the file = names of the ROIs
    """
    with open(filename, "r") as inf:
        ls = inf.readlines()
        nomenclature = []
        for l in ls:
            nomenclature.append(l.split()[0])  # first column
    return nomenclature
