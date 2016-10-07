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
This script does the following:
* read the file to create a list with the content of the file
* select the number corresponding to a name of ROI from a given file

Main dependencies:

Author: Sandrine Lefranc, 2015
"""

#----------------------------Imports-------------------------------------------


#----------------------------Functions-----------------------------------------


def read_file(filename, mode=0):
    """Read the file to create a list with the content of the file.

    Parameter
    ---------
    filename : str (mandatory)
        the file gives the nomenclature of the ROIs (name and number)
        extension of this file: .txt
    mode: int (mandatory)
        there are two modes:
            - 0 to read the entire file
            - <int> to read the column corresponding to <int>th column of file
              example: mode=1, the first column is read.
        by default, mode=0

    Return
    ------
    nomenclature : list
        list of the content of the file (entire file or one column of the file)
    """
    with open(filename, "r") as inf:
        ls = inf.readlines()
        if mode == 0:
            nomenclature = {}
            for l in ls:
                nomenclature[l.split()[0]] = l.split()[1]
        else:
            nomenclature = []
            for l in ls:
                nomenclature.append(l.split()[mode - 1])
    return nomenclature


def select_ROI_number(filename, ROIname):
    """Select the number corresponding to a name of ROI from a given file.

    Parameters
    ----------
    filename: str (mandatory)
        the file with all names and labels of ROIs
    ROIname: str (mandatory)
        the selected ROI name in the file

    Return
    ------
    ROIlabel: int
    """
    s = read_file(filename)

    ROIlabel = s.keys()[s.values().index(ROIname)]
    #ROIlabel = list(s.keys())[list(s.values()).index(ROIname) # python3
    return ROIlabel


def delete_regions_in_nomenclature(nomenclature, regions, up_nomenclature):
    """ Remove some regions in the cortical region nomenclature.
    
    Parameters
    ----------
    nomenclature: string (mandatory)
        The filename with the correspondance between the name region and its
        number.
    regions: list (mandatory)
        List of the cortical regions to delete.
    up_nomenclature: string (mandatory)
        The name (with directory) of the nomenclature without the regions.
    
    Returns
    -------
    up_nomenclature: list
        The nomenclature without the regions.
    """
    with open(nomenclature) as f:
        #content = f.readlines()
        with open(up_nomenclature, "w") as output:
            for line in f:
                if line.split()[0] not in regions:
                    output.write(line)

    return up_nomenclature


def add_region_in_nomenclature(nomenclature, new_region, number_region):
    """Add a new name in the cortical region nomenclature.
    
    Parameters
    ----------
    nomenclature: string (mandatory)
        The filename with the correspondance between the name region and its
        number.
    new_region: string (mandatory)
        The new name for inclusion into the nomenclature.
    number_region: int (mandatory)
        The number associated to the new region.
    """
    with open(nomenclature, "a") as f:
        new_line = str(number_region) + " " + new_region + "\n"
        f.writelines(new_line)
    

    
