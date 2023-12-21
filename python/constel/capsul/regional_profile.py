###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

import subprocess
from capsul.api import Process
from soma.controller import File, field
from constel.lib.utils.filetools import select_ROI_number


class RegionalProfile(Process):
    """
    Run the command 'constelMeanConnectivityProfileFromMatrix'.
    Computes the connectivity profile of a region from a connectivity matrix.
    """
    atlas: field(type_=File,
                 doc="Surface atlas used as base parcellation. Number of vertices of the choosen surface referential.")
    nomenclature: field(type_=File,
                        doc="Regions nomenclature of the surface atlas.")
    region: field(type_=str,
                  doc="Region of interest.")
    matrix: field(type_=File,
                  doc="Connectivity matrix. Shape: (n, n) with n being the size of the choosen surface referential.")
    erase_matrix: field(type_=bool,
                        doc="Erase the input matrix if its too large for storage capacity. Especially useful for smoothed matrix.")
    profile: field(type_=File,
                   output=True,
                   write=True,
                   doc="Regional connectivity profile obtained by averaging connectivity profiles of vertices of the selected region. Size of the choosen surface referential.")

    def execute(self, context=None):
        # selects the label corresponding to region name
        label = select_ROI_number(self.nomenclature,
                                  self.region)

        cmd = ["constelMeanConnectivityProfileFromMatrix",
               "-connfmt", "binar_sparse",
               "-connmatrixfile", self.matrix,
               "-outconntex", self.profile,
               "-seedregionstex", self.atlas,
               "-seedlabel", str(label),
               "-type", "seed_mean_connectivity_profile",
               "-normalize", "0",
               "-verbose", "1"]
        subprocess.check_call(cmd)
        if self.erase_matrix:
            self.matrix.eraseFiles(remove_from_database=True)
