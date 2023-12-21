###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

from capsul.api import Process
from soma.controller import File, Directory, Literal, field


class SmoothMatrix(Process):
    """Run the command 'AimsSparseMatrixSmoothing'.
       Smoothing of the individual connectivity matrix."""
    atlas: field(type_=File,
                 doc="Surface atlas used as base parcellation. Number of vertices of the choosen surface referential.")
    nomenclature: field(type_=File,
                        doc="Regions nomenclature of the surface atlas.")
    region: field(type_=str,
                  doc="Region of interest.")
    individual_matrix: field(type_=File,
                             doc="Complete individual connectivity matrix. This matrix should be sparse. Shape: (n, n) with n being the size of the choosen surface referential.")
    individual_white_mesh: field(type_=File,
                                 doc="White-grey individual interface. Number of vertices of the choosen surface referential.")
    smoothing_value: field(type_=float,
                           doc="Smoothing value to apply in millimeters. Default to 3.0 mm.")
    smoothed_matrix: field(type_=File,
                           output=True,
                           write=True,
                           doc="Smoothed individual connectivity matrix. Shape: (n, n) with n being the size of the choosen surface referential.")

    def execute(self, context=None):
        import subprocess
        from constel.lib.utils.filetools import select_ROI_number
        from constel.lib.utils.matrixtools import replace_negative_values

        # selects the label number corresponding to label name
        label_number = select_ROI_number(self.nomenclature,
                                         self.region)

        # matrix smoothing: -s in millimetres
        cmd = ["AimsSparseMatrixSmoothing",
               "-i", self.individual_matrix,
               "-m", self.individual_white_mesh,
               "-o", str(self.smoothed_matrix),
               "-s", str(self.smoothing_value),
               "-l", self.atlas,
               "-p", str(label_number)]
        print('call')
        subprocess.check_call(cmd)
        print('done')
        #replace_negative_values(self.individual_matrix)
