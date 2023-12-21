from capsul.api import Capsul

atlas = '/casa/host/build/share/constellation-5.2/constellation_atlas_hcp_200s/freesurfer_gyri/group_analysis/1113S/average_brain/bh.annot.averagebrain.gii'
nomenclature = '/volatile/clanglet/constel_utils/nomenclature_desikan_freesurfer.txt'
matrix = '/volatile/clanglet/100206/100206_hcp_lh.bankssts_complete_matrix_smooth0.0_20.0to500.0mm.imas'
mesh = '/volatile/clanglet/100206/bh.r.aims.white.gii'
smatrix = '/volatile/clanglet/100206/100206_hcp_lh.bankssts_complete_matrix_smooth3.0_20.0to500.0mm.imas'


capsul = Capsul()
executable = capsul.executable('constel.capsul.smooth_matrix.SmoothMatrix')
with capsul.engine() as engine:
    engine.run(executable,
               atlas=atlas,
               nomenclature=nomenclature,
               region='lh.bankssts',
               individual_matrix=matrix,
               individual_white_mesh=mesh,
               smoothing_value=3.0,
               smoothed_matrix=smatrix)
