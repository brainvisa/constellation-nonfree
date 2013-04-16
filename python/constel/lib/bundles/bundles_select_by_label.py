from soma import aims
import constel.lib.misctools as OT
import os

def selectPatchBundles(patch_label, patchsMeshTexture_filename, subjectPatchsShakeUpBundles_filename, outputSubjectPatchBundles_filename, verbose = False, addlabel_list = [""] ):
  """
  reduces connection profiles for a given patch on a subject mesh:
  inputs:
  patch_label
  patchsTexture_filename
  subject_mesh
  subjectPatchsShakeUpBundles_filename
  outputSubjectPatchBundles_filename
  
  outputs:
  outputSubjectPatchBundles_filename
  """
  patchs_tex = aims.read(patchsMeshTexture_filename)
  patchs_nb = int(patchs_tex[0].arraydata().max())
  names_list, names_dict = OT.generateIntPairsNames(patchs_nb)
  names_patch = names_dict[patch_label]
  names = " ".join( names_patch )

  names += str(-1) + "_" + str(patch_label) + " "
  if verbose:
    print names
  names += ' '.join( [ str(x) for x in addlabel_list ] )
  command_BundlesPatchSelection = "selectBundlesFromNames -i " + subjectPatchsShakeUpBundles_filename + " -o " + outputSubjectPatchBundles_filename + " -names " + names  + " \n\n"

  if verbose:
    print command_BundlesPatchSelection

  os.system(command_BundlesPatchSelection)
  return True

