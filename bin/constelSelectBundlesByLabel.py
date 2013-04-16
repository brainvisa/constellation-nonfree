#!/usr/bin/env python
import sys
from optparse import OptionParser
import constel.lib.bundles.bundles_select_by_label as SB


def parseOpts(argv):
  description = 'Select bundles with at least label --label in the bundlename.'
  parser = OptionParser(description)
  parser.add_option('-i', '--inputparcelstex', dest='parcels_texname',
    metavar = 'FILE', help='input texture name of a labeled texture representing the parcel regions, between 0 and parcelRegionsNb, 0 = background, 1 to parcelRegionsNb = parcel regions', default = None)
  parser.add_option('-l', '--label', dest='parcel_label',
    metavar = 'FILE',help= 'parcel label, between 1 to parcelRegionsNb')
  parser.add_option('-b', '--bundles', dest='input_patch_shake_up_bundles',
    metavar = 'FILE',help= 'input bundles filename, with bundles regrouped accroding to the inpur parcellation texture')
  parser.add_option('-o', '--obundles', dest ='out_label_bundles',
    metavar = 'FILE', action='store', help ='select all the fibers which intersect the parcel of label --label and regroup them into this bundle')
  parser.add_option('-s', '--addbundlesname', dest='addBundlesName',
    metavar = 'FILE',help= 'BundlesNames to add to the selected bundles', default = "")
  parser.add_option('-v', '--verbose', dest ='verbose', default = None)
  
  
  return parser, parser.parse_args(argv)


def main():
  parser, (options, args) = parseOpts(sys.argv)
  patch_label = int(options.parcel_label)
  SB.selectPatchBundles(int(options.parcel_label), options.parcels_texname, options.input_patch_shake_up_bundles, options.out_label_bundles, options.verbose, addlabel_list = [options.addBundlesName])
  
if __name__ == "__main__" : main()





