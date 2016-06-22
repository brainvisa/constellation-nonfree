#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

# system modules
import optparse
import sys

# constel module
from constel.lib.utils.fibertools import load_fiber_tracts
from constel.lib.misctools import check_no_empty_list


def parseOpts(argv):
    desc = """usage: %prog [options] filename"
    Load fiber tracts from a directory.
    """

    parser = optparse.OptionParser(desc)

    parser.add_option("-d", "--directory", dest="directory",
                      help="directory containing the fiber tracts")
    parser.add_option("-f", "--formats", dest="formats",
                      help="format of fiber tracts: .bundles or .trk")
    parser.add_option("-o", "--fibertracts", dest="fibertracts", type=str,
                      help="List of fiber tracts for one subject.")

    return parser, parser.parse_args(argv)

def main():
    parser, (options, args) = parseOpts(sys.argv)
    
    fibertracts = load_fiber_tracts(options.directory, options.formats)
    options.fibertracts = str(fibertracts)
    print type(options.fibertracts)
    check_no_empty_list(options.fibertracts)

    return options.fibertracts

if __name__ == "__main__":
    main()
