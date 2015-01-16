#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

# System module
import optparse
import sys

# constel modules
from constel.lib.texturetools import management_internal_connections
from constel.lib.texturetools import normalize_profile

# soma module
from soma import aims


def parseOpts(argv):
    desc = """usage: %prog [options] filename"
    Remove the patch internal connections of the cortical connectivity profile.
    Then, the profile is normalized by its total number of connections.
    """
    
    parser = optparse.OptionParser(desc)
     
    parser.add_option("-p", "--profile", dest="profile",
                      help="a coritcal connectivity profile")
    parser.add_option("-g", "--gyriseg", dest="gyriseg",
                      help="a labeling of gyri cortical segmentation")
    parser.add_option("-c", action="store_true", dest="internal_connections")
    parser.add_option("-q", action="store_false", dest="internal_connections")
    parser.add_option("-n", "--nprofile", dest="normprofile",
                      help="normalized connectivity profile")

    return parser, parser.parse_args(argv)


if __name__ == "__main__":
    parser, (options, args) = parseOpts(sys.argv)
    
    new_profile = management_internal_connections(
        options.gyriseg, options.profile, options.internal_connections)
    
    normalized_profile = normalize_profile(new_profile)
    
    aims.write(normalized_profile, options.normprofile)
