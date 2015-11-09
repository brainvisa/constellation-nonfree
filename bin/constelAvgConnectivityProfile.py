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
*
*

Main dependencies: PyAims library

Author: Sandrine Lefranc, 2015
"""


#----------------------------Imports-------------------------------------------


# python system module
import sys
import json
import numpy
import argparse
import textwrap
import subprocess
import exceptions

# aims module
from soma import aims


#----------------------------Functions-----------------------------------------


def mylist(string):
    return json.loads(string)


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        "normprofiles", type=mylist, help="")
    parser.add_argument("oprofiles", type=mylist, help="")
    parser.add_argument("meanprofile", type=str, help="")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    # Load the arguments of parser (delete script name: sys.arg[0])
    arguments = (json.dumps(eval(sys.argv[1])),
                 json.dumps(eval(sys.argv[2])),
                 sys.argv[3])
    parser, args = parse_args(arguments)

    #for normprofile, oprofile in zip(args.normprofiles, args.oprofiles):
        #aimsprofile = aims.read(normprofile)
        #profile = aimsprofile[0].arraydata()
        #if profile.sum() > 0:
            #z = 1./ profile.sum()
            #newprofile = profile*z
        #aimstex = aims.TimeTexture_FLOAT()
        #aimstex[0].assign(newprofile)
        #aims.write(aimstex, oprofile)
        
    # sum all the profiles
    for idx, tex_filename in enumerate(args.normprofiles):
        tex = aims.read(tex_filename)[0].arraydata()
        if idx == 0:
            t = tex
            #shutil.copy2(str(tex_filename), str(args.meanprofile))
        else:
            t = numpy.add(t,tex)
            #cmd = ["AimsLinearComb",
                   #"-i", tex_filename,
                   #"-j", args.meanprofile,
                   #"-o", args.meanprofile]
            #subprocess.check_call(cmd)
            
    # mean profile
    meanprofile = t / len(args.normprofiles)

    # create a time texture object by assigning the mean profile
    aimstex = aims.TimeTexture_FLOAT()
    aimstex[0].assign(meanprofile)

    # write the file (mask) on the disk
    aims.write(aimstex, args.meanprofile)

    #meanprofile = aims.read(args.meanprofile)[0].arraydata()
    #new_profile = meanprofile / len(args.normprofiles)
       
    #aimstex = aims.TimeTexture_FLOAT()
    #aimstex[0].assign(new_profile)
    #aims.write(aimstex, args.meanprofile)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()
