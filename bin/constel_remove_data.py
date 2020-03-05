#! /usr/bin/env python
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

Author: Sandrine Lefranc, 2016
"""

#----------------------------Imports-------------------------------------------


# system module
from __future__ import print_function
from __future__ import absolute_import
import os
import re
import sys
import argparse
import textwrap


#----------------------------Functions-----------------------------------------


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Constellation Remove Data
            ~~~~~~~~~~~~~~~~~~~~~~~~~
            Run this command to remove the Constellation output files for
            all the subjects.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument(
        "ctdir",
        help="the Constellation processing home directory.")
    parser.add_argument(
        "method", type=str,
        help="the method to inspect.")
    parser.add_argument(
        "studyname", type=str,
        help="the study name to inspect.")
    parser.add_argument(
        "roi", type=str,
        help="the ROI name to inspect.")
    parser.add_argument(
        "-s", "--subjectid", dest="subjectid",
        help="the subject identifier.")
    parser.add_argument(
        "ofile",
        help="output file to write the result of removed data")
    parser.add_argument(
        "--ifibers", dest="ifibers", action="store_true",
        help="if activated, remove the fibers near cortex for all subjects")
    parser.add_argument(
        "--ofibers", dest="ofibers", action="store_true",
        help="if activated, remove the outside fibers of cortex for all"
        " subjects")
    parser.add_argument(
        "--ovfibers", dest="ovfibers", action="store_true",
        help="if activated, remove the oversampled fibers for all subjects")
    parser.add_argument(
        "--imatrix", dest="imatrix", action="store_true",
        help="if activated, remove the fibers near cortex matrix for all"
        " subjects")
    parser.add_argument(
        "--omatrix", dest="omatrix", action="store_true",
        help="if activated, remove the outside fibers matrix for all subjects")
    parser.add_argument(
        "--smatrix", dest="smatrix", action="store_true",
        help="if activated, remove the smoothed matrix for all subjects")
    parser.add_argument(
        "--ifiberp", dest="ifiberp", action="store_true",
        help="if activated, remove the fibers near cortex profile for all"
        " subjects")
    parser.add_argument(
        "--ofiberp", dest="ofiberp", action="store_true",
        help="if activated, remove the outside fibers profile for all"
        " subjects")
    parser.add_argument(
        "--meanp", dest="meanp", action="store_true",
        help="if activated, remove the mean profile for all subjects")
    parser.add_argument(
        "--normp", dest="normp", action="store_true",
        help="if activated, remove the normed profile for all subjects")
    parser.add_argument(
        "-v", "--verbose", dest="verbose", type=int,
        choices=[0, 1, 2], default=0,
        help="increase the verbosity level: 0 silent, [1, 2] verbose.")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    """
    """
    # load the arguments of parser (delete script name: sys.arg[0])
    arguments = (sys.argv[1:])
    parser, args = parse_args(arguments)

    # define the directory with method, study and roi name.
    subdir = os.path.join(
        args.ctdir,
        "diffusion/default_acquisition/default_analysis/"
        "default_tracking_session/connectivity_parcellation/",
        args.method, args.studyname, args.roi)

    with open(args.ofile, 'w') as f:
        if args.ifibers:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "filteredTracts"))
                if re.match('.*_fibersNearCortex_*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.ofibers:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "filteredTracts"))
                if re.match('.*_outsideFibersOfCortex_*mm.bundles*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.ovfibers:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "filteredTracts"))
                if re.match('.*_oversampled*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.imatrix:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "matrix"))
                if re.match('.*connectivityMatrix_fibers*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.omatrix:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "matrix"))
                if re.match('.*connectivityMatrix_outside*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.smatrix:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "matrix"))
                if re.match('.*connectivityMatrixSmooth*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.ifiberp:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "matrix"))
                if re.match('.*meanConnectivityProfile_fibers*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.ofiberp:
            filenames = [
                i for i in os.listdir(os.path.join(subdir, "matrix"))
                if re.match('.*meanConnectivityProfile_outside*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.meanp:
            filenames = [
                i for i in os.listdir(os.path.join(subdir,
                                      "clustering/smooth3.0"))
                if re.match('.*meanConnectivityProfile*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))
        if args.normp:
            filenames = [
                i for i in os.listdir(os.path.join(subdir,
                                      "clustering/smooth3.0"))
                if re.match('.*normedMeanConnectivityProfile*', i)]
            for filename in filenames:
                if os.path.exists(filename):
                    os.remove(filename)
                    f.write("REMOVE: {0}".format(filename))


if __name__ == "__main__":
    main()
