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
import os
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
            Constellation Data Check
            ~~~~~~~~~~~~~~~~~~~~~~~~
            Run this command to check if your input Constellation processing
            home directory contains the expected Constellation output files for
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
        help="output file to write the result of datacheck")
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
        "diffusion/default_acquisition/default_analysis/"
        "default_tracking_session/connectivity_parcellation/",
        args.method, args.studyname, args.roi)
    # define the expected organization of the Constellation output files
    constel_struct = {
        "matrix": 5,
        "filteredTracts": 6,
        "clustering": 1,
        os.path.join("clustering", "smooth3.0"): 4,
        "extrapaths": [""]}

    # detect each subject in the Constellation processing homedir "consteldir"
    status = {}
    consteldir = args.ctdir + "/subjects"
    if args.subjectid:
        ctdir_subfolders = [args.subjectid]
    else:
        ctdir = os.listdir(consteldir)
        ctdir_subfolders = [d for d in ctdir if not d.endswith(".minf")]
        if "group_analysis" in ctdir_subfolders:
            ctdir_subfolders.remove("group_analysis")

    for sid in ctdir_subfolders:
        # Store the subject tree folder-folder files counts
        status[sid] = {"extrapaths": []}
        siddir = os.path.join(consteldir, sid, subdir)
        for path, dirs, files in os.walk(siddir):
            rpath = path.replace(siddir, "").lstrip(os.sep)
            if rpath in constel_struct:
                filedir = os.listdir(os.path.join(siddir, rpath))
                files = [d for d in filedir if not d.endswith(".minf")]
                status[sid][rpath] = len(files)
            else:
                status[sid]["extrapaths"].append(rpath)

    # depending on the verbosity, display the number of complete/failed
    # processings, or the list of subjects that didn't run properly.
    failed_sids = []
    total_success = 0
    total_failed = 0
    total = len(status)
    for sid, sid_status in status.items():
        if sid_status == constel_struct:
            total_success += 1
        else:
            total_failed += 1
            failed_sids.append(sid)
    with open(args.ofile, 'w') as f:
        f.write("SUCCESS: {0}/{1}".format(total_success, total))
        f.write("\nFAILED: {0}/{1}".format(total_failed, total))
        f.write("\nFAILED SIDS: {0}".format(failed_sids))
        if args.subjectid:
            f.write("\nTREE '{0}': ".format(args.subjectid))
            for key, value in status[args.subjectid].items():
                f.write("{0}: {1} (observed) - {2} (reference)".format(
                    key, value, constel_struct[key]))


if __name__ == "__main__":
    main()
