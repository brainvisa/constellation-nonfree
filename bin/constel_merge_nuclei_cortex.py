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
Merge the white and nuclei meshes.
"""


#----------------------------Imports-------------------------------------------

from __future__ import print_function

# python system module
import os
import sys
import glob
import argparse
import textwrap
import numpy

# aims module
from soma import aims


#----------------------------Functions-----------------------------------------


def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Merge the white and nuclei meshes.
            Create the output parcellation adding the labels for the
            subcortical structures.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument("white_dir",
                        type=str,
                        help="White mesh directory.")
    parser.add_argument("nuclei_dir",
                        type=str,
                        help="Nuclei mesh directory.")
    parser.add_argument("texture_dir",
                        type=str,
                        help="Cortical parcellation directory.")
    parser.add_argument("out_tex_filename",
                        type=str,
                        help="Output parcellation fielname.")
    parser.add_argument("out_mesh_filename",
                        type=str,
                        help="Output mesh filename.")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    """
    """
    # Load the arguments of parser (ignore script name: sys.arg[0])
    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)


#white_dir = "/neurospin/archi-public/DataBaseArchi/Analysis/4-FreeSurfer/fs_archi_v5.3.0/" "%s/surf"
#nuclei_dir = "/neurospin/archi-private/Users/Analysis/13-Nuclei/inverse_nuclei_mesh/%s"
#texture_dir = "/neurospin/archi-public/DataBaseArchi/Analysis/4-FreeSurfer/fs_archi_v5.3.0/" "%s/label"
#out_mesh_filename = "/neurospin/archi-public/Users/lefranc/white_nuclei_%s.gii"
#out_tex_filename = "/neurospin/archi-public/Users/lefranc/white_nuclei_tex_%s.gii"

    # Defines the inputs and outputs arguments
    white_dir = args.white_dir + "/%s/surf"
    print(white_dir)
    nuclei_dir = args.nuclei_dir + "/%s"
    print(nuclei_dir)
    texture_dir = args.texture_dir + "/%s/label"
    print(texture_dir)
    out_mesh_filename = args.out_mesh_filename + "/white_nuclei_mesh_%s.gii"
    print(out_mesh_filename)
    out_tex_filename = args.out_tex_filename + "/white_nuclei_tex_%s.gii"
    print(out_tex_filename)
    subjects = sorted([os.path.basename(element) for element in glob.glob(args.nuclei_dir + "/*")])
    print(subjects)

    # Labels from the volbrain nomenclature.
    # The ventricules are not used.
    nuclei = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

# TODO: add the palette in the header of the parcellation
#cols = [
#    [1., 0., 0., 1.],
#    [0., 1., 0., 1.],
#    [0., 0., 1., 1.],
#    [1., 1., 0., 1.],
#    [1., 0., 1., 1.],
#    [0., 1., 1., 1.],
#    [1., 1., 1., 1.],
#    [0.5, 0., 0., 1.],
#    [0., 0.5, 0., 1.],
#    [0., 0., 0.5, 1.],
#    [0.5, 0.5, 0., 1.],
#    [0.5, 0., 0.5, 1.],
#    [0., 0.5, 0.5, 1.],
#    [0.5, 0.5, 0.5, 1.],
#    [0.5, 1., 0., 1.],
#    [0.5, 0., 1., 1.],
#]

    # Merge all meshes and create the output parcellation for each subject.
    for subject in subjects:
        tex = aims.read(
            os.path.join(texture_dir % subject, "bh.r.aparc.annot.gii"))
        print(tex)
        label_max = numpy.max(tex[0].arraydata())
        new_tex = numpy.array(tex[0].arraydata(), copy=True)

        wmesh_file = os.path.join(
            white_dir % subject, "bh.r.aims.white.gii")
        print(wmesh_file)
        mesh_files = [
            os.path.join(nuclei_dir % subject, "wnuclei_atlas_%d_0.mesh" % i)
            for i in nuclei]
        print(mesh_files)
        out_mesh = aims.read(wmesh_file)
        for i, meshf in enumerate(mesh_files):
            mesh = aims.read(meshf)
            nb_vert = len(mesh.vertex())
            aims.SurfaceManip.meshMerge(out_mesh, mesh)
            new_tex = numpy.concatenate(
                [new_tex, numpy.ones((nb_vert,)) * (label_max + i + 1)])
        # TODO: add the palette in the header of the parcellation
        #        tex.header()["GIFTI_labels_table"][label_max + i] = {
        #            'Label' : 'nucleus_%d' % i,
        #            'RGB' : cols[i]
        #          }
        tex[0].assign(new_tex)
        del tex.header()["vertex_number"]
        aims.write(out_mesh, out_mesh_filename % subject)
        aims.write(tex, out_tex_filename % subject)


#----------------------------Main program--------------------------------------


if __name__ == "__main__":
    main()

