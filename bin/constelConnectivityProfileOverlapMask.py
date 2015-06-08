#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

# soma module
from soma import aims

# system module
import optparse
import numpy
import sys


def parseOpts(argv):
    desc="""Creation of a mask : while it is estimated that there are at least
    "fibersnbmin" connexions on one gyrus mean connectivity profile
    (for more than half of the subjects or equivalent )"""

    parser = optparse.OptionParser(desc)
    sys.stdout.flush()
    parser.add_option('-p', '--profiles', dest='connectivity_profile',
                      action='append', help='connectivity_profile' )
    parser.add_option('-o', '--output', dest='mask', help='mask')

    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)

    files = options.connectivity_profile

    nbOfSubjects_threshold = int(round(len(files)/2))
    tex = aims.TimeTexture_S16()
    vertex_nb = 0
    subjects_count = 0
    final_valid_vertex_sum = numpy.ones((0))
    listOfTextures= []
    for texture in files:
        listOfTextures = aims.read(texture)
        listOfTextures_array = listOfTextures[0].arraydata()
        threshold_value = listOfTextures_array.max() / 100. #numpy.sum( listOfTextures_array ) / 1500.
        print threshold_value
        print listOfTextures_array
        if subjects_count == 0:
            vertex_nb = listOfTextures[0].nItem()
            tex[0].reserve(vertex_nb)
            for i in xrange(vertex_nb):
                tex[0].push_back(0)
            final_valid_vertex_sum = tex[0].arraydata()
        valid_vertex = listOfTextures_array >= threshold_value
        print valid_vertex
        final_valid_vertex_sum[valid_vertex]+=1
        subjects_count += 1
    #context.write( 'type final_valid_vertex_sum:', final_valid_vertex_sum.dtype, ', min/max:', final_valid_vertex_sum.min(),  final_valid_vertex_sum.max() )
    #context.write( 'size:', final_valid_vertex_sum.shape )
    texsum= aims.TimeTexture('S16')
    t = texsum[0]
    t.reserve(final_valid_vertex_sum.shape[0])
    for i in xrange(final_valid_vertex_sum.shape[0]):
        t.push_back(int(final_valid_vertex_sum[i]))
    # tex[0].assign( final_valid_vertex_sum )
    #context.write( 'tex min/max:', texsum[0].arraydata().min(), texsum[0].arraydata().max() )
    # aims.write(texsum, '/tmp/final_valid_vertex_sum.gii')

    final_valid_vertex_sum[numpy.where(final_valid_vertex_sum < nbOfSubjects_threshold)] = 0
    final_valid_vertex_sum[numpy.where(final_valid_vertex_sum >= nbOfSubjects_threshold)] = 1
    aims.write(tex, options.mask)

if __name__ == "__main__":
    main()
