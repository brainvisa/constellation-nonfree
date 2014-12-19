#!/usr/bin/env python

# Aims module
from soma import aims

# system module
import optparse
import sys
try:
    from constel.lib.texturetools import clean_gyri_texture, find_wrong_labels
except:
    pass


def validation():
    try:
        from constel.lib.texturetools import clean_gyri_texture, find_wrong_labels
    except:
        raise ValidationError('constellation module is not here.')


def parseOpts(argv):
    description = 'Threshold an aimsTimeTexture' \
                  '(Warning: with one time step only)'
    parser = optparse.OptionParser(description)
    parser.add_option('-i', '--itex', dest='tex', metavar='FILE',
                      help='input gyri texture')
    parser.add_option('-m', '--mesh', dest='mesh', metavar='FILE',
                      help='input mesh')
    parser.add_option('-o', '--otex', dest='otex', metavar='FILE',
                      action='store', default=None,
                      help='output gyri texture with one connected component' \
                           'per gyrus' \
                           '(default : don\'t store resulting texture)')
    return parser, parser.parse_args(argv)


def main():
    parser, (options, args) = parseOpts(sys.argv)
    tex = aims.read(options.tex)
    mesh = aims.read(options.mesh)
    cont = True
    count = 0
    while cont and count < 10:
        print "clean up number ", str(count + 1), " :"
        tex = clean_gyri_texture(mesh, tex)
        wrong_labels = find_wrong_labels(mesh, tex)
        cont = False
        count += 1
        if len(wrong_labels) > 0:
            cont = True
    aims.write(tex, options.otex)

if __name__ == "__main__":
    main()