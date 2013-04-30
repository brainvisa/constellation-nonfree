from soma import aims
import numpy
import sys
from constel.lib import texturetools

def usage():
  print "Create average texture labels"
  print "usage : python average_texture_labels.py output.gii(.tex) subject1.gii(.tex) ... subjectN.gii(.tex)"


if __name__ == "__main__":
    if len(sys.argv)<=3:
        usage()
        sys.exit(1)

    li_labels = []
    output = sys.argv[1]
    inputs = sys.argv[2:]
    texturetools.average_texture_labels(output, inputs)