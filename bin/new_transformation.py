from soma import aims
import optparse

def parseOpts(argv):
  desc="""transformation between landmark native and normalized landmark"""

  parser = optparse.OptionParser( desc )
  parser.add_option( '-i', '--img', dest='rawT1image', help='native rawT1image' )
  parser.add_option( '-n', '--normimg', dest='norm_rawT1image', help='normalised rawT1image' )
  parser.add_option( '-t', '--transfo', dest='dw_to_t1', help='dw_to_t1 emerging from connectomist' )
  parser.add_option( '-o', '--output', dest='transformation', help='transformation' )

  return parser, parser.parse_args(argv)

def main():
  parser, (options, args) = parseOpts(sys.argv)

  f = aims.Finder()

  f.check( options.rawT1image )
  hd_nat = f.header

  f.check( options.norm_rawT1image )
  hd_nrom = f.header

  tf_nat = aims.AffineTransformation3d( hd_nat[ 'transformations' ][-1] )
  tf_norm = aims.AffineTransformation3d( hd_norm[ 'transformations' ][-1] )
  tf_connect = aims.read( options.dw_to_t1 )

  tf_final = tf_nat.inverse() * tf_norm * tf_connect

  aims.write( tf_final, options.transformation )

if __name__ == "__main__" : main()

  