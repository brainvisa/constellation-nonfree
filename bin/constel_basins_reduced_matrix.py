# python system module
from __future__ import absolute_import
from __future__ import print_function
import numpy
import argparse

# aims module
from soma import aims
from constel.lib.utils.filetools import select_ROI_name

def parse_args(argv):
    """Parses the given list of arguments."""

    # creating a parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            -------------------------------------------------------------------
            Calculate an individual reduced matrix from gyri-based matrices.
            The obtained reduced matrix represents connectivity from base
            parcellation to Constellation watershed basins.
            -------------------------------------------------------------------
            """))

    # adding arguments
    parser.add_argument("reduced_matrices",
                        type=str,
                        help="List of gyrus reduced matrices of a subject."
                             "Size if (nb_vertex_in_gyrus, nb_basins)"
                             "Format: nii.gz")

    parser.add_argument("regions",
                        type=str,
                        help="List of gyrus corresponding to the matrices")

    parser.add_argument('reduced_matrix_out',
                        type=str,
                        help="Path of the output reduced matrix"
                             "Order of previous matrices will be kept")

    # parsing arguments
    return parser, parser.parse_args(argv)


def main():
    """Execute the command to calculate a connectome and its snapshot.
    """
    # Load the arguments of parser (ignore script name: sys.arg[0])
    arguments = sys.argv[1:]
    parser, args = parse_args(arguments)

    reduced_matrices = eval(args.reduced_matrices)
    regions = eval(args.regions)

    gyrus_reduced_matrices = []

    # Import and fully reduce the matrices
    for matrix_path in reduced_matrices:
        vol = aims.read(matrix_path)
        matrix = np.asarray(vol)[:, :, 0, 0]

        if matrix.shape[0] < matrix.shape[1]:
            matrix = np.transpose(matrix)

        gyrus_reduced_matrix = np.sum(matrix, axis=0)
        gyrus_reduced_matrices.append(gyrus_reduced_matrix)

    # Create the (base parcellation to basins) reduced matrix
    s = 0
    for mat in gyrus_reduced_matrices:
        s += mat.size

    reduced_matrix = np.zeros((len(gyrus_reduced_matrix), s))

    row = 0
    column_offset = 0
    for mat in gyrus_reduced_matrices:
        for i in range(mat.size):
            reduced_matrix[row][i + column_offset] = mat[i]
        column_offset += rmat.size
        row += 1

    # Write reduced matrix on disk
    aims.write(reduced_matrix, args.reduced_matrix_out)
