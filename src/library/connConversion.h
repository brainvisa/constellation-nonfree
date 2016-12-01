#ifndef CONSTELLATION_CONNCONVERSION_H
#define CONSTELLATION_CONNCONVERSION_H

#include <aims/sparsematrix/sparseMatrix.h>
#include <constellation/connectivities.h>

namespace constel {

  aims::SparseMatrix *connectivitiesToSparseMatrix(
      const Connectivities &conn);

  Connectivities *sparseMatrixToConnectivities(
      const aims::SparseMatrix &mat);

  void sparseMatrixToConnectivities(
      const aims::SparseMatrix &mat, Connectivities &conn);

  void sparseMatrixToConnectivities(
      const aims::boost_sparse_matrix &mat, Connectivities &conn);

  void writeConnectivities(
      const Connectivities &conn, const std::string &filename,
      bool ascii=false);

}
#endif

