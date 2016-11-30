#include <constellation/connConversion.h>

using namespace aims;
using namespace std;


namespace constel {

  //--------------------------------
  //  connectivitiesToSparseMatrix
  //--------------------------------
  aims::SparseMatrix *connectivitiesToSparseMatrix(
      const Connectivities &conn) {
    size_t size1 = conn.size();
    size_t size2 = conn[0].size();

    SparseMatrix *mp = new SparseMatrix(size1,size2);
    SparseMatrix & m = *mp;
    Connectivities::const_iterator il, el = conn.end();
    size_t i, j;
    Connectivity::sparse_const_iterator ic, ec;
    for (il=conn.begin(), i=0; il!=el; ++il, ++i) {
      for (ic=il->sparse_begin(), ec=il->sparse_end(); ic!=ec; ++ic) {
        j = ic->first;
        double val = ic->second;
        if( val!= 0) {
          m(i,j) = val;
        }
      }
    }
    return mp;
  }

  //--------------------------------
  //  sparseMatrixToConnectivities
  //--------------------------------
  void sparseMatrixToConnectivities(
      const aims::boost_sparse_matrix &mat, Connectivities &conn) {
    conn.clear();
    conn.resize(mat.size1(), Connectivity(mat.size2()));
    boost_sparse_matrix::const_iterator1 il, el = mat.end1();
    boost_sparse_matrix::const_iterator2 ic, ec;
    unsigned i;
    for (il = mat.begin1(); il != el; ++il) {
      i = il.index1();
      for (ic = il.begin(), ec = il.end(); ic != ec; ++ic) {
        conn[i][ic.index2()] = *ic;
      }
    }
  }

  //--------------------------------
  //  sparseMatrixToConnectivities
  //--------------------------------
  void sparseMatrixToConnectivities(
      const aims::SparseMatrix &mat, Connectivities &conn) {
    conn.clear();
    conn.resize(mat.getSize1(), Connectivity(mat.getSize2()));
    SparseMatrix::const_iterator1 il, el = mat.end1();
    SparseMatrix::const_iterator2 ic, ec;
    unsigned i;
    for (il = mat.begin1(); il != el; ++il) {
      i = il.index1();
      for (ic = il.begin(), ec = il.end(); ic != ec; ++ic) {
        conn[i][ic.index2()] = *ic;
      }
    }
  }

  //--------------------------------
  //  sparseMatrixToConnectivities
  //--------------------------------
  Connectivities *sparseMatrixToConnectivities(
      const aims::SparseMatrix &mat) {
    Connectivities* conn = new Connectivities;
    sparseMatrixToConnectivities(mat, *conn);
    return conn;
  }

  //--------------------------------
  //  sparseMatrixToConnectivities
  //--------------------------------
  void writeConnectivities(
      const Connectivities &conn, const string &filename, bool ascii) {
    aims::SparseMatrix *m = connectivitiesToSparseMatrix(conn);
    try {
      if (ascii) m->write(filename, "ascii");
      else m->write(filename);
    }
    catch( ... ) {
      delete m;
      throw;
    }
    delete m;
  }

} // namespace constel

