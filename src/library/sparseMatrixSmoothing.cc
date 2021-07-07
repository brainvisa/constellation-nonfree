
// to enable shallow_array_adaptor in boost
#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <constellation/sparseMatrixSmoothing.h>
#include <constellation/connConversion.h>
#include <aims/mesh/curv.h>
#include <aims/sparsematrix/sparseordensematrix.h>
#include <cartobase/smart/rcptrtrick.h>
#include <boost/numeric/ublas/matrix.hpp> // zero_matrix is needed in operation
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <aims/utility/stl_conversion.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/utility/converter_texture.h>
#include <float.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;
using namespace boost;


namespace
{

  /* MatrixTraits and MatrixProxy: allows common API for SparseMatrix and
     Connectivities
  */

  template <typename Matrix>
  class MatrixTraits {
   public:
    typedef typename Matrix::iterator iterator;
    typedef typename Matrix::iterator columniterator;
  };

  template <>
  class MatrixTraits<SparseMatrix> {
   public:
    typedef SparseMatrix::iterator1 iterator;
    typedef SparseMatrix::iterator2 columniterator;
  };

  template <>
  class MatrixTraits<Connectivities> {
   public:
    typedef Connectivities::iterator iterator;
    typedef constel::Connectivity::iterator columniterator;
  };

  template <typename Matrix>
  class MatrixProxy {
   public:
    typedef boost::numeric::ublas::shallow_array_adaptor<double> dstorage;
    typedef boost::numeric::ublas::matrix<
      double, boost::numeric::ublas::row_major, dstorage> boost_matrix;

    MatrixProxy(rc_ptr<Matrix> matrix);
    size_t nlines() const;
    size_t ncols() const;
    rc_ptr<SparseMatrix> toSparseMatrix();
    SparseOrDenseMatrix::DenseMatrixType toDenseMatrix();
    void fromMatrix(rc_ptr<SparseMatrix> matrix);
    void fromMatrix(SparseOrDenseMatrix::DenseMatrixType matrix);
    void fromMatrix(boost_sparse_matrix &matrix);
    void fromMatrix(boost_matrix &matrix);
    void clear();

   private:
    rc_ptr<Matrix> _matrix;
  };

  template <typename Matrix>
  inline
  MatrixProxy<Matrix>::MatrixProxy(rc_ptr<Matrix> matrix) : _matrix(matrix) {}

  // SparseMatrix
  template <>
  inline
  size_t MatrixProxy<SparseMatrix>::nlines() const {
    return _matrix->getSize1();
  }

  template <>
  inline
  size_t MatrixProxy<SparseMatrix>::ncols() const {
    return _matrix->getSize2();
  }

  template <>
  inline
  rc_ptr<SparseMatrix> MatrixProxy<SparseMatrix>::toSparseMatrix() {
    return _matrix;
  }

  template <>
  inline
  SparseOrDenseMatrix::DenseMatrixType
  MatrixProxy<SparseMatrix>::toDenseMatrix() {
    SparseOrDenseMatrix sd;
    sd.setMatrix(_matrix);
    return sd.asDense();
  }

  template <>
  inline
  void MatrixProxy<SparseMatrix>::fromMatrix(rc_ptr<SparseMatrix> matrix) {
    _matrix = matrix;
  }

  template <>
  inline
  void MatrixProxy<SparseMatrix>::fromMatrix(
      SparseOrDenseMatrix::DenseMatrixType matrix) {
    SparseOrDenseMatrix sd;
    sd.setMatrix(matrix);
    _matrix = sd.asSparse();
  }

  template <>
  inline
  void MatrixProxy<SparseMatrix>::fromMatrix(boost_sparse_matrix &matrix) {
    if (&matrix != &_matrix->boostMatrix()) {
      _matrix->boostMatrix() = matrix;
    }
  }

  template <>
  inline
  void MatrixProxy<SparseMatrix>::fromMatrix(boost_matrix &matrix) {
    _matrix->boostMatrix() = matrix;
  }

  template <>
  inline
  void MatrixProxy<SparseMatrix>::clear() {
    _matrix->boostMatrix().clear();
  }

  // SparseOrDenseMatrix
  template <>
  inline
  size_t MatrixProxy<SparseOrDenseMatrix>::nlines() const {
    return _matrix->getSize1();
  }

  template <>
  inline
  size_t MatrixProxy<SparseOrDenseMatrix>::ncols() const {
    return _matrix->getSize2();
  }

  template <>
  inline
  rc_ptr<SparseMatrix> MatrixProxy<SparseOrDenseMatrix>::toSparseMatrix() {
    return _matrix->asSparse();
  }

  template <>
  inline
  SparseOrDenseMatrix::DenseMatrixType
  MatrixProxy<SparseOrDenseMatrix>::toDenseMatrix() {
    return _matrix->asDense();
  }

  template <>
  inline
  void MatrixProxy<SparseOrDenseMatrix>::fromMatrix(
      rc_ptr<SparseMatrix> matrix) {
    _matrix->setMatrix( matrix );
  }

  template <>
  inline
  void MatrixProxy<SparseOrDenseMatrix>::fromMatrix(
      SparseOrDenseMatrix::DenseMatrixType matrix) {
    _matrix->setMatrix( matrix );
  }

  template <>
  inline
  void MatrixProxy<SparseOrDenseMatrix>::fromMatrix(
      boost_sparse_matrix &matrix) {
    if (!_matrix->isDense()
        && &_matrix->sparseMatrix()->boostMatrix() == &matrix)
      return;
    _matrix->setMatrix(rc_ptr<SparseMatrix>(new SparseMatrix));
    _matrix->sparseMatrix()->boostMatrix() = matrix;
  }

  template <>
  inline
  void MatrixProxy<SparseOrDenseMatrix>::fromMatrix(boost_matrix &matrix) {
    if (_matrix->isDense()
        && &_matrix->denseMatrix()->at(0) == &matrix(0, 0))
      return;
    if (!_matrix->isDense())
      _matrix->reallocate(matrix.size1(), matrix.size2());
    memcpy(&_matrix->denseMatrix()->at(0),
           &matrix, matrix.size1() * matrix.size2() * sizeof(double));
  }

  template <>
  inline
  void MatrixProxy<SparseOrDenseMatrix>::clear() {
    _matrix->setMatrix(rc_ptr<SparseMatrix>(new SparseMatrix));
  }

  // Connectivities
  template <>
  inline
  size_t MatrixProxy<Connectivities>::nlines() const {
    return _matrix->size();
  }

  template <>
  inline
  size_t MatrixProxy<Connectivities>::ncols() const {
    return _matrix->begin()->size();
  }

  template <>
  inline
  rc_ptr<SparseMatrix> MatrixProxy<Connectivities>::toSparseMatrix() {
    return rc_ptr<SparseMatrix>(connectivitiesToSparseMatrix(*_matrix));
  }

  template <>
  inline
  SparseOrDenseMatrix::DenseMatrixType
  MatrixProxy<Connectivities>::toDenseMatrix() {
    SparseOrDenseMatrix sd;
    sd.setMatrix(rc_ptr<SparseMatrix>(connectivitiesToSparseMatrix(*_matrix)));
    return sd.asDense(); // FIXME: inefficient, double conversion
  }

  template <>
  inline
  void MatrixProxy<Connectivities>::fromMatrix(rc_ptr<SparseMatrix> matrix) {
    sparseMatrixToConnectivities(*matrix, *_matrix);
  }

  template <>
  inline
  void MatrixProxy<Connectivities>::fromMatrix(
    SparseOrDenseMatrix::DenseMatrixType matrix) {
    SparseOrDenseMatrix sd;
    sd.setMatrix(matrix); // FIXME: inefficient, double conversion
    sparseMatrixToConnectivities(*sd.asSparse(), *_matrix);
  }

  template <>
  inline
  void MatrixProxy<Connectivities>::fromMatrix(boost_sparse_matrix &matrix) {
    sparseMatrixToConnectivities(matrix, *_matrix);
  }

  template <>
  inline
  void MatrixProxy<Connectivities>::fromMatrix(boost_matrix & /* matrix */) {
    // FIXME
    cout << "MatrixProxy<Connectivities>::fromMatrix( boost_matrix &matrix )\n\
            - NOT IMPLEMENTED\n";
    // sparseMatrixToConnectivities( matrix, *_matrix );
  }

  template <>
  inline
  void MatrixProxy<Connectivities>::clear() {
    _matrix->clear();
  }

  /*
  size_t countElements(const boost_sparse_matrix & matrix) {
    size_t count = 0;
    boost_sparse_matrix::const_iterator1 il, el = matrix.end1();
    boost_sparse_matrix::const_iterator2 ic, ec;
    for(il=matrix.begin1(); il!=el; ++il)
      for(ic=il.begin(), ec=il.end(); ic!=ec; ++ic) {
        ++count;
      }
    cout << "countElements in mat " << matrix.size1() << " x "
      << matrix.size2() << ": " << count << endl;
    return count;
  }
  */

  boost_sparse_matrix* subMatrixLaplacian(
      const LaplacianWeights & matrix, const vector<int32_t> & items,
      int32_t patch) {
    size_t i, j, n = items.size();
    vector<size_t> index(items.size(), 0);

    for (i=0, j=0; i<n; ++i)
      if (items[i] == patch)
        index[i] = j++;
    n = j;

    boost_sparse_matrix * res = new boost_sparse_matrix(n, n);
    boost_sparse_matrix & rmat = *res;
    LaplacianWeights::const_iterator il, el = matrix.end();
    size_t line;
    // size_t count = 0;
    double weight;

    for (il=matrix.begin(), i=0; il!=el; ++il) {
      line = il->first;
      if (items[line] == patch) {
        set<pair<unsigned,float> >::const_iterator
          ic = il->second.begin(), ec = il->second.end();
        size_t col;
        weight = 0;
        for ( ; ic!=ec; ++ic) {
          col = ic->first;
          if (items[col] != patch) { // not in patch, get weight for correction
            weight += ic->second;
          } else {
            rmat(i, index[col]) = ic->second;
            // ++count;
          }
        }
        // weight has been dropped and has to be regained
        rmat(i, i) += weight;
        ++i;
      }
    }
    // cout << "subMatrixLaplacian items: " << count << endl;
    return res;
  }

/*
  boost_sparse_matrix* subMatrix(
      const LaplacianWeights & matrix, const vector<size_t> & items) {
    boost_sparse_matrix * res
      = new boost_sparse_matrix(items.size(), items.size());
    boost_sparse_matrix & rmat = *res;
    LaplacianWeights::const_iterator il, el = matrix.end();
    size_t i, j, n = items.size();
    size_t line;
    size_t count = 0;
    for (i=0; i<n; ++i) {
      line = items[i];
      il = matrix.find(line);
      if (il != el) {
        set<pair<unsigned,float> >::const_iterator
          ic = il->second.begin(), ec = il->second.end();
        size_t col;
        for (j=0; j!=n; ++j) {
          col = items[j];
          while (ic != ec && ic->first < col)
            ++ic;
          if (ic != ec && ic->first == col) {
            rmat(i, j) = ic->second;
            ++count;
          }
        }
      }
    }
    cout << "subMatrix items: " << count << endl;
    return res;
  }
*/

  boost_sparse_matrix* toTransposedBoostSparseMatrix(
      LaplacianWeights & matrix) {
    boost_sparse_matrix * res
      = new boost_sparse_matrix(matrix.size(), matrix.size());
    boost_sparse_matrix & rmat = *res;
    LaplacianWeights::const_iterator il, el = matrix.end();
    size_t line;
    // size_t count = 0;
    for (il=matrix.begin(); il!=el; ++il) {
      line = il->first;
      set<pair<unsigned,float> >::const_iterator
        ic, ec = il->second.end();
      for (ic=il->second.begin(); ic!=ec; ++ic) {
        rmat(ic->first, line) = ic->second;
        // ++count;
      }
    }
    // cout << "boost_sparse_matrix items: " << count << endl;
    return res;
  }

  //
  template <typename Matrix>
  void _sparseMatrixDiffusionSmoothing(
      rc_ptr<Matrix> matrix,
      const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
      const TimeTexture<int32_t> & patches, int32_t patch) {
    cout << "Sparse Matrix Diffusion Smoothing: weights calculation...\n";
    MatrixProxy<Matrix> matp(matrix);
    double wmax = 0.98; // arbitrary, as in AimsTextureSmoothing
    // Compute the weights
    std::map<unsigned, std::set< std::pair<unsigned,float> > >
      weightLapl = AimsMeshWeightFiniteElementLaplacian(mesh.begin()->second,
                                                       wmax);

    unsigned nl = matp.nlines(), nc = matp.ncols();
    // compute the diffusion duration equivalent to a gaussian sigma
    // (probably a standard formula, I just found it here)
    double dur = sigma * sigma / (8. * ::log(2));
    cout << "dur: " << dur << endl;
    float dt = 0.001; // FIXME
    unsigned niter = unsigned(rint(dur/dt));
    vector<double> outrow;
    cout << "laplacian smoothing... niter: " << niter << ", dt: " << dt
      << endl;
    cout << "lines: " << nl << ", cols: " << nc << endl;

    LaplacianWeights* laplmat
      = makeLaplacianSmoothingCoefficients(weightLapl, niter, dt, wthresh);

    cout << "Laplacian coefs done.\n";

    /* Smoothing along lines is done using the coefs in laplmat.
      More precisely, the transposed of laplmat should be applied to the
      input matrix ( matrix * tr( laplmat ) ) as matricial operation.

      Now we also want to smooth along columns, inside the starting points
      region.
      Laplacian coefficients are corrected according to vertices outside the
      patch also having weights in the laplacian matrix.
      The sub-matrix of laplmat corresponding to the region vertices is
      extracted using this correction.
      So, we take the sub-matrix submat, and apply it to the columns
      of martrix: submat * matrix.

      Combining the two operations, we have:
      submat * martrix * tr( laplmat )
    */

    boost_sparse_matrix *submat = subMatrixLaplacian(*laplmat,
      patches.begin()->second.data(), patch);
    cout << "Submatrix done. Getting SparseMatrix\n";
    rc_ptr<SparseMatrix> smat = matp.toSparseMatrix();
    cout << "Applying for columns...\n";

    // typename MatrixProxy<Matrix>::boost_matrix colsmoothed( nl, nc );
    boost_sparse_matrix colsmoothed(nl, nc);
    // smooth along cols
    sparse_prod(*submat, smat->boostMatrix(), colsmoothed);
    delete submat;
    cout << "getting SparseMatrix for coefs\n";
    boost_sparse_matrix* tlapl = toTransposedBoostSparseMatrix(*laplmat);
    delete laplmat;
    cout << "Applying for lines...\n";
    smat->boostMatrix().clear();
    // smooth along lines
    smat.reset(0); // release smat
    matp.clear(); // release input matrix
    SparseOrDenseMatrix::DenseMatrixType res(nc, nl);
    typename MatrixProxy<Matrix>::dstorage darray(
      res->getSizeY() * res->getSizeX(), &res->at(0));
    typename MatrixProxy<Matrix>::boost_matrix
      mdata(res->getSizeY(), res->getSizeX(), darray);
    sparse_prod(colsmoothed, *tlapl, mdata);
    cout << "done.\n";
    delete tlapl;
    colsmoothed.clear();
    // back to input matrix
    matp.fromMatrix(res);
    return;
  }
}


namespace constel
{

  void sparseMatrixDiffusionSmoothing(
      rc_ptr<SparseMatrix> matrix,
      const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
      const TimeTexture<int32_t> & patches, int32_t patch)
  {
    _sparseMatrixDiffusionSmoothing(matrix, mesh, wthresh, sigma, patches,
                                    patch);
  }

  void sparseMatrixDiffusionSmoothing(
      rc_ptr<SparseMatrix> matrix,
      const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
      const TimeTexture<int16_t> & patches, int32_t patch)
  {
    TimeTexture<int32_t> t32;
    Converter<TimeTexture<int16_t>, TimeTexture<int32_t> > c;
    c.convert(patches, t32);
    _sparseMatrixDiffusionSmoothing(matrix, mesh, wthresh, sigma, t32,
                                    patch);
  }

  void sparseMatrixDiffusionSmoothing(
      rc_ptr<SparseOrDenseMatrix> matrix,
      const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
      const TimeTexture<int32_t> & patches, int32_t patch)
  {
    _sparseMatrixDiffusionSmoothing( matrix, mesh, wthresh, sigma, patches,
                                     patch );
  }

  void sparseMatrixDiffusionSmoothing(
      rc_ptr<SparseOrDenseMatrix> matrix,
      const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
      const TimeTexture<int16_t> & patches, int32_t patch)
  {
    TimeTexture<int32_t> t32;
    Converter<TimeTexture<int16_t>, TimeTexture<int32_t> > c;
    c.convert(patches, t32);
    _sparseMatrixDiffusionSmoothing(matrix, mesh, wthresh, sigma, t32,
                                    patch);
  }

  void sparseMatrixDiffusionSmoothing(
      rc_ptr<Connectivities> conn_ptr,
      const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
      const TimeTexture<int32_t> & patches, int32_t patch)
  {
    _sparseMatrixDiffusionSmoothing(conn_ptr, mesh, wthresh, sigma,
                                    patches, patch);
  }

  void sparseMatrixDiffusionSmoothing(
      rc_ptr<Connectivities> conn_ptr,
      const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
      const TimeTexture<int16_t> & patches, int32_t patch)
  {
    TimeTexture<int32_t> t32;
    Converter<TimeTexture<int16_t>, TimeTexture<int32_t> > c;
    c.convert(patches, t32);
    _sparseMatrixDiffusionSmoothing(conn_ptr, mesh, wthresh, sigma,
                                    t32, patch);
  }

  void sparseMatrixGaussianSmoothing(
      SparseMatrix & matrix, const AimsSurfaceTriangle & inAimsMesh,
      float distthresh, float /* wthresh */)
  {
    /*
    Smoothing of a connectivity matrix according to the aims mesh and to the
    neighbourhood distance (distthresh), threshold of the resulting matrix by
    wthresh
    */

    // Comuting geomap : neighborhood map
    //if distthresh!= 0
    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity
    double square_sigma = distthresh*distthresh;
    std::cout << "Computing geomap..." << std::flush;
    size_t nv = inAimsMesh.vertex().size();

    std::vector< std::map<size_t, float> > res;

    meshdistance::pairwiseDistanceMaps( inAimsMesh.begin()->second, res,
                                        distthresh );
    std::cout << "OK" << std::endl;
    
    //Sparse matrix smoothing
    SparseMatrix * smoothed_matrix_ptr = new SparseMatrix(matrix.getSize1(),
                                                          matrix.getSize2());
    SparseMatrix & smoothed_matrix = *smoothed_matrix_ptr;
    SparseMatrix::iterator1 s1;
    SparseMatrix::iterator2 s2;
    if( ( (int32_t) nv == matrix.getSize1() )
        && ( (int32_t) nv == matrix.getSize2() ) )
    {
      for (s1 = matrix.begin1(); s1 != matrix.end1(); s1++)
      {
        std::size_t A = s1.index1();
        for (s2 = s1.begin(); s2 != s1.end(); s2++)
        {
          std::size_t B = s2.index2();
          double w = 0;
          for (map<size_t, float>::const_iterator B_neigh = res[B].begin();
               B_neigh != res[B].end(); ++B_neigh)
          {
            for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
                 A_neigh != res[A].end(); ++A_neigh)
            {
              double e1 = square(A_neigh->second);
              double e2 = square(B_neigh->second);
              //"smoothing coefficient"
              w = (*s2)*std::exp(-(e1 + e2) / (2*distthresh*distthresh))
                  /(square_sigma * two_pi);
              if (w > G_THRESH)
              {
                smoothed_matrix(A_neigh->first,B_neigh->first)+=w;
              }
            }
          }
        }
      }
    }
    else if ((int32_t) nv == matrix.getSize2())
    {
      for (s2 = matrix.begin2(); s2 != matrix.end2(); s2++)
      {
        std::size_t A = s2.index2();
        for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
             A_neigh != res[A].end(); ++A_neigh)
        {
          double e1 = square(A_neigh->second);
          
          double smooth_coef = std::exp(-e1  / ( 2*distthresh*distthresh));
          for (s1 = s2.begin(); s1 != s2.end(); s1++)
          {
            std::size_t B = s1.index1();
            float w = (*s1)*smooth_coef;
            if (w > G_THRESH)
            {
              smoothed_matrix(B,A_neigh->first) += w;
            }
          }
        }
      }
    }
    matrix.setZero();
    matrix = smoothed_matrix;
    delete smoothed_matrix_ptr;

  }// sparseMatrixGaussianSmoothing


  void sparseMatrixGaussianSmoothingNormed(
      SparseMatrix & matrix, const AimsSurfaceTriangle & inAimsMesh,
      float distthresh, float /* wthresh */) {
    /*
    Smoothing of a connectivity matrix according to the aims mesh and to the
    neighbourhood distance (distthresh), threshold of the resulting matrix by
    wthresh. Each connection contribution is normalized to 1
    */

    // Comuting geomap : neighborhood map
    //if distthresh!= 0
    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity
    double square_sigma = distthresh*distthresh;
    std::cout << "Computing geomap..." << std::flush;
    size_t nv = inAimsMesh.vertex().size();

    std::vector< std::map<size_t, float> > res;

    meshdistance::pairwiseDistanceMaps( inAimsMesh.begin()->second, res,
                                          distthresh );
    std::cout << "OK" << std::endl;
    
    //Sparse matrix smoothing
    SparseMatrix * smoothed_matrix_ptr = new SparseMatrix(matrix.getSize1(),
                                                          matrix.getSize2());
    SparseMatrix & smoothed_matrix = *smoothed_matrix_ptr;
    SparseMatrix::iterator1 s1;
    SparseMatrix::iterator2 s2;
    if (((int32_t) nv == matrix.getSize1())
        && ((int32_t) nv == matrix.getSize2()))
    {
      for (s1 = matrix.begin1(); s1 != matrix.end1(); s1++)
      {
        std::size_t A = s1.index1();
        for (s2 = s1.begin(); s2 != s1.end(); s2++)
        {
          std::size_t B = s2.index2();
          double w = 0;
          for (map<size_t, float>::const_iterator B_neigh = res[B].begin();
               B_neigh != res[B].end(); ++B_neigh)
          {
            for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
                 A_neigh != res[A].end(); ++A_neigh)
            {
              double e1 = square(A_neigh->second);
              double e2 = square(B_neigh->second);
              //"smoothing coefficient"
              w = (*s2)*std::exp(-(e1 + e2)  / (2*distthresh*distthresh ))
                  /(square_sigma * two_pi);
              if (w > G_THRESH)
              {
                smoothed_matrix(A_neigh->first,B_neigh->first)+=w;
              }
            }
          }
        }
      }
    }
    else if ((int32_t) nv == matrix.getSize2())
    {
      for (s2 = matrix.begin2(); s2 != matrix.end2(); s2++)
      {
        size_t A = s2.index2();
        size_t a = res[A].size();
        std::vector<float> weightsVector(a);
        for (size_t i = 0; i <a; i++)
        {
          weightsVector[i]=0;
        }
//         std::cout << weightsVector[0] << std::endl;
        for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
             A_neigh != res[A].end(); ++A_neigh)
        {
          double e1 = square(A_neigh->second);
          
          double smooth_coef = std::exp(-e1 / (2*distthresh*distthresh));
//           /(square_sigma * two_pi);//"smoothing coefficient"
          for (s1 = s2.begin(); s1 != s2.end(); s1++)
          {
            std::size_t B = s1.index1();
            float w = (*s1)*smooth_coef;
            if (w > G_THRESH)
            {
              smoothed_matrix(B,A_neigh->first)+=w;
            }
          }
        }
        for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
             A_neigh != res[A].end(); ++A_neigh)
        {
          double e1 = square(A_neigh->second);
          
          double smooth_coef = std::exp(-e1 / (2*distthresh*distthresh));
//           /(square_sigma * two_pi);//"smoothing coefficient"
          for (s1 = s2.begin(); s1 != s2.end(); s1++)
          {
            std::size_t B = s1.index1();
            float w = (*s1)*smooth_coef;
            if (w > G_THRESH)
            {
              smoothed_matrix(B,A_neigh->first)+=w;
            }
          }
        }
      }
    }
    matrix.setZero();
    matrix = smoothed_matrix;
    delete smoothed_matrix_ptr;
  
  }// sparseMatrixGaussianSmoothingNormed

} // namespace constel
