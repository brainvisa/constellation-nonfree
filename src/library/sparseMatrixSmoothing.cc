#include <constellation/sparseMatrixSmoothing.h>
#include <aims/mesh/curv.h>
#include <float.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;


namespace
{

  /* MatrixTraits and MatrixProxy: allows common API for SparseMatrix and
     Connectivities
  */

  template <typename Matrix>
  class MatrixTraits
  {
  public:
    typedef typename Matrix::iterator iterator;
    typedef typename Matrix::iterator columniterator;
  };


  template <>
  class MatrixTraits<SparseMatrix>
  {
  public:
    typedef SparseMatrix::iterator1 iterator;
    typedef SparseMatrix::iterator2 columniterator;
  };


  template <>
  class MatrixTraits<Connectivities>
  {
  public:
    typedef Connectivities::iterator iterator;
    typedef til::SparseVector<double>::sparse_iterator columniterator;
  };


  template <typename Matrix>
  class MatrixProxy
  {
  public:
    typedef typename MatrixTraits<Matrix>::iterator iterator;
    typedef typename MatrixTraits<Matrix>::columniterator columniterator;

    MatrixProxy( Matrix & matrix );
    size_t nlines() const;
    size_t ncols() const;
    iterator begin();
    iterator end();
    static size_t column( const columniterator & );
    static double & value( const columniterator & );
    static columniterator colbegin( const iterator & );
    static columniterator colend( const iterator & );
    void assign( size_t line, size_t col, double value, iterator & it );

    Matrix & _matrix;
  };


  template <typename Matrix>
  inline
  MatrixProxy<Matrix>::MatrixProxy( Matrix & matrix )
    : _matrix( matrix )
  {
  }

  // SparseMatrix

  template <>
  inline
  size_t MatrixProxy<SparseMatrix>::nlines() const
  {
    return _matrix.getSize1();
  }


  template <>
  inline
  size_t MatrixProxy<SparseMatrix>::ncols() const
  {
    return _matrix.getSize2();
  }


  template <>
  inline
  MatrixProxy<SparseMatrix>::iterator MatrixProxy<SparseMatrix>::begin()
  {
    return _matrix.begin1();
  }


  template <>
  inline
  MatrixProxy<SparseMatrix>::iterator MatrixProxy<SparseMatrix>::end()
  {
    return _matrix.end1();
  }


  template <>
  inline
  size_t MatrixProxy<SparseMatrix>::column(
    const MatrixProxy<SparseMatrix>::columniterator & iter )
  {
    return iter.index2();
  }


  template <>
  inline
  double & MatrixProxy<SparseMatrix>::value(
    const MatrixProxy<SparseMatrix>::columniterator & iter )
  {
    return *iter;
  }


  template <>
  inline
  MatrixProxy<SparseMatrix>::columniterator
    MatrixProxy<SparseMatrix>::colbegin(
      const MatrixProxy<SparseMatrix>::iterator & iter )
  {
    return iter.begin();
  }


  template <>
  inline
  MatrixProxy<SparseMatrix>::columniterator
    MatrixProxy<SparseMatrix>::colend(
      const MatrixProxy<SparseMatrix>::iterator & iter )
  {
    return iter.end();
  }


  template <>
  inline
  void MatrixProxy<SparseMatrix>::assign( size_t line, size_t col,
    double value, MatrixProxy<SparseMatrix>::iterator & )
  {
    _matrix( line, col ) = value;
  }


  // Connectivities


  template <>
  inline
  size_t MatrixProxy<Connectivities>::nlines() const
  {
    return _matrix.size();
  }


  template <>
  inline
  size_t MatrixProxy<Connectivities>::ncols() const
  {
    return _matrix.begin()->size();
  }


  template <>
  inline
  MatrixProxy<Connectivities>::iterator MatrixProxy<Connectivities>::begin()
  {
    return _matrix.begin();
  }


  template <>
  inline
  MatrixProxy<Connectivities>::iterator MatrixProxy<Connectivities>::end()
  {
    return _matrix.end();
  }


  template <>
  inline
  size_t MatrixProxy<Connectivities>::column(
    const MatrixProxy<Connectivities>::columniterator & iter )
  {
    return iter->first;
  }


  template <>
  inline
  double & MatrixProxy<Connectivities>::value(
    const MatrixProxy<Connectivities>::columniterator & iter )
  {
    return iter->second;
  }


  template <>
  inline
  MatrixProxy<Connectivities>::columniterator
    MatrixProxy<Connectivities>::colbegin(
      const MatrixProxy<Connectivities>::iterator & iter )
  {
    return iter->sparse_begin();
  }


  template <>
  inline
  MatrixProxy<Connectivities>::columniterator
    MatrixProxy<Connectivities>::colend(
      const MatrixProxy<Connectivities>::iterator & iter )
  {
    return iter->sparse_end();
  }


  template <>
  inline
  void MatrixProxy<Connectivities>::assign( size_t, size_t col,
    double value, MatrixProxy<Connectivities>::iterator & iter )
  {
    (*iter)[ col ] = value;
  }

  //

  template <typename Matrix>
  void _sparseMatrixDiffusionSmoothing( Matrix & matrix,
    const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma )
  {
    cout << "sparseMatrixDiffusionSmoothing. Weights calculation...\n";
    MatrixProxy<Matrix> matp( matrix );
    double wmax = 0.98; // arbitrary, as in AimsTextureSmoothing
    std::map<unsigned, std::set< std::pair<unsigned,float> > >  weightLapl= AimsMeshWeightFiniteElementLaplacian( mesh.begin()->second, wmax );

    unsigned line, nl = matp.nlines(), col, nc = matp.ncols();
    double dur = sigma * sigma / ( 8. * ::log(2) );
    cout << "dur: " << dur << endl;
    float dt = 0.001; // FIXME
    unsigned it, niter = rint(dur/dt);
    vector<double> outrow;
    cout << "laplacian smoothing... niter: " << niter << ", dt: " << dt
      << endl;
    cout << "lines: " << nl << ", cols: " << nc << endl;
    double tmp;

#define use_matpow
// use_matpow: use coef matrix ^^ niter, calculated once
#ifdef use_matpow
    LaplacianWeights* laplmat
      = makeLaplacianSmoothingCoefficients( weightLapl, niter, dt, wthresh );
#endif

    double tminglob = FLT_MAX, tmaxglob = -FLT_MAX,
      otmin=FLT_MAX, otmax=-FLT_MAX;
    typename MatrixProxy<Matrix>::iterator iline = matp.begin();

    for( line=0; line<nl; ++line, ++iline )
    {
      cout << "\rsmoothing line " << line << " / " << nl << "..." << flush;
      vector<double> row( nc, 0. );
      typename MatrixProxy<Matrix>::columniterator ic = matp.colbegin( iline ),
        ec = matp.colend( iline );
      double tmin = FLT_MAX, tmax = -FLT_MAX;

      for( ; ic!=ec; ++ic )
      {
        double value = matp.value( ic );
        row[ matp.column( ic ) ] = value;
        if( value < tmin )
        {
          tmin = value;
          if( tmin < tminglob )
            tminglob = tmin;
        }
        if( value > tmax )
        {
          tmax = value;
          if( tmax > tmaxglob )
            tmaxglob = tmax;
        }
      }
      if( tmin > 0. )
        tmin = 0.;
      if( tmax < 0. )
        tmax = 0.;
      if( tminglob > 0. )
        tminglob = 0.;
      if( tmaxglob < 0. )
        tmaxglob = 0.;

#ifdef use_matpow
      applyLaplacianMatrix( row, outrow, *laplmat );
#else

//       otex =  AimsMeshLaplacian(itex,weightLapl);
//       cout << "Estime dt\n";
//       dt = AimsMeshFiniteElementDt(itex,otex,H);
//       cout << "dt estimated : " << dt << endl;

      for( it=0; it<niter; ++it )
      {
        AimsMeshLaplacian( row, outrow, weightLapl );
        for( col=0; col<nc; ++col)
          row[col] += dt * outrow[col];
      }
#endif

      for( col=0; col<nc; ++col )
      {
#ifdef use_matpow
        tmp = outrow[col];
#else
        tmp = row[col];
#endif
        if( tmp != 0 ) // could use thresholding to ensure sparsity
        {
          if( tmp < otmin )
            otmin = tmp;
          if( tmp > otmax )
            otmax = tmp;
          // clamp values between initial bounds
          if( tmp < tmin )
            tmp = tmin;
          else if( tmp > tmax )
            tmp = tmax;
          matp.assign( line, col, tmp, iline );
        }
      }
    }

    cout << endl;
    cout << "input matrix bounds: " << tminglob << " / " << tmaxglob << endl;
    cout << "output matrix bounds: " << otmin << " / " << otmax << endl;
  }

}


namespace constel
{

  void sparseMatrixDiffusionSmoothing( SparseMatrix & matrix,
    const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma )
  {
    _sparseMatrixDiffusionSmoothing( matrix, mesh, wthresh, sigma );
  }


  void sparseMatrixDiffusionSmoothing( Connectivities * conn_ptr,
    const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma )
  {
    _sparseMatrixDiffusionSmoothing( *conn_ptr, mesh, wthresh, sigma );
  }


} // namespace constel
