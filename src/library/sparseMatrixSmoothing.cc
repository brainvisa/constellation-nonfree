#include <constellation/sparseMatrixSmoothing.h>
#include <constellation/tildefs.h>
#include <aims/mesh/curv.h>
#include <cathier/triangle_mesh_geodesic_map.h>
#include <boost/numeric/ublas/matrix.hpp> // zero_matrix is needed in operation
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <float.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;
using namespace boost;


namespace constel
{
  // currently in sparseMatrixTools
  aims::SparseMatrix* connectivitiesToSparseMatrix(
    const Connectivities & conn );
}


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
    SparseMatrix* toSparseMatrix();

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


  template <>
  inline
  SparseMatrix* MatrixProxy<SparseMatrix>::toSparseMatrix()
  {
    return &_matrix;
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


  template <>
  inline
  SparseMatrix* MatrixProxy<Connectivities>::toSparseMatrix()
  {
    return connectivitiesToSparseMatrix( _matrix );
  }

  //

  size_t countElements( const SparseMatrix & matrix )
  {
    const boost_sparse_matrix & bmat = matrix.boostMatrix();
    size_t count = 0;
    boost_sparse_matrix::const_iterator1 il, el = bmat.end1();
    boost_sparse_matrix::const_iterator2 ic, ec;
    for( il=bmat.begin1(); il!=el; ++il )
    {
      for( ic=il.begin(), ec=il.end(); ic!=ec; ++ic )
      {
        ++count;
      }
    }
    cout << "countElements in mat " << bmat.size1() << " x " << bmat.size2() << ": " << count << endl;
    return count;
  }


  SparseMatrix* subMatrix(
    const LaplacianWeights & matrix, const vector<size_t> & items )
  {
    SparseMatrix * res = new SparseMatrix( items.size(), items.size() );
    boost_sparse_matrix & rmat = res->boostMatrix();
    LaplacianWeights::const_iterator il, el = matrix.end();
    size_t i, j, n = items.size();
    size_t line;
    size_t count = 0;
    for( i=0; i<n; ++i )
    {
      line = items[i];
      il = matrix.find( line );
      if( il != el )
      {
        set<pair<unsigned,float> >::const_iterator
          ic = il->second.begin(), ec = il->second.end();
        size_t col;
        for( j=0; j!=n; ++j )
        {
          col = items[j];
          while( ic != ec && ic->first < col )
            ++ic;
          if( ic != ec && ic->first == col )
          {
            rmat.insert_element( i, j, ic->second );
            ++count;
          }
        }
      }
    }
    cout << "subMatrix items: " << count << endl;
    return res;
  }


  SparseMatrix* toSparseMatrix( LaplacianWeights & matrix )
  {
    SparseMatrix * res = new SparseMatrix( matrix.size(), matrix.size() );
    boost_sparse_matrix & rmat = res->boostMatrix();
    LaplacianWeights::const_iterator il, el = matrix.end();
    size_t line;
    size_t count = 0;
    for( il=matrix.begin(); il!=el; ++il )
    {
      line = il->first;
      set<pair<unsigned,float> >::const_iterator
        ic, ec = il->second.end();
      for( ic=il->second.begin(); ic!=ec; ++ic )
      {
        rmat( line, ic->first ) = ic->second;
        ++count;
      }
    }
    cout << "SparseMatrix items: " << count << endl;
    return res;
  }


  SparseMatrix* trans( const SparseMatrix & matrix )
  {
    SparseMatrix * res = new SparseMatrix( matrix.getSize2(),
                                           matrix.getSize1() );
    const boost_sparse_matrix & bmat = matrix.boostMatrix();
    boost_sparse_matrix & omat = res->boostMatrix();
    boost_sparse_matrix::const_iterator1 il, el = bmat.end1();
    boost_sparse_matrix::const_iterator2 ic, ec;
    for( il=bmat.begin1(); il!=el; ++il )
    {
      for( ic=il.begin(), ec=il.end(); ic!=ec; ++ic )
        omat( ic.index2(), ic.index1() ) = *ic;
    }
    return res;
  }


  //

  template <typename Matrix>
  void _sparseMatrixDiffusionSmoothing( Matrix & matrix,
    const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
    const vector<size_t> & indices )
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

    cout << "Laplacian coefs done.\n";
    SparseMatrix *submat = subMatrix( *laplmat, indices );
    countElements( *submat );
    cout << "Submatrix done. Getting SparseMatrix\n";
    SparseMatrix *smat = matp.toSparseMatrix();
    cout << "OK, applying for columns\n";
    countElements( *smat );

    SparseMatrix ax( nl, nc );
    sparse_prod( submat->boostMatrix(), smat->boostMatrix(), ax.boostMatrix() ); // smooth along cols
    cout << "ax:\n";
    countElements( ax );
    delete submat;
    cout << "getting SparseMatrix for coefs\n";
    SparseMatrix* lmat = toSparseMatrix( *laplmat );
    delete laplmat;
    cout << "OK, transposing coefs\n";
    SparseMatrix *tlapl = trans( *lmat ); // transpose A
    delete lmat;
    cout << "tlapl:\n";
    countElements( *tlapl );
    cout << "OK, applying for lines\n";
    smat->boostMatrix().clear();
    sparse_prod( ax.boostMatrix(), tlapl->boostMatrix(), smat->boostMatrix()  ); // smooth along lines
    cout << "done.\n";
    countElements( *smat );
    delete tlapl;
    ax.boostMatrix().clear();

    // delete smat;
    return;

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
      delete laplmat;
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
    const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
    const vector<size_t> & indices )
  {
    _sparseMatrixDiffusionSmoothing( matrix, mesh, wthresh, sigma, indices );
  }


  void sparseMatrixDiffusionSmoothing( Connectivities * conn_ptr,
    const AimsTimeSurface<3,Void> & mesh, double wthresh, double sigma,
    const vector<size_t> & indices )
  {
    _sparseMatrixDiffusionSmoothing( *conn_ptr, mesh, wthresh, sigma,
                                     indices );
  }


  void sparseMatrixGaussianSmoothing(SparseMatrix & matrix,
    const AimsSurfaceTriangle & inAimsMesh, float distthresh, float wthresh )
  {
    /*
    Smoothing of a connectivity matrix according to the aims mesh and to the neighbourhood distance (distthresh), threshold of the resulting matrix by wthresh
    */
    //Convert inAimsMesh to Mesh mesh (Pascal Cathier Format)
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);

// Generating kdtree
    std::cout << "Generating kdtree" << std::endl;
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    std::cout << "getVertices(mesh):" << getVertices(mesh)[0] << ", " << getVertices(mesh)[1]  << std::endl;

    // Comuting geomap : neighborhood map
    //if distthresh!= 0
    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity
    double square_sigma = distthresh*distthresh;
    std::cout << "Computing geomap..." << std::flush;
    std::vector<QuickMap> res(getVertices(mesh).size());
    til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
    shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
    til::Triangle_mesh_geodesic_map<Mesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
        geomap(getVertices(mesh), *pneighc, stopGhost);
    std::vector<std::size_t> startPoints(1);
    std::vector<double> dist(1, 0.0);
    std::vector<std::size_t> nneigh(til::size(getVertices(mesh)));
    {
      for (std::size_t i = 0; i < til::size(getVertices(mesh)); ++i)
      {
        startPoints[0] = i;
        geomap.init(startPoints, dist);
        geomap.process();
        shared_ptr<til::sparse_vector<double> > tmp = geomap.distanceMap();
        res[i].resize(tmp->getMap().size());
        {
          using namespace til::expr;
          til::detail::loop_xx(castTo(*_1, *_2), res[i], tmp->getMap());
        }

        nneigh[i] = res[i].size();
      }
    }
    std::cout << "OK" << std::endl;
    
    //Sparse matrix smoothing
    SparseMatrix * smoothed_matrix_ptr = new SparseMatrix(matrix.getSize1(),matrix.getSize2());
    SparseMatrix & smoothed_matrix = *smoothed_matrix_ptr;
    SparseMatrix::iterator1 s1;
    SparseMatrix::iterator2 s2;
    if (getVertices(mesh).size()==matrix.getSize1() and getVertices(mesh).size()==matrix.getSize2())
    {
      for ( s1 = matrix.begin1(); s1 != matrix.end1(); s1++ )
      {
        std::size_t A = s1.index1();
        for ( s2 = s1.begin(); s2 != s1.end(); s2++ )
        {
          std::size_t B = s2.index2();
          double w = 0;
          for (QuickMap::const_iterator B_neigh = res[B].begin(); B_neigh != res[B].end(); ++B_neigh)
          {
            for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
            {
              double e1 = til::square(A_neigh->second);
              double e2 = til::square(B_neigh->second);
              w = (*s2)*std::exp( -(e1 + e2)  / ( 2*distthresh*distthresh ))/(square_sigma * two_pi);//"smoothing coefficient"
              if (w > G_THRESH)
              {
                smoothed_matrix(A_neigh->first,B_neigh->first)+=w;
              }
            }
          }
        }
      }
    }
    else if (getVertices(mesh).size()==matrix.getSize2())
    {
      for ( s2 = matrix.begin2(); s2 != matrix.end2(); s2++ )
      {
        std::size_t A = s2.index2();
        
        for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
        {
          double e1 = til::square(A_neigh->second);
          
          double smooth_coef = std::exp( -e1  / ( 2*distthresh*distthresh ));
//           /(square_sigma * two_pi);//"smoothing coefficient"
//           if (e1==0)
//           {
//             std::cout << "zero dist:" << e1 << ", weight:"<< smooth_coef <<std::endl;
//           }
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

  }// sparseMatrixGaussianSmoothing


  void sparseMatrixGaussianSmoothingNormed(SparseMatrix & matrix, const AimsSurfaceTriangle & inAimsMesh, float distthresh, float wthresh)
  {
    /*
    Smoothing of a connectivity matrix according to the aims mesh and to the neighbourhood distance (distthresh), threshold of the resulting matrix by wthresh
    Each connection contribution is normalized to 1
    */
    //Convert inAimsMesh to Mesh mesh (Pascal Cathier Format)
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);

// Generating kdtree
    std::cout << "Generating kdtree" << std::endl;
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    std::cout << "getVertices(mesh):" << getVertices(mesh)[0] << ", " << getVertices(mesh)[1]  << std::endl;

    // Comuting geomap : neighborhood map
    //if distthresh!= 0
    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity
    double square_sigma = distthresh*distthresh;
    std::cout << "Computing geomap..." << std::flush;
    std::vector<QuickMap> res(getVertices(mesh).size());
    til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
    shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
    til::Triangle_mesh_geodesic_map<Mesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
        geomap(getVertices(mesh), *pneighc, stopGhost);
    std::vector<std::size_t> startPoints(1);
    std::vector<double> dist(1, 0.0);
    std::vector<std::size_t> nneigh(til::size(getVertices(mesh)));
    {
      for (std::size_t i = 0; i < til::size(getVertices(mesh)); ++i)
      {
        startPoints[0] = i;
        geomap.init(startPoints, dist);
        geomap.process();
        shared_ptr<til::sparse_vector<double> > tmp = geomap.distanceMap();
        res[i].resize(tmp->getMap().size());
        {
          using namespace til::expr;
          til::detail::loop_xx(castTo(*_1, *_2), res[i], tmp->getMap());
        }

        nneigh[i] = res[i].size();
      }
    }
    std::cout << "OK" << std::endl;
    
    //Sparse matrix smoothing
    SparseMatrix * smoothed_matrix_ptr = new SparseMatrix(matrix.getSize1(),matrix.getSize2());
    SparseMatrix & smoothed_matrix = *smoothed_matrix_ptr;
    SparseMatrix::iterator1 s1;
    SparseMatrix::iterator2 s2;
    if (getVertices(mesh).size()==matrix.getSize1() and getVertices(mesh).size()==matrix.getSize2())
    {
      for ( s1 = matrix.begin1(); s1 != matrix.end1(); s1++ )
      {
        std::size_t A = s1.index1();
        for ( s2 = s1.begin(); s2 != s1.end(); s2++ )
        {
          std::size_t B = s2.index2();
          double w = 0;
          for (QuickMap::const_iterator B_neigh = res[B].begin(); B_neigh != res[B].end(); ++B_neigh)
          {
            for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
            {
              double e1 = til::square(A_neigh->second);
              double e2 = til::square(B_neigh->second);
              w = (*s2)*std::exp( -(e1 + e2)  / ( 2*distthresh*distthresh ))/(square_sigma * two_pi);//"smoothing coefficient"
              if (w > G_THRESH)
              {
                smoothed_matrix(A_neigh->first,B_neigh->first)+=w;
              }
            }
          }
        }
      }
    }
    else if (getVertices(mesh).size()==matrix.getSize2())
    {
      for ( s2 = matrix.begin2(); s2 != matrix.end2(); s2++ )
      {
        std::size_t A = s2.index2();
        float a = res[A].size();
        std::vector<float> weightsVector(a);
        for (std::size_t i = 0; i <a; i++)
        {
          weightsVector[i]=0;
        }
//         std::cout << weightsVector[0] << std::endl;
        std::size_t neighACount = 0;
        for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
        {
          double e1 = til::square(A_neigh->second);
          
          double smooth_coef = std::exp( -e1  / ( 2*distthresh*distthresh ));
//           /(square_sigma * two_pi);//"smoothing coefficient"
//           if (e1==0)
//           {
//             std::cout << "zero dist:" << e1 << ", weight:"<< smooth_coef <<std::endl;
//           }
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
        for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
        {
          double e1 = til::square(A_neigh->second);
          
          double smooth_coef = std::exp( -e1  / ( 2*distthresh*distthresh ));
//           /(square_sigma * two_pi);//"smoothing coefficient"
//           if (e1==0)
//           {
//             std::cout << "zero dist:" << e1 << ", weight:"<< smooth_coef <<std::endl;
//           }
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
