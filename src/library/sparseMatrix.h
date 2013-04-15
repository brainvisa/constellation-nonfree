#ifndef AIMS_SPARSEMATRIX_H
#define AIMS_SPARSEMATRIX_H

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/version.hpp>
#include <string>

#include "cathier/SparseVector.h"
namespace aims
{


  class SparseCMatrix;
#if BOOST_VERSION >= 103300
  typedef boost::numeric::ublas::mapped_matrix< double > boost_sparse_matrix;
#else
  typedef boost::numeric::ublas::sparse_matrix< double > boost_sparse_matrix;
#endif


  class SparseMatrix
  {

    public:

      typedef boost_sparse_matrix::reference
          reference;
      typedef boost_sparse_matrix::const_reference
          const_reference;

      typedef boost_sparse_matrix::iterator1 iterator1;
      typedef boost_sparse_matrix::iterator2 iterator2;

      typedef boost_sparse_matrix::const_iterator1
          const_iterator1;
      typedef boost_sparse_matrix::const_iterator2
          const_iterator2;

      SparseMatrix( int32_t size1 = 1, int32_t size2 = 1 );
      SparseMatrix( const SparseMatrix& other );
//       SparseMatrix( const Matrix& other );
      virtual ~SparseMatrix();

      SparseMatrix& operator=( const SparseMatrix& other );
//       SparseMatrix& operator=( const Matrix& other );

      void reallocate( int32_t size1, int32_t size2 );

      iterator1 begin1();
      iterator1 end1();

      const_iterator1 begin1() const;
      const_iterator1 end1() const;

      iterator2 begin2();
      iterator2 end2();

      const_iterator2 begin2() const;
      const_iterator2 end2() const;

      int32_t getSize1() const;
      int32_t getSize2() const;
      int32_t getNonZeroElementCount() const;

      bool hasElement( int32_t i, int32_t j ) const;

      const_reference operator()( int32_t i, int32_t j ) const; 
      reference operator()( int32_t i, int32_t j ); 

      bool isSquare() const;

      void setZero();
      void setIdentity();//bizarre ?
      void setRow( int32_t s1, const std::vector<double>& row );
      void setColumn( int32_t s2, const std::vector<double>& column );
      void fill( const double& value );
      void fill( int32_t offset1,
                 int32_t offset2,
                 int32_t size1,
                 int32_t size2,
                 const double& value );
      void fill( int32_t offset1,
                 int32_t offset2,
                 const SparseMatrix& other,
                 int32_t size1 = 0,
                 int32_t size2 = 0 );
//       void fill( int32_t offset1,
//                  int32_t offset2,
//                  const Matrix& other,
//                  int32_t size1 = 0,
//                  int32_t size2 = 0 );
      void setDiagonal( const double& value );
      void transpose();

      double getTrace() const;
      SparseMatrix getTransposition() const;
      SparseMatrix getComposition( const SparseMatrix& other ) const;
      std::vector<double> getComposition( const std::vector<double>& other ) const;
//       SparseCMatrix getComposition( const SparseCMatrix& other ) const;
//       Cstd::vector getComposition( const Cstd::vector& other ) const;
      std::vector<double> toVector() const;
//       Cstd::vector toVector() const;
//       Matrix toMatrix() const;

      std::vector<double> getRow( int32_t i ) const;
      std::vector<double> getColumn( int32_t j ) const;
      til::SparseVector<double> getSparseRow( int32_t i ) const;
      til::SparseVector<double> getSparseColumn( int32_t j ) const;
      void erase_element( int32_t i, int32_t j );

      void read( const std::string& filename,
                 const std::string & openmode = "binar", 
                 bool bswap = false );
      void write( const std::string& filename,
                  const std::string & openmode = "binar", 
                  bool bswap = false ) const;
      aims::SparseMatrix &
        operator += ( const aims::SparseMatrix& thing );
      aims::SparseMatrix &
        operator -= ( const aims::SparseMatrix& thing );
      aims::SparseMatrix &
        operator *= ( double x );
      aims::SparseMatrix &
        operator /= ( double x );

    protected:

      boost_sparse_matrix _matrix;

  };


}


//
// unary operators + & -
//

// + SparseMatrix
aims::SparseMatrix 
    operator + ( const aims::SparseMatrix& thing );

// - SparseMatrix
aims::SparseMatrix 
    operator - ( const aims::SparseMatrix& thing );

//
// binary operator +
//

// SparseMatrix + SparseMatrix
aims::SparseMatrix
    operator + ( const aims::SparseMatrix& thing1,
                 const aims::SparseMatrix& thing2 );

// SparseMatrix + SparseCMatrix
// aims::SparseCMatrix
//     operator + ( const aims::SparseMatrix& thing1,
//                  const aims::SparseCMatrix& thing2 );


//
// binary operator -
//

// SparseMatrix - SparseMatrix
aims::SparseMatrix
    operator - ( const aims::SparseMatrix& thing1,
                 const aims::SparseMatrix& thing2 );

// SparseMatrix - SparseCMatrix
// aims::SparseCMatrix
//     operator - ( const aims::SparseMatrix& thing1,
//                  const aims::SparseCMatrix& thing2 );


//
// binary operator *
//

// SparseMatrix * SparseMatrix
aims::SparseMatrix
    operator * ( const aims::SparseMatrix& thing1,
                 const aims::SparseMatrix& thing2 );

// SparseMatrix * SparseCMatrix
// aims::SparseCMatrix
//     operator * ( const aims::SparseMatrix& thing1,
//                  const aims::SparseCMatrix& thing2 );

// SparseMatrix * std::vector
std::vector<double>
    operator * ( const aims::SparseMatrix& thing1,
                 const std::vector<double>& thing2 );

// SparseMatrix * Cstd::vector
// aims::Cstd::vector
//     operator * ( const aims::SparseMatrix& thing1,
//                  const aims::Cstd::vector& thing2 );

// SparseMatrix * double
aims::SparseMatrix
    operator * ( const aims::SparseMatrix& thing1,
                 const double& thing2 );

// SparseMatrix * complex< double >
// aims::SparseCMatrix
//     operator * ( const aims::SparseMatrix& thing1,
//                  const std::complex< double >& thing2 );

//
// binary operator /
//

// SparseMatrix / double
aims::SparseMatrix
    operator / ( const aims::SparseMatrix& thing1,
                 const double& thing2 );

// SparseMatrix / complex< double >
// aims::SparseCMatrix
//     operator / ( const aims::SparseMatrix& thing1,
//                  const std::complex< double >& thing2 );


namespace std
{

  std::ostream& operator << ( ostream& os, const aims::SparseMatrix& thing );
  bool operator == ( const aims::SparseMatrix& thing1,
                     const aims::SparseMatrix& thing2 );

}


#endif
