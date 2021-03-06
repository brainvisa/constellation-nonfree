#include <constellation/connMatrixTools.h>
#include <constellation/tildefs.h>
#include <constellation/connectivities.h>
#include <aims/io/writer.h>
#include <aims/distancemap/meshdistance.h>

using namespace aims;
using namespace carto;
using namespace std;

namespace constel
{

  typedef Volume<float> Matrix;

  //------------------------------
  //  connMatrixTargetsToTargets
  //------------------------------
  void connMatrixTargetsToTargets(
      const Fibers &fibers, const AimsSurfaceTriangle &inAimsMesh,
      Motion motion, const TimeTexture<short> &targetRegionsTex,
      string file_name)
  {
    /*  Computing connectivity matrix regrouped according to targetRegions
        given in targetRegionsTexture (labeled texture, between 0 and
        targetRegionsNb, 0 = background, 1 to targetRegionsNb = target regions
    */
    const Texture<short> &targetTex0 = targetRegionsTex.begin()->second;
    const vector<short> &targetTex0Data = targetTex0.data();
    const short maxTargetTex0Data
      = *(max_element(targetTex0Data.begin(),targetTex0Data.end()));
    int targetRegionsNb = int(maxTargetTex0Data);
    cout <<"targetRegionsNb:" << targetRegionsNb << endl;

    //Connectectivity matrix initialization:
    size_t colNb = targetRegionsNb;
    size_t rowsNb = targetRegionsNb;
    Matrix matrix(rowsNb, colNb, 1);
    matrix.setVoxelSize(100.0/rowsNb , 80.0/colNb);

    // Generating kdtree
    KDTreeVertices m = kdt_vertices( inAimsMesh );
    KDTree kdt( m.begin(), m.end() );

    int countRemoved = 0;
    size_t fiberCount = 0, nFibers = fibers.size();
    size_t five_count = nFibers / 20.0;
    Point3df p1, p2;

    for (Fibers::const_iterator iFiber = fibers.begin();
         iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (five_count != 0) {
        if (fiberCount % five_count == 0)
        {
          cout << 5 * int(fiberCount / five_count) << "%..." << flush;
        }
      }
      else
      {
        cout << "...100%..." << flush;
      }
      //Looking for closest points
      // Transform iFiber->front() and iFiber->back() from t2 to anat space
      // (with motion)
      p1 = motion.transform(iFiber->front()[0],
                            iFiber->front()[1],
                            iFiber->front()[2]);
      p2 = motion.transform(iFiber->back()[0],
                            iFiber->back()[1],
                            iFiber->back()[2]);
      size_t A = kdt.find_nearest( make_pair( 0L, p1 ) ).first->first;
      size_t B = kdt.find_nearest( make_pair( 0L, p2 ) ).first->first;

      if( dist2( p1, inAimsMesh.vertex()[A] ) > 25.0 ||
          dist2( p2, inAimsMesh.vertex()[B] ) > 25.0)
      {
        ++countRemoved;
        continue;
      }
      // Filling the connectivity matrix
      int labelA = targetTex0[A];
      int labelB = targetTex0[B];
      bool validLabelA = (1 <= labelA && labelA <= targetRegionsNb);
      bool validLabelB = (1 <= labelB && labelB <= targetRegionsNb);
      if (validLabelA && validLabelB)
      {
        matrix(labelA-1,labelB-1,0) += 1;
        matrix(labelB-1,labelA-1,0) += 1;
      }
    }
    Writer<Matrix> w(file_name);
    w.write(matrix);
  }

  //------------------------
  //  connMatrixSeedRegion
  //------------------------
  void connMatrixSeedRegion(
      const Fibers &fibers, const AimsSurfaceTriangle &inAimsMesh,
      Motion motion, const TimeTexture<short> &seedRegionsTex,
      size_t seedRegionLabel, string connmatrix_filename,
      string /* connTexture_filename */ )
  {
    /*  Computing connectivity matrix of a given seed Region
        defined by  seedRegionLabel vertices in seedRegionTex
        (labeled texture, 0 = background)
        without smoothing
    */

    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    size_t meshVertexNb = seedTex0.nItem();
    //Identify seedRegion:
    size_t seedRegionVertexNb = 0;
    vector<size_t > seedVertexIndexInMatrix(meshVertexNb);
    size_t label= 0;
    for (size_t i = 0; i < meshVertexNb; ++i)
    {
      label = seedTex0.item(i);
      if (label == seedRegionLabel)
      {
        seedVertexIndexInMatrix[i]=seedRegionVertexNb;
        seedRegionVertexNb++;
      }
    }

    //Connectivity matrix initialization:
    size_t colNb = meshVertexNb;
    size_t rowsNb = seedRegionVertexNb;
    Matrix matrix(rowsNb, colNb, 1);
    matrix.setVoxelSize(100.0/rowsNb , 80.0/colNb);

    // Generating kdtree
    cout << "Generating kdtree" << endl;
    KDTreeVertices m = kdt_vertices( inAimsMesh );
    KDTree kdt( m.begin(), m.end() );

    //Connectivity matrix filling:
    int countRemoved = 0;
    size_t fiberCount = 0, nFibers = fibers.size();
    size_t five_count = nFibers / 20.0;
    Point3df p1, p2;

    for (Fibers::const_iterator iFiber = fibers.begin();
         iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          cout << 5 * int(fiberCount / five_count) << "%..." << flush;
        }
      }
      else
      {
        cout << "...100%..." << flush;
      }
      //Looking for closest points
      // Transform iFiber->front() and iFiber->back() from t2 to anat space
      // (with motion)
      p1 = motion.transform(iFiber->front()[0],
                            iFiber->front()[1],
                            iFiber->front()[2]);
      p2 = motion.transform(iFiber->back()[0],
                            iFiber->back()[1],
                            iFiber->back()[2]);
      size_t A = kdt.find_nearest( make_pair( 0U, p1 ) ).first->first;
      size_t B = kdt.find_nearest( make_pair( 0U, p2 ) ).first->first;

      if( dist2( p1, inAimsMesh.vertex()[A] ) > 25.0 ||
          dist2( p2, inAimsMesh.vertex()[B] ) > 25.0)
      {
        ++countRemoved;
        continue;
      }
      // Filling the connectivity matrix
      size_t labelA = seedTex0[A];
      size_t labelB = seedTex0[B];
      bool validLabelA = labelA == seedRegionLabel;
      bool validLabelB = labelB == seedRegionLabel;
      if (validLabelA)
      {
        matrix(seedVertexIndexInMatrix[A],B,0) += 1;
      }
      if (validLabelB)
      {
        matrix(seedVertexIndexInMatrix[B],A,0) += 1;
      }

    }
    Writer<Matrix> w(connmatrix_filename);
    w.write(matrix);
  }

  //--------------------------------
  //  connMatrixSeedRegionSmoothed
  //--------------------------------
  void connMatrixSeedRegionSmoothed(
      const Fibers &fibers, const AimsSurfaceTriangle &inAimsMesh,
      Motion motion, const TimeTexture<short> &seedRegionsTex,
      size_t seedRegionLabel, float distthresh, float wthresh,
      string connmatrix_filename, string /* connTexture_filename */,
      bool logOption)
  {
    /*  Computing connectivity matrix of a given seed Region
        defined by  seedRegionLabel vertices in seedRegionTex
        (labeled texture, 0=background) with smoothing according to distthresh
    */

    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    size_t meshVertexNb = seedTex0.nItem();
    //Identify seedRegion:
    size_t seedRegionVertexNb = 0;
    vector<size_t > seedVertexIndexInMatrix(meshVertexNb);
    size_t label = 0;
    for (size_t i = 0; i < meshVertexNb; ++i)
    {
      label = seedTex0.item(i);
      if (label == seedRegionLabel)
      {
        seedVertexIndexInMatrix[i]=seedRegionVertexNb;
        seedRegionVertexNb++;
      }
    }

    //Connectivity matrix initialization:
    int32_t size1 = seedRegionVertexNb;
    int32_t size2 = meshVertexNb;
    SparseMatrix matrix(size1,size2);

    // Generating kdtree
    cout << "Generating kdtree" << endl;
    KDTreeVertices m = kdt_vertices( inAimsMesh );
    KDTree kdt( m.begin(), m.end() );

    // Comuting geomap : neighborhood map
    //if distthresh!= 0
    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity
    double square_sigma = distthresh*distthresh;
    cout << "Computing geomap..." << flush;

    vector<map<size_t, float> > res;

    meshdistance::pairwiseDistanceMaps( inAimsMesh.begin()->second, res,
                                        distthresh );

    //Connectivity matrix filling:
    int countRemoved = 0;
    size_t fiberCount = 0, nFibers = fibers.size();
    size_t five_count = nFibers / 20.0;
    Point3df p1, p2;

    for (Fibers::const_iterator iFiber = fibers.begin();
         iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          cout << 5 * int(fiberCount / five_count) << "%..." << flush;
        }
      }
      else
      {
        cout << "...100%..." << flush;
      }
      //Looking for closest points
      // Transform iFiber->front() and iFiber->back() from t2 to anat space
      // (with motion)
      p1 = motion.transform(iFiber->front()[0],
                            iFiber->front()[1],
                            iFiber->front()[2]);
      p2 = motion.transform(iFiber->back()[0],
                            iFiber->back()[1],
                            iFiber->back()[2]);
      size_t A = kdt.find_nearest( make_pair( 0U, p1 ) ).first->first;
      size_t B = kdt.find_nearest( make_pair( 0U, p2 ) ).first->first;

      if( dist2( p1, inAimsMesh.vertex()[A] ) > 25.0 ||
          dist2( p2, inAimsMesh.vertex()[B] ) > 25.0)
      {
        ++countRemoved;
        continue;
      }
      // Filling the connectivity matrix
      for (map<size_t, float>::const_iterator B_neigh = res[B].begin();
           B_neigh != res[B].end(); ++B_neigh)
      {
        for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
             A_neigh != res[A].end(); ++A_neigh)
        {
          size_t labelA = seedTex0[A_neigh->first];
          size_t labelB = seedTex0[B_neigh->first];
          bool validLabelA = labelA == seedRegionLabel;
          bool validLabelB = labelB == seedRegionLabel;
          double w = 0;
          if (validLabelA || validLabelB)
          {
            double e1 = square(A_neigh->second);
            double e2 = square(B_neigh->second);
            w = exp( -(e1 + e2)  / ( 2*distthresh*distthresh ))
                /(square_sigma * two_pi);//"smoothing coefficient"
            if (w > G_THRESH)
            {
              if (validLabelA)
              {
                matrix(seedVertexIndexInMatrix[A_neigh->first],B_neigh->first)
                  += w;
              }
              if (validLabelB)
              {
                matrix(seedVertexIndexInMatrix[B_neigh->first],A_neigh->first)
                  += w;
              }
            }
          }
        }
      }
    }
    if (logOption)
    {
      SparseMatrix::iterator1 s1;
      SparseMatrix::iterator2 s2;
      for (s1 = matrix.begin1(); s1 != matrix.end1(); s1++)
      {
        for (s2 = s1.begin(); s2 != s1.end(); s2++)
        {
          *s2 = log(1+*s2);
          if (*s2 <= wthresh)
          {
            matrix.erase_element(s1.index1(), s2.index2());
          }
        }
      }
    }
    matrix.write(connmatrix_filename);
  }

  //-------------------------------------------
  //  matrixSmoothingWithVertexAreaCorrection
  //-------------------------------------------
  void matrixSmoothingWithVertexAreaCorrection(
      SparseMatrix * matrix_ptr, const AimsSurfaceTriangle &inAimsMesh,
      float distthresh, float /* wthresh */)
  {
    /*
    Smoothing of a connectivity matrix according to the aims mesh and to the
    neighbourhood distance (distthresh), threshold of the resulting matrix by
    wthresh if the sampling of the mesh is not uniform:
    for each tract (non zero entry of the matrix), for each of its extremities
    (A and B): A for example: 
        - compute the area of its neighbourhood
        - for each neighbour vertex A_neigh_v of A:
          (neighboor = below a certain distance from A)
          - compute the distance from A
          - the "area" surrounding this vertex: sum of the areas of the
            neighbourhood polygons of A_neigh_v: A_neigh_v_area
        - Compute the total area of the neigbourhood: A_neigh_total_area
    */

    // Comuting geomap : neighborhood map
    //if distthresh!= 0
    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity
    double square_sigma = distthresh*distthresh;
    cout << "Computing geomap..." << flush;

    size_t nv = inAimsMesh.vertex().size();
    vector<map<size_t, float> > res;

    meshdistance::pairwiseDistanceMaps( inAimsMesh.begin()->second, res,
                                        distthresh );

    //Sparse matrix smoothing
    SparseMatrix & matrix = *matrix_ptr;
    SparseMatrix * smoothed_matrix_ptr
      = new SparseMatrix(matrix.getSize1(),matrix.getSize2());
    SparseMatrix & smoothed_matrix = *smoothed_matrix_ptr;
    SparseMatrix::iterator1 s1;
    SparseMatrix::iterator2 s2;
    if (( (int32_t) nv == matrix.getSize1() )
        && ((int32_t) nv ==matrix.getSize2()))
    {
      for (s1 = matrix.begin1(); s1 != matrix.end1(); s1++)
      {
        size_t A = s1.index1();
        for (s2 = s1.begin(); s2 != s1.end(); s2++)
        {
          size_t B = s2.index2();
          double w = 0;
          for (map<size_t, float>::const_iterator B_neigh = res[B].begin();
               B_neigh != res[B].end(); ++B_neigh)
          {
            for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
                 A_neigh != res[A].end(); ++A_neigh)
            {
              double e1 = square(A_neigh->second);
              double e2 = square(B_neigh->second);
              w = (*s2)*exp( -(e1 + e2)  / ( 2*distthresh*distthresh )) 
                  /(square_sigma * two_pi);//"smoothing coefficient"
              if (w > G_THRESH)
              {
                smoothed_matrix(A_neigh->first,B_neigh->first)+=w;
              }
            }
          }
        }
      }
    }
    else if ((int32_t) nv ==matrix.getSize2())
    {
      for ( s2 = matrix.begin2(); s2 != matrix.end2(); s2++ )
      {
        size_t A = s2.index2(); 
        for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
             A_neigh != res[A].end(); ++A_neigh)
        {
          double e1 = square(A_neigh->second);
          double smooth_coef = exp( -e1  / ( 2*distthresh*distthresh ));
          for (s1 = s2.begin(); s1 != s2.end(); s1++)
          {
            size_t B = s1.index1();
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
  }
} // namespace constel
