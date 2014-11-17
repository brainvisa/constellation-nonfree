#include <constellation/connMatrixTools.h>
#include <constellation/tildefs.h>
#include <constellation/connectivities.h>
// includes from CATHIER
#include <cathier/triangle_mesh_geodesic_map.h>


using namespace aims;
using namespace carto;
using namespace std;

namespace constel {

  typedef AimsData< float > Matrix;

  void connMatrixTargetsToTargets(const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh, Motion motion, const TimeTexture<short> & targetRegionsTex, std::string file_name)
  {
    /*  Computing connectivity matrix regrouped according to targetRegions given in targetRegionsTexture (labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to targetRegionsNb = target regions
    */

    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    const std::vector<short> & targetTex0Data = targetTex0.data();
    const short maxTargetTex0Data = *(std::max_element(targetTex0Data.begin(),targetTex0Data.end()));
    int targetRegionsNb = int(maxTargetTex0Data);
    std::cout <<"targetRegionsNb:" << targetRegionsNb << std::endl;
    

    //Connectectivity matrix initialization:
    std::size_t colNb = targetRegionsNb;
    std::size_t rowsNb = targetRegionsNb;
    Matrix matrix(rowsNb, colNb, 1);
    matrix.setSizeXYZT(100.0/rowsNb , 80.0/colNb);

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

    int countRemoved = 0;
    std::size_t fiberCount = 0, nFibers = fibers.size();
    std::size_t five_count = nFibers / 20.0;
    Point3df p1, p2;

    for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          std::cout << 5 * int(fiberCount / five_count) << "%..." << std::flush;
        }
      }
      else /*if(five_count == 0)*/
      {
        std::cout << "...100%..." << std::flush;
      }
      //Looking for closest points
      til::Find_closest< double, KDTree > fc(kdt);
      // Transform iFiber->front() and iFiber->back() from t2 to anat space (with motion)
      p1 = motion.transform(iFiber->front()[0],iFiber->front()[1],iFiber->front()[2]);
      til::numeric_array<float, 3> p1na(p1[0], p1[1], p1[2]);
      p2 = motion.transform(iFiber->back()[0],iFiber->back()[1],iFiber->back()[2]);
      til::numeric_array<float, 3> p2na(p2[0], p2[1], p2[2]);
      std::size_t A = fc(p1na);
//       std::cout << "p1na:" << p1na << "," << "p2na:" << p2na << std::endl;
//       std::cout << "mesh_closest_point:" << A << std::endl;
      std::size_t B = fc(p2na);

      if (til::dist2(p1na, getVertices(mesh)[A], til::prec<float>()) > 25.0 ||
          til::dist2(p2na, getVertices(mesh)[B], til::prec<float>()) > 25.0)
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
    Writer<Matrix> w( file_name );
    w.write( matrix );
  
  } // connMatrixTargetsToTargets

  void connMatrixSeedRegion(const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh, Motion motion, const TimeTexture<short> & seedRegionsTex, std::size_t seedRegionLabel, std::string connmatrix_filename, std::string connTexture_filename)
  {
    /*  Computing connectivity matrix of a given seed Region
        defined by  seedRegionLabel vertices in seedRegionTex
        (labeled texture, 0 = background)
        without smoothing
    */

    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    std::size_t meshVertexNb = seedTex0.nItem();
    //Identify seedRegion:
    std::size_t seedRegionVertexNb = 0;
    std::vector<std::size_t > seedVertexIndexInMatrix(meshVertexNb);
    std::size_t label= 0;
    for (std::size_t i = 0; i < meshVertexNb; ++i)
    {
      label = seedTex0.item(i);
      if (label == seedRegionLabel)
      {
        seedVertexIndexInMatrix[i]=seedRegionVertexNb;
        seedRegionVertexNb++;
      }
    }

    //Connectivity matrix initialization:
    std::size_t colNb = meshVertexNb;
    std::size_t rowsNb = seedRegionVertexNb;
    Matrix matrix(rowsNb, colNb, 1);
    matrix.setSizeXYZT(100.0/rowsNb , 80.0/colNb);

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

    //Connectivity matrix filling:
    int countRemoved = 0;
    std::size_t fiberCount = 0, nFibers = fibers.size();
    std::size_t five_count = nFibers / 20.0;
    Point3df p1, p2;

    for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          std::cout << 5 * int(fiberCount / five_count) << "%..." << std::flush;
        }
      }
      else /*if(five_count == 0)*/
      {
        std::cout << "...100%..." << std::flush;
      }
      //Looking for closest points
      til::Find_closest< double, KDTree > fc(kdt);
      // Transform iFiber->front() and iFiber->back() from t2 to anat space (with motion)
      p1 = motion.transform(iFiber->front()[0],iFiber->front()[1],iFiber->front()[2]);
      til::numeric_array<float, 3> p1na(p1[0], p1[1], p1[2]);
      p2 = motion.transform(iFiber->back()[0],iFiber->back()[1],iFiber->back()[2]);
      til::numeric_array<float, 3> p2na(p2[0], p2[1], p2[2]);
      std::size_t A = fc(p1na);
//       std::cout << "p1na:" << p1na << "," << "p2na:" << p2na << std::endl;
//       std::cout << "mesh_closest_point:" << A << std::endl;
      std::size_t B = fc(p2na);

      if (til::dist2(p1na, getVertices(mesh)[A], til::prec<float>()) > 25.0 ||
          til::dist2(p2na, getVertices(mesh)[B], til::prec<float>()) > 25.0)
      {
        ++countRemoved;
        continue;
      }
      // Filling the connectivity matrix
      std::size_t labelA = seedTex0[A];
      std::size_t labelB = seedTex0[B];
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
    Writer<Matrix> w( connmatrix_filename );
    w.write( matrix );

  } // connMatrixSeedRegion

  void connMatrixSeedRegionSmoothed(const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh, Motion motion, const TimeTexture<short> & seedRegionsTex, std::size_t seedRegionLabel, float distthresh, float wthresh, std::string connmatrix_filename, std::string connTexture_filename, bool logOption)
  {
    /*  Computing connectivity matrix of a given seed Region
        defined by  seedRegionLabel vertices in seedRegionTex
        (labeled texture, 0 = background) with smoothing according to distthresh
    */

    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    std::size_t meshVertexNb = seedTex0.nItem();
    //Identify seedRegion:
    std::size_t seedRegionVertexNb = 0;
    std::vector<std::size_t > seedVertexIndexInMatrix(meshVertexNb);
    std::size_t label = 0;
    for (std::size_t i = 0; i < meshVertexNb; ++i)
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
    boost::shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
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
        boost::shared_ptr<til::sparse_vector<double> > tmp = geomap.distanceMap();
        res[i].resize(tmp->getMap().size());
        {
          using namespace til::expr;
          til::detail::loop_xx(castTo(*_1, *_2), res[i], tmp->getMap());
        }

        nneigh[i] = res[i].size();
      }
    }
    std::cout << "OK" << std::endl;

    //Connectivity matrix filling:
    int countRemoved = 0;
    std::size_t fiberCount = 0, nFibers = fibers.size();
    std::size_t five_count = nFibers / 20.0;
    Point3df p1, p2;

    for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          std::cout << 5 * int(fiberCount / five_count) << "%..." << std::flush;
        }
      }
      else /*if(five_count == 0)*/
      {
        std::cout << "...100%..." << std::flush;
      }
      //Looking for closest points
      til::Find_closest< double, KDTree > fc(kdt);
      // Transform iFiber->front() and iFiber->back() from t2 to anat space (with motion)
      p1 = motion.transform(iFiber->front()[0],iFiber->front()[1],iFiber->front()[2]);
      til::numeric_array<float, 3> p1na(p1[0], p1[1], p1[2]);
      p2 = motion.transform(iFiber->back()[0],iFiber->back()[1],iFiber->back()[2]);
      til::numeric_array<float, 3> p2na(p2[0], p2[1], p2[2]);
      std::size_t A = fc(p1na);
      std::size_t B = fc(p2na);

      if (til::dist2(p1na, getVertices(mesh)[A], til::prec<float>()) > 25.0 ||
          til::dist2(p2na, getVertices(mesh)[B], til::prec<float>()) > 25.0)
      {
        ++countRemoved;
        continue;
      }
      // Filling the connectivity matrix
      for (QuickMap::const_iterator B_neigh = res[B].begin(); B_neigh != res[B].end(); ++B_neigh)
      {
        for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
        {
          std::size_t labelA = seedTex0[A_neigh->first];
          std::size_t labelB = seedTex0[B_neigh->first];
          bool validLabelA = labelA == seedRegionLabel;
          bool validLabelB = labelB == seedRegionLabel;
          double w = 0;
          if (validLabelA || validLabelB)
          {
            double e1 = til::square(A_neigh->second);
            double e2 = til::square(B_neigh->second);
            w = std::exp( -(e1 + e2)  / ( 2*distthresh*distthresh ))/(square_sigma * two_pi);//"smoothing coefficient"
            if (w > G_THRESH)
            {
              if (validLabelA)
              {
                matrix(seedVertexIndexInMatrix[A_neigh->first],B_neigh->first) += w;
              }
              if (validLabelB)
              {
                matrix(seedVertexIndexInMatrix[B_neigh->first],A_neigh->first) += w;
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
      for ( s1 = matrix.begin1(); s1 != matrix.end1(); s1++ )
      {

        for ( s2 = s1.begin(); s2 != s1.end(); s2++ )
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
    

  } // connMatrixSeedRegionSmoothed


  void matrixSmoothingWithVertexAreaCorrection(SparseMatrix * matrix_ptr, const AimsSurfaceTriangle & inAimsMesh, float distthresh, float wthresh)
  {
    /*
    Smoothing of a connectivity matrix according to the aims mesh and to the neighbourhood distance (distthresh), threshold of the resulting matrix by wthresh
    if the sampling of the mesh is not uniform:
    for each tract (non zero entry of the matrix), for each of its extremities (A and B):
      A for example: 
        - compute the area of its neighbourhood
        - for each neighbour vertex A_neigh_v of A: (neighboor = below a certain distance from A)
          - compute the distance from A
          - the "area" surrounding this vertex: sum of the areas of the neighbourhood polygons of A_neigh_v: A_neigh_v_area
        - Compute the total area of the neigbourhood: A_neigh_total_area
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
    boost::shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
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
        boost::shared_ptr<til::sparse_vector<double> > tmp = geomap.distanceMap();
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
    SparseMatrix & matrix = *matrix_ptr;
    std::cout << "coucou1" << std::endl;
    SparseMatrix * smoothed_matrix_ptr = new SparseMatrix(matrix.getSize1(),matrix.getSize2());
    std::cout << "coucou1" << std::endl;
    SparseMatrix & smoothed_matrix = *smoothed_matrix_ptr;
    std::cout << "coucou1" << std::endl;
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
  
  }// matrixSmoothingWithVertexAreaCorrection

} // namespace constel
