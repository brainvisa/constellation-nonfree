#include <constellation/connMatrixTools.h>
#include <constellation/textureAndMeshTools.h>
#include <aims/mesh/curv.h>
#include <til/numeric_array.h>
#include <til/numeric_array_tools.h>
#include <float.h>

using namespace comist;
using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;

namespace constel
{

  void fillconnMatrix(Connectivities * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, double distanceThreshold, uint connectionLength)
  {
    double two_pi = 2*3.1415926535897931;
    /*
    inputs:
              conn_ptr: connectivity matrix to fill, symmetric in this case
              fiberExtremity1NeighMeshVertex: std::vector< std::pair<std::size_t, double> > of pairs of:  < i1 = mesh vertex index near the first fiber extremity, distance between i1 and first fiber extremity closest vertex in the mesh (or intersection between fiber and mesh)>
              fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
              connectivityThreshold: threshold: connectivity contributions below are excluded
              distanceThreshold: distance for the smoothing of the connectivity matrix

    output:
              void, (*conn_ptr has changed: filled according to input mesh neighborhoods (for fiber extrimities))
    */
    double wtotal = 0;
    Connectivities & conn = *conn_ptr;
    list<pair<pair<size_t, size_t>, double> > weights;
    double square_sigma = distanceThreshold*distanceThreshold;
    double sm_coef = 1. / (square_sigma * two_pi); // should be 1./(sigma*sqrt(two_pi)) ?
    double exp_coef = 1. / ( 2*square_sigma );
    for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex.begin(); neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2)
    {
      for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex.begin(); neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1)
      {
        double e1 = til::square(neigh_1->second);
        double e2 = til::square(neigh_2->second);
        // double w = std::exp( -(e1 + e2)  / ( 2*square_sigma ))/(square_sigma * two_pi);//"smoothing coefficient"
        double w1 = std::exp( -e1 * exp_coef ) * sm_coef;//"smoothing coefficient"
        double w2 = std::exp( -e2 * exp_coef ) * sm_coef;
        if( w1 < connectivityThreshold )
          w1 = 0;
        if( w2 < connectivityThreshold )
          w2 = 0;
        double w = w1 + w2;
        if( w > 0 )
        {
          wtotal += w;
          weights.push_back( make_pair( make_pair( neigh_1->first, neigh_2->first ), w ) );
        }
      }
    }
//         std::cout << w << std::endl;
//         std::cout <<"e1 :" << e1 << ", e2:" << e2 << ", w:" << w << std::endl;
        //if (w > connectivityThreshold)
    list<pair<pair<size_t, size_t>, double> >::iterator iw, ew = weights.end();
    for ( iw=weights.begin(); iw!=ew; ++iw )
    {
      double w = connectionLength * iw->second / wtotal;
      conn[iw->first.second][iw->first.first] += w;
      conn[iw->first.first][iw->first.second] += w;
//           std::cout <<"w > connectivityThreshold:" << w << std::endl;
//           std::cout << "(i,j): (" << neigh_1->first << "," << neigh_2->first << ")" << std::endl;
    }
  }//fillconnMatrix


  void fillconnMatrixNoSmoothing(Connectivities * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex)
  {
    /*
    inputs:
              conn_ptr: connectivity matrix to fill, symmetric in this case
              fiberExtremity1NeighMeshVertex: std::vector< std::pair<std::size_t, double> > of pairs of:  < i1 = mesh vertex index near the first fiber extremity, distance between i1 and first fiber extremity closest vertex in the mesh (or intersection between fiber and mesh)>
              fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
              connectivityThreshold: threshold: connectivity contributions below are excluded

    output:
              void, (*conn_ptr has changed: filled according to input mesh neighborhoods (for fiber extrimities))
    */
    
    Connectivities & conn = *conn_ptr;
    double isz = 1. / ( fiberExtremity1NeighMeshVertex.size()
      * fiberExtremity2NeighMeshVertex.size() );
    for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex.begin(); neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2)
    {
      for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex.begin(); neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1)
      {
        conn[neigh_1->first][neigh_2->first] += isz;
        conn[neigh_2->first][neigh_1->first] += isz;
//           std::cout <<"w > connectivityThreshold:" << w << std::endl;
//           std::cout << "(i,j): (" << neigh_1->first << "," << neigh_2->first << ")" << std::endl;
      }
    }
  }//fillconnMatrixNoSmoothing


  void fillconnMatrix(aims::SparseMatrix * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, double distanceThreshold, uint connectionLength)
    {
      double two_pi = 2*3.1415926535897931;
      /*
      inputs:
                conn_ptr: connectivity matrix to fill, symmetric in this case
                fiberExtremity1NeighMeshVertex: std::vector< std::pair<std::size_t, double> > of pairs of:  < i1 = mesh vertex index near the first fiber extremity, distance between i1 and first fiber extremity closest vertex in the mesh (or intersection between fiber and mesh)>
                fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
                connectivityThreshold: threshold: connectivity contributions below are excluded
                distanceThreshold: distance for the smoothing of the connectivity matrix
  
      output:
                void, (*conn_ptr has changed: filled according to input mesh neighborhoods (for fiber extrimities))
      */
      
      aims::SparseMatrix & conn = *conn_ptr;
      for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex.begin(); neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2)
      {
        for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex.begin(); neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1)
        {
          double e1 = til::square(neigh_1->second);
          double e2 = til::square(neigh_2->second);
          double square_sigma = distanceThreshold*distanceThreshold;
          double w = std::exp( -(e1 + e2)  / ( 2*square_sigma ))/(square_sigma*two_pi);//"smoothing coefficient"
//           std::cout <<"e1 :" << e1 << ", e2:" << e2 << ", w:" << w << std::endl;
          if (w > connectivityThreshold)
          {
            w = connectionLength *w;
            conn(neigh_1->first,neigh_2->first) = conn(neigh_1->first,neigh_2->first) + w;
            conn(neigh_2->first,neigh_1->first) = conn(neigh_2->first,neigh_1->first) + w; 
//             std::cout <<"w > connectivityThreshold:" << w << std::endl;
  //           std::cout << "(i,j): (" << neigh_1->first << "," << neigh_2->first << ")" << std::endl;
          }
        }
      }
    }//fillconnMatrixSparse

  void fillconnMatrix(aims::SparseMatrix * conn_ptr, aims::SparseMatrix * conn_ptr2, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, double distanceThreshold, std::size_t rowIndex_min, std::size_t rowIndex_max, uint connectionLength)
    {
      double two_pi = 2*3.1415926535897931;
      /*
      inputs:
                conn_ptr: connectivity matrix to fill, symmetric in this case
                fiberExtremity1NeighMeshVertex: std::vector< std::pair<std::size_t, double> > of pairs of:  < i1 = mesh vertex index near the first fiber extremity, distance between i1 and first fiber extremity closest vertex in the mesh (or intersection between fiber and mesh)>
                fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
                connectivityThreshold: threshold: connectivity contributions below are excluded
                distanceThreshold: distance for the smoothing of the connectivity matrix
  
      output:
                void, (*conn_ptr has changed: filled according to input mesh neighborhoods (for fiber extrimities))
      */
      
      aims::SparseMatrix & conn = *conn_ptr;
      aims::SparseMatrix & conn2 = *conn_ptr2;
      for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex.begin(); neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2)
      {
        for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex.begin(); neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1)
        {
          double e1 = til::square(neigh_1->second);
          double e2 = til::square(neigh_2->second);
          double square_sigma = distanceThreshold*distanceThreshold;
          double w = std::exp( -(e1 + e2)  / ( 2*square_sigma ))/(square_sigma*two_pi);//"smoothing coefficient"
//           std::cout <<"e1 :" << e1 << ", e2:" << e2 << ", w:" << w << std::endl;
          if (w > connectivityThreshold)
          {
            w = connectionLength *w;
            if (rowIndex_max == 0 or (not conn_ptr2))
            {

              conn(neigh_1->first,neigh_2->first) = conn(neigh_1->first,neigh_2->first) + w;
              conn(neigh_2->first,neigh_1->first) = conn(neigh_2->first,neigh_1->first) + w;

            }
            else
            {
              if (neigh_1->first >= rowIndex_min and neigh_1->first < rowIndex_max)
              {
                conn(neigh_1->first,neigh_2->first) = conn(neigh_1->first,neigh_2->first) + w;
              }
              else
              {
                conn2(neigh_1->first,neigh_2->first) = conn2(neigh_1->first,neigh_2->first) + w;
              }
              if (neigh_2->first >= rowIndex_min and neigh_2->first < rowIndex_max)
              {
                conn(neigh_2->first,neigh_1->first) = conn(neigh_2->first,neigh_1->first) + w;
              }
              else
              {
                conn2(neigh_2->first,neigh_1->first) = conn2(neigh_2->first,neigh_1->first) + w;
              }
            }
            
//             std::cout <<"w > connectivityThreshold:" << w << std::endl;
  //           std::cout << "(i,j): (" << neigh_1->first << "," << neigh_2->first << ")" << std::endl;
          }
        }
      }
    }//fillconnTwoMatrixSparse

//   void fillconnMatrixWithoutSmoothing(Connectivities * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, uint connectionLength)
//   {
//     /*
//     inputs:
//               conn_ptr: connectivity matrix to fill, symmetric in this case
//               fiberExtremity1NeighMeshVertex: std::vector< std::pair<std::size_t, double> > of pairs of:  < i1 = mesh vertex index near the first fiber extremity, distance between i1 and first fiber extremity closest vertex in the mesh (or intersection between fiber and mesh)>
//               fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
//               connectivityThreshold: threshold: connectivity contributions below are excluded
//               distanceThreshold: distance for the smoothing of the connectivity matrix
// 
//     output:
//               void, (*conn_ptr has changed: filled according to input mesh neighborhoods (for fiber extrimities))
//     */
//     double meanDistNeigh1 = 0.0;
//     for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex.begin(); neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1)
//     {
//       meanDistNeigh1 += neigh_1->second;
//     }
//     double meanDistNeigh2 = 0.0;
//     for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex.begin(); neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2)
//     {
//       meanDistNeigh2 += neigh_2->second;
//     }
//     std::cout << "meanDistNeigh1:" << meanDistNeigh1 << ", meanDistNeigh2:" << meanDistNeigh2 << std::endl;
//     Connectivities & conn = *conn_ptr;
//     for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex.begin(); neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2)
//     {
//       for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex.begin(); neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1)
//       {
//         double e1 = (neigh_1->second)/meanDistNeigh1;
//         double e2 = (neigh_2->second);
//         double w = std::exp( -(e1 + e2)  / ( 2*distanceThreshold*distanceThreshold ));//"smoothing coefficient"
//         std::cout <<"e1 :" << e1 << ", e2:" << e2 << ", w:" << w << std::endl;
//         if (w > connectivityThreshold)
//         {
//           w = connectionLength *w;
//           conn[neigh_1->first][neigh_2->first] += w;
//           conn[neigh_2->first][neigh_1->first] += w;
//           std::cout <<"w > connectivityThreshold:" << w << std::endl;
// //           std::cout << "(i,j): (" << neigh_1->first << "," << neigh_2->first << ")" << std::endl;
//         }
//       }
//     }
//   }//fillconnMatrixWithoutSmoothing
  void fillconnMatrixWeightingLength(Connectivities * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, double distanceThreshold, float connectionLength)
  {
    /*
    inputs:
              conn_ptr: connectivity matrix to fill, symmetric in this case
              fiberExtremity1NeighMeshVertex: std::vector< std::pair<std::size_t, double> > of pairs of:  < i1 = mesh vertex index near the first fiber extremity, distance between i1 and first fiber extremity closest vertex in the mesh (or intersection between fiber and mesh)>
              fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
              connectivityThreshold: threshold: connectivity contributions below are excluded
              distanceThreshold: distance for the smoothing of the connectivity matrix

    output:
              void, (*conn_ptr has changed: filled according to input mesh neighborhoods (for fiber extrimities))
    */
    
    Connectivities & conn = *conn_ptr;
    for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex.begin(); neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2)
    {
      for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex.begin(); neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1)
      {
        double e1 = til::square(neigh_1->second);
        double e2 = til::square(neigh_2->second);
        double w = std::exp( -(e1 + e2)  / ( 2*distanceThreshold*distanceThreshold ));//"smoothing coefficient"
        if (w > connectivityThreshold)
        {
          float connectionLenghtWeight = connectionLength;
          w = connectionLenghtWeight *w;
          conn[neigh_1->first][neigh_2->first] += w;
          conn[neigh_2->first][neigh_1->first] += w;
//           std::cout << "(i,j): (" << neigh_1->first << "," << neigh_2->first << ")" << std::endl;
        }
      }
    }
  }//fillconnMatrixWeightinfLength

  void fillNonSymetricConnMatrix(Connectivities * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex_rows, QuickMap & fiberExtremity2NeighMeshVertex_cols, double connectivityThreshold, double distanceThreshold)
  {
    /*
    inputs:
    conn_ptr: connectivity matrix to fill, non symmetric in this case : shape = (cortexMeshVertex_nb, thalamusMeshVertex_nb) for example
    fiberExtremity1NeighMeshVertex: std::vector< std::pair<std::size_t, double> > of pairs of:  < i1 = mesh vertex index near the first fiber extremity, distance between i1 and first fiber extremity closest vertex in the mesh (or intersection between fiber and mesh)>
    fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
    connectivityThreshold: threshold: connectivity contributions below are excluded
    distanceThreshold: distance for the smoothing of the connectivity matrix

    output:
    void, (*conn_ptr has changed: filled according to input mesh neighborhoods (for fiber extrimities))
    */
    Connectivities & conn = *conn_ptr;
    for (QuickMap::const_iterator neigh_2 = fiberExtremity2NeighMeshVertex_cols.begin(); neigh_2 != fiberExtremity2NeighMeshVertex_cols.end(); ++neigh_2)
    {
      for(QuickMap::const_iterator neigh_1 = fiberExtremity1NeighMeshVertex_rows.begin(); neigh_1 != fiberExtremity1NeighMeshVertex_rows.end(); ++neigh_1)
      {
        double e1 = til::square(neigh_1->second);
        double e2 = til::square(neigh_2->second);
        double w = std::exp( -(e1 + e2)  / ( 2*distanceThreshold*distanceThreshold ));//"smoothing coefficient"
        if (w > connectivityThreshold)
        {
          conn[neigh_1->first][neigh_2->first] += w;
        }
      }
    }
  }//fillNonSymetricConnMatrix

  template<int D, class T>
      bool computeIntersectionPointFiberSegmentAndMesh(const AimsTimeSurface<D,T> & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, QuickMap ** polygonVerticesDistMap)
  {
    //fill polygonVerticesDistMap with the D vertex of the intersection polygon (if there is an intersection)
    const std::vector<Point3df>  & meshVertices  = aimsMesh.vertex();
    const std::vector< AimsVector<uint,D> > & meshPolygons = aimsMesh.polygon();
    std::set<uint> closest_polygons = polygonsByVertex_Index[meshClosestPoint];
    std::set<uint>::const_iterator set_it;
    AimsVector<uint,D> currentAims_polygon;
    std::vector< Point3df > current_polygon(D);
    Point3df * intersection_point = 0;
    bool inter = false;
    for (set_it = closest_polygons.begin(); set_it != closest_polygons.end(); set_it++)
    {
      currentAims_polygon = meshPolygons[*set_it];
      for ( uint i = 0; i < D; ++i)
      {
        current_polygon[i] = meshVertices[currentAims_polygon[i]];
      }
      inter = computeIntersectionSegmentPolygon( fiberPoint1, fiberPoint2, current_polygon, & intersection_point );
      if (inter)
      {
        uint currentPolygonVerticeIndex;
        Point3df currentPolygonVertice;
        * polygonVerticesDistMap = new QuickMap(D);
        til::numeric_array<float, 3> intersection_point_na((*intersection_point)[0], (*intersection_point)[1], (*intersection_point)[2]);
        for  ( uint j = 0; j < D; ++j)
        {
          currentPolygonVerticeIndex = currentAims_polygon[j];
          currentPolygonVertice = meshVertices[currentPolygonVerticeIndex];
          til::numeric_array<float, 3> currentPolygonVertice_na(currentPolygonVertice[0], currentPolygonVertice[1], currentPolygonVertice[2]);
          std::pair<std::size_t, double > pair_verticeIndex_distanceToIntersection = std::make_pair(currentAims_polygon[j],sqrt(dist2(currentPolygonVertice_na, intersection_point_na)));
          (**polygonVerticesDistMap)[j] = pair_verticeIndex_distanceToIntersection;
        }
        break;
      }
    }
    delete intersection_point;
    return inter;
  }//computeIntersectionPointFiberSegmentAndMesh
  
  template<int D, class T>
      bool computeIntersectionPointFiberSegmentAndMesh2(const AimsTimeSurface<D,T> & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, QuickMap ** closestPointsTointersectionWeightMap)
  {
    //fill closestPointsTointersectionWeightMap with the polygon vertice closest to the intersection point, and the weight value = 1.0
    const std::vector<Point3df>  & meshVertices  = aimsMesh.vertex();
    const std::vector< AimsVector<uint,D> > & meshPolygons = aimsMesh.polygon();
    std::set<uint> closest_polygons = polygonsByVertex_Index[meshClosestPoint];
    std::set<uint>::const_iterator set_it;
    AimsVector<uint,D> currentAims_polygon;
    std::vector< Point3df > current_polygon(D);
    Point3df * intersection_point = 0;
    bool inter = false;
    for (set_it = closest_polygons.begin(); set_it != closest_polygons.end(); set_it++)
    {
      currentAims_polygon = meshPolygons[*set_it];
      for ( uint i = 0; i < D; ++i)
      {
        current_polygon[i] = meshVertices[currentAims_polygon[i]];
      }
      inter = computeIntersectionSegmentPolygon( fiberPoint1, fiberPoint2, current_polygon, & intersection_point );
      if (inter)
      {
        uint currentPolygonVerticeIndex;
        Point3df currentPolygonVertice;
        * closestPointsTointersectionWeightMap = new QuickMap(1);
        til::numeric_array<float, 3> intersection_point_na((*intersection_point)[0], (*intersection_point)[1], (*intersection_point)[2]);
        
        std::vector <float> polygonVerticesIntersection_distances(D);
        float polygonVerticesIntersection_minDistance = 1000;//mm
        uint closestVerticeNb_On_D;
        std::size_t closestVerticeIndex_On_mesh;
        for  ( uint j = 0; j < D; ++j)
        {
          currentPolygonVerticeIndex = currentAims_polygon[j];
          currentPolygonVertice = meshVertices[currentPolygonVerticeIndex];
          til::numeric_array<float, 3> currentPolygonVertice_na(currentPolygonVertice[0], currentPolygonVertice[1], currentPolygonVertice[2]);
          polygonVerticesIntersection_distances[j] = dist2(currentPolygonVertice_na, intersection_point_na);
          if ( polygonVerticesIntersection_distances[j] < polygonVerticesIntersection_minDistance)
          {
            polygonVerticesIntersection_minDistance = polygonVerticesIntersection_distances[j];
            closestVerticeNb_On_D = j;
            closestVerticeIndex_On_mesh = currentPolygonVerticeIndex;
          }
        }
        
        std::pair<std::size_t, double > pair_verticeIndex_distanceToIntersection = std::make_pair(closestVerticeIndex_On_mesh,1.0);
        (**closestPointsTointersectionWeightMap)[0] = pair_verticeIndex_distanceToIntersection;
        break;
      }
    }
    delete intersection_point;
    return inter;
  }//computeIntersectionPointFiberSegmentAndMesh2
  
  template<int D, class T>
      bool computeIntersectionPointNeighborhoodFiberSegmentAndMesh(const AimsTimeSurface<D,T> & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, std::vector<QuickMap> & distanceThresholdNeighborhoodByVertex, QuickMap ** polygonVerticesDistMap_2ptr)
  {
    //fill polygonVerticesDistMap with the vertex closest to the intersection point and its neighborhood according to distanceThreshold value, if there is an intersection (a changer: fastmarching a partir des sommets du polygone mieux pour construire le voisinage)
    const std::vector<Point3df>  & meshVertices  = aimsMesh.vertex();
    const std::vector< AimsVector<uint,D> > & meshPolygons = aimsMesh.polygon();
    std::set<uint> closest_polygons = polygonsByVertex_Index[meshClosestPoint];
    std::set<uint>::const_iterator set_it;
    AimsVector<uint,D> currentAims_polygon;
    std::vector< Point3df > current_polygon(D);
    Point3df * intersection_point = 0;
    bool inter = false;
    for (set_it = closest_polygons.begin(); set_it != closest_polygons.end(); set_it++)
    {
      currentAims_polygon = meshPolygons[*set_it];
      for ( uint i = 0; i < D; ++i)
      {
        current_polygon[i] = meshVertices[currentAims_polygon[i]];
      }
      inter = computeIntersectionSegmentPolygon( fiberPoint1, fiberPoint2, current_polygon, & intersection_point );
      if (inter)
      {
        uint currentPolygonVerticeIndex;
        Point3df currentPolygonVertice;
        * polygonVerticesDistMap_2ptr = new QuickMap(D);
        til::numeric_array<float, 3> intersection_point_na((*intersection_point)[0], (*intersection_point)[1], (*intersection_point)[2]);
        
        constel::QuickMap & polygonVerticesDistMap = **polygonVerticesDistMap_2ptr;
        std::vector <float> polygonVerticesIntersection_distances(D);
        float polygonVerticesIntersection_minDistance = 1000;//mm
        uint closestVerticeNb_On_D;
        for  ( uint j = 0; j < D; ++j)
        {
          currentPolygonVerticeIndex = currentAims_polygon[j];
          currentPolygonVertice = meshVertices[currentPolygonVerticeIndex];
          til::numeric_array<float, 3> currentPolygonVertice_na(currentPolygonVertice[0], currentPolygonVertice[1], currentPolygonVertice[2]);
          polygonVerticesIntersection_distances[j] = dist2(currentPolygonVertice_na, intersection_point_na);
          if ( polygonVerticesIntersection_distances[j] < polygonVerticesIntersection_minDistance)
          {
            polygonVerticesIntersection_minDistance = polygonVerticesIntersection_distances[j];
            closestVerticeNb_On_D = j;
          }
//           polygonVerticesDistMap = distanceThresholdNeighborhoodByVertex[currentAims_polygon[j]];
//           std::cout << "closestVerticeNb_On_D:" << closestVerticeNb_On_D << std::endl;
//           std::cout << "j:" << j << std::endl;
//           std::pair<std::size_t, double > pair_verticeIndex_distanceToIntersection = std::make_pair(currentAims_polygon[j],dist2(currentPolygonVertice_na, intersection_point_na));
//           (**polygonVerticesDistMap_2ptr)[j] = pair_verticeIndex_distanceToIntersection;
        }
        polygonVerticesDistMap = distanceThresholdNeighborhoodByVertex[currentAims_polygon[closestVerticeNb_On_D]];
        break;
      }
    }
    delete intersection_point;
    return inter;
  }//computeIntersectionPointNeighborhoodFiberSegmentAndMesh

  void fillconnMatrixWithConnections(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold)
  {
    size_t connectionsCount = 0;
    if (distanceThreshold>0.0)
    {
      for (BundleConnections::const_iterator iConnection = connections.begin(); iConnection != connections.end(); ++iConnection, ++connectionsCount)
      {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr, ConnectionMapFront, ConnectionMapBack, connectivityThreshold, distanceThreshold);
  //       std::cout << "connection !" << std::endl;
      }
    }
   else
   {
      for (BundleConnections::const_iterator iConnection = connections.begin(); iConnection != connections.end(); ++iConnection, ++connectionsCount)
      {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrixNoSmoothing(conn_ptr, ConnectionMapFront, ConnectionMapBack);
  //       std::cout << "connection !" << std::endl;
      }
   }
    std::cout << "connectionsCount:" << connectionsCount << std::endl;
  }//fillconnMatrixWithConnections
  

  void fillconnMatrixWithConnections(aims::SparseMatrix * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, std::size_t rowIndex_min, std::size_t rowIndex_max,  aims::SparseMatrix * conn_ptr2)
  {
    size_t connectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin(); iConnection != connections.end(); ++iConnection, ++connectionsCount)
    {
      constel::QuickMap ConnectionMapFront = iConnection->front();
      constel::QuickMap ConnectionMapBack = iConnection->back();
      if (not conn_ptr2)
      {
        fillconnMatrix(conn_ptr, ConnectionMapFront, ConnectionMapBack, connectivityThreshold, distanceThreshold);
      }
      else
      {
        fillconnMatrix(conn_ptr, conn_ptr2, ConnectionMapFront, ConnectionMapBack, connectivityThreshold, distanceThreshold, rowIndex_min, rowIndex_max);
      }
//       std::cout << "connection !" << std::endl;
    }
    std::cout << "connectionsCount:" << connectionsCount << std::endl;
  }//fillconnMatrixWithConnections



  void fillconnMatrixWithConnectionsPlusLength(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, uint length_min, uint length_max, ConnectionsLength & connectionsLength)
  {
    size_t connectionsCount = 0;
    size_t chosenConnectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin(); iConnection != connections.end(); ++iConnection, ++connectionsCount)
    {
      uint connection_length = connectionsLength[connectionsCount];
      if ( connection_length < length_max and connection_length >= length_min )
      {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr, ConnectionMapFront, ConnectionMapBack, connectivityThreshold, distanceThreshold);
        ++chosenConnectionsCount;
      }
    }
    std::cout << "connectionsCount:" << connectionsCount << std::endl;
    std::cout << "chosen connectionsCount:" << chosenConnectionsCount << std::endl;
  }//fillconnMatrixWithConnectionsPlusLength

  void fillconnMatrixWithConnectionsPlusLengthWeight(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, uint length_min, uint length_max, ConnectionsLength & connectionsLength)
  {
    size_t connectionsCount = 0;
    size_t chosenConnectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin(); iConnection != connections.end(); ++iConnection, ++connectionsCount)
    {
      uint connection_length = connectionsLength[connectionsCount];
      if ( connection_length < length_max and connection_length >= length_min )
      {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr, ConnectionMapFront, ConnectionMapBack, connectivityThreshold, distanceThreshold, connection_length);
        ++chosenConnectionsCount;
      }
    }
    std::cout << "connectionsCount:" << connectionsCount << std::endl;
    std::cout << "chosen connectionsCount:" << chosenConnectionsCount << std::endl;
  }//fillconnMatrixWithConnectionsPlusLengthWeight
  
  void fillconnMatrixWithConnectionsPlusFloatLengthWeight(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, float length_min, float length_max, ConnectionsFloatLength & connectionsLength)
  {
    size_t connectionsCount = 0;
    size_t chosenConnectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin(); iConnection != connections.end(); ++iConnection, ++connectionsCount)
    {
      float connection_length = connectionsLength[connectionsCount];
      if ( connection_length < length_max and connection_length >= length_min )
      {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr, ConnectionMapFront, ConnectionMapBack, connectivityThreshold, distanceThreshold, connection_length);
        ++chosenConnectionsCount;
      }
    }
    std::cout << "connectionsCount:" << connectionsCount << std::endl;
    std::cout << "chosen connectionsCount:" << chosenConnectionsCount << std::endl;
  }//fillconnMatrixWithConnectionsPlusFloatLengthWeight
  
  void fillNonSymetricConnMatrixWithConnections(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold)
  {
    size_t connectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin(); iConnection != connections.end(); ++iConnection, ++connectionsCount)
    {
      constel::QuickMap ConnectionMapFront = iConnection->front();
      constel::QuickMap ConnectionMapBack = iConnection->back();
      fillNonSymetricConnMatrix(conn_ptr, ConnectionMapFront, ConnectionMapBack, connectivityThreshold, distanceThreshold);
    }
    std::cout << "connectionsCount:" << connectionsCount << std::endl;
  }//fillNonSymetricConnMatrixWithConnections
  template bool
      computeIntersectionPointFiberSegmentAndMesh(const AimsSurfaceTriangle & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, QuickMap ** polygonVerticesDistMap);
  
  template bool
      computeIntersectionPointFiberSegmentAndMesh2(const AimsSurfaceTriangle & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, QuickMap ** polygonVerticesDistMap);

  template bool
      computeIntersectionPointNeighborhoodFiberSegmentAndMesh(const AimsSurfaceTriangle & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, std::vector<QuickMap> & distanceThresholdNeighborhoodByVertex, QuickMap ** polygonVerticesDistMap_2ptr);


  aims::SparseMatrix* connectivitiesToSparseMatrix(
    const Connectivities & conn )
  {
    size_t size1 = conn.size();
    size_t size2 = conn[0].size();

    SparseMatrix *mp = new SparseMatrix(size1,size2);
    SparseMatrix & m = *mp;
    Connectivities::const_iterator il, el = conn.end();
    size_t i, j;
    Connectivity::sparse_const_iterator ic, ec;
    for( il=conn.begin(), i=0; il!=el; ++il, ++i )
    {
      for( ic=il->sparse_begin(), ec=il->sparse_end(); ic!=ec; ++ic )
      {
        j = ic->first;
        double val = ic->second;
        if( val!= 0)
        {
          m(i,j) = val;
        }
      }
    }

    return mp;
  }


  Connectivities* sparseMatrixToConnectivities(
    const aims::SparseMatrix & mat )
  {
    Connectivities* conn = new Connectivities( mat.getSize1(),
      Connectivity( mat.getSize2() ) );
    SparseMatrix::const_iterator1 il, el = mat.end1();
    SparseMatrix::const_iterator2 ic, ec;
    unsigned i;
    for( il=mat.begin1(); il!=el; ++il )
    {
      i = il.index1();
      for( ic=il.begin(), ec=il.end(); ic!=ec; ++ic )
        (*conn)[ i ][ ic.index2() ] = *ic;
    }

    return conn;
  }


  void writeConnectivities( const Connectivities & conn,
                            const string & filename, bool ascii )
  {
    aims::SparseMatrix *m = connectivitiesToSparseMatrix( conn );
    try
    {
      if( ascii )
        m->write( filename, "ascii" );
      else
        m->write( filename );
    }
    catch( ... )
    {
      delete m;
      throw;
    }
    delete m;
  }


} // namespace constel
