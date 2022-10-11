#include <constellation/connMatrixTools.h>
#include <constellation/textureAndMeshTools.h>
#include <aims/mesh/curv.h>
#include <float.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;

namespace
{
  float dist2( const Point3df & p1, const Point3df & p2 )
  {
    return (p2 - p1).norm2();
  }
}

namespace constel {

  //------------------
  //  fillconnMatrix
  //------------------
  /*
  inputs:
      conn_ptr: connectivity matrix to fill, symmetric in this case
      fiberExtremity1NeighMeshVertex: vector< pair<size_t, double> >
      of pairs of:  < i1 = mesh vertex index near the first fiber extremity,
      distance between i1 and first fiber extremity closest vertex in the mesh
      (or intersection between fiber and mesh)>
      fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
      connectivityThreshold: threshold: connectivity contributions below are
      excluded
      distanceThreshold: distance for the smoothing of the connectivity matrix
   output:
      void (*conn_ptr has changed: filled according to input mesh neighborhoods
      (for fiber extrimities))
  */
  void fillconnMatrix(
      Connectivities *conn_ptr, QuickMap &fiberExtremity1NeighMeshVertex,
      QuickMap &fiberExtremity2NeighMeshVertex, double connectivityThreshold,
      double distanceThreshold, unsigned connectionLength)
  {
    double two_pi = 2*3.1415926535897931;
    double wtotal = 0;
    Connectivities &conn = *conn_ptr;
    list<pair<pair<size_t, size_t>, double> > weights;
    double square_sigma = distanceThreshold*distanceThreshold;
    // should be 1./(sigma*sqrt(two_pi)) ?
    double sm_coef = 1. / (square_sigma * two_pi);
    double exp_coef = 1. / ( 2*square_sigma );
    for (QuickMap::const_iterator neigh_2
           = fiberExtremity2NeighMeshVertex.begin();
         neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2) {
      for (QuickMap::const_iterator neigh_1
             = fiberExtremity1NeighMeshVertex.begin();
           neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1) {
        double e1 = square(neigh_1->second);
        double e2 = square(neigh_2->second);
        double w1 = exp( -e1 * exp_coef ) * sm_coef;
        double w2 = exp( -e2 * exp_coef ) * sm_coef;
        if (w1 < connectivityThreshold) w1 = 0;
        if (w2 < connectivityThreshold) w2 = 0;
        double w = w1 + w2;
        if (w > 0) {
          wtotal += w;
          weights.push_back(
              make_pair(make_pair( neigh_1->first, neigh_2->first), w));
        }
      }
    }
    list<pair<pair<size_t, size_t>, double> >::iterator iw, ew
      = weights.end();
    for (iw=weights.begin(); iw!=ew; ++iw) {
      double w = connectionLength * iw->second / wtotal;
      conn[iw->first.second][iw->first.first] += w;
      conn[iw->first.first][iw->first.second] += w;
    }
  }

  //-----------------------------
  //  fillconnMatrixNoSmoothing
  //-----------------------------
  /*
  inputs:
      conn_ptr: connectivity matrix to fill, symmetric in this case
      fiberExtremity1NeighMeshVertex: vector<pair<size_t, double> >
      of pairs of:  < i1 = mesh vertex index near the first fiber extremity,
      distance between i1 and first fiber extremity closest vertex in the mesh
      (or intersection between fiber and mesh)>
      fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
      connectivityThreshold: threshold: connectivity contributions below are
      excluded

  output:
      void (*conn_ptr has changed: filled according to input mesh neighborhoods
      (for fiber extrimities))
  */
  void fillconnMatrixNoSmoothing(
      Connectivities *conn_ptr, QuickMap &fiberExtremity1NeighMeshVertex,
      QuickMap & fiberExtremity2NeighMeshVertex) {
    Connectivities &conn = *conn_ptr;
    double isz
      = 1. / (fiberExtremity1NeighMeshVertex.size()
              * fiberExtremity2NeighMeshVertex.size());
    for (QuickMap::const_iterator neigh_2
           = fiberExtremity2NeighMeshVertex.begin();
         neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2) {
      for (QuickMap::const_iterator neigh_1
             = fiberExtremity1NeighMeshVertex.begin();
           neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1) {
        conn[neigh_1->first][neigh_2->first] += isz;
        conn[neigh_2->first][neigh_1->first] += isz;
      }
    }
  }

  void fillWeightedConnMatrixNoSmoothing(
      Connectivities *conn_ptr, QuickMap &fiberExtremity1NeighMeshVertex,
      QuickMap & fiberExtremity2NeighMeshVertex, double weight) {
    Connectivities &conn = *conn_ptr;
    double isz
      = weight / (fiberExtremity1NeighMeshVertex.size()
              * fiberExtremity2NeighMeshVertex.size());
    for (QuickMap::const_iterator neigh_2
           = fiberExtremity2NeighMeshVertex.begin();
         neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2) {
      for (QuickMap::const_iterator neigh_1
             = fiberExtremity1NeighMeshVertex.begin();
           neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1) {
        conn[neigh_1->first][neigh_2->first] += isz;
        conn[neigh_2->first][neigh_1->first] += isz;
      }
    }
  }

  //------------------
  //  fillconnMatrix
  //------------------
  /*
  inputs:
      conn_ptr: connectivity matrix to fill, symmetric in this case
      fiberExtremity1NeighMeshVertex: vector<pair<size_t, double> >
      of pairs of:  < i1 = mesh vertex index near the first fiber extremity,
      distance between i1 and first fiber extremity closest vertex in the mesh
      (or intersection between fiber and mesh)>
      fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
      connectivityThreshold: threshold: connectivity contributions below are
      excluded
      distanceThreshold: distance for the smoothing of the connectivity matrix

  output:
      void (*conn_ptr has changed: filled according to input mesh neighborhoods
      (for fiber extrimities))
  */
  void fillconnMatrix(
      aims::SparseMatrix *conn_ptr, QuickMap &fiberExtremity1NeighMeshVertex,
      QuickMap &fiberExtremity2NeighMeshVertex, double connectivityThreshold,
      double distanceThreshold, unsigned connectionLength) {
    double two_pi = 2*3.1415926535897931;
    aims::SparseMatrix & conn = *conn_ptr;
    for (QuickMap::const_iterator neigh_2
           = fiberExtremity2NeighMeshVertex.begin();
         neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2) {
      for (QuickMap::const_iterator neigh_1
             = fiberExtremity1NeighMeshVertex.begin();
           neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1) {
        double e1 = square(neigh_1->second);
        double e2 = square(neigh_2->second);
        double square_sigma = distanceThreshold*distanceThreshold;
        double w = exp( -(e1 + e2)
                   / ( 2*square_sigma ))/(square_sigma*two_pi);
        if (w > connectivityThreshold) {
          w = connectionLength *w;
          conn(neigh_1->first,neigh_2->first)
            = conn(neigh_1->first,neigh_2->first) + w;
          conn(neigh_2->first,neigh_1->first)
            = conn(neigh_2->first,neigh_1->first) + w;
        }
      }
    }
  }

  //------------------
  //  fillconnMatrix
  //------------------
  /*
  inputs:
      conn_ptr: connectivity matrix to fill, symmetric in this case
      fiberExtremity1NeighMeshVertex: vector<pair<size_t, double> >
      of pairs of:  < i1 = mesh vertex index near the first fiber extremity,
      distance between i1 and first fiber extremity closest vertex in the
      mesh (or intersection between fiber and mesh)>
      fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
      connectivityThreshold: threshold: connectivity contributions below are
      excluded
      distanceThreshold: distance for the smoothing of the connectivity
      matrix
  output:
      void (*conn_ptr has changed: filled according to input mesh
      neighborhoods (for fiber extrimities))
  */
  void fillconnMatrix(
      aims::SparseMatrix *conn_ptr, aims::SparseMatrix *conn_ptr2,
      QuickMap &fiberExtremity1NeighMeshVertex,
      QuickMap &fiberExtremity2NeighMeshVertex, double connectivityThreshold,
      double distanceThreshold, size_t rowIndex_min, size_t rowIndex_max,
      unsigned connectionLength) {
    double two_pi = 2*3.1415926535897931;
    aims::SparseMatrix & conn = *conn_ptr;
    aims::SparseMatrix & conn2 = *conn_ptr2;
    for (QuickMap::const_iterator neigh_2
           = fiberExtremity2NeighMeshVertex.begin();
         neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2) {
      for (QuickMap::const_iterator neigh_1
             = fiberExtremity1NeighMeshVertex.begin();
           neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1) {
        double e1 = square(neigh_1->second);
        double e2 = square(neigh_2->second);
        double square_sigma = distanceThreshold*distanceThreshold;
        double w = exp( -(e1 + e2)  / ( 2*square_sigma ))
                   /(square_sigma*two_pi);//"smoothing coefficient"
        if (w > connectivityThreshold) {
          w = connectionLength *w;
          if (rowIndex_max == 0 or (not conn_ptr2)) {
            conn(neigh_1->first,neigh_2->first)
              = conn(neigh_1->first,neigh_2->first) + w;
            conn(neigh_2->first,neigh_1->first)
              = conn(neigh_2->first,neigh_1->first) + w;
          } else {
            if (neigh_1->first >= rowIndex_min
                and neigh_1->first < rowIndex_max) {
              conn(neigh_1->first,neigh_2->first)
                = conn(neigh_1->first,neigh_2->first) + w;
            } else {
              conn2(neigh_1->first,neigh_2->first)
                = conn2(neigh_1->first,neigh_2->first) + w;
            }
            if (neigh_2->first >= rowIndex_min
                and neigh_2->first < rowIndex_max) {
              conn(neigh_2->first,neigh_1->first)
                = conn(neigh_2->first,neigh_1->first) + w;
            } else {
              conn2(neigh_2->first,neigh_1->first)
                = conn2(neigh_2->first,neigh_1->first) + w;
            }
          }
        }
      }
    }
  }

  //---------------------------------
  //  fillconnMatrixWeightingLength
  //---------------------------------
  /*
  inputs:
      conn_ptr: connectivity matrix to fill, symmetric in this case
      fiberExtremity1NeighMeshVertex: vector<pair<size_t, double> >
      of pairs of:  < i1 = mesh vertex index near the first fiber extremity,
      distance between i1 and first fiber extremity closest vertex in the mesh
      (or intersection between fiber and mesh)>
      fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
      connectivityThreshold: threshold: connectivity contributions below are
      excluded
      distanceThreshold: distance for the smoothing of the connectivity matrix

  output:
      void (*conn_ptr has changed: filled according to input mesh neighborhoods
      (for fiber extrimities))
  */
  void fillconnMatrixWeightingLength(
      Connectivities *conn_ptr, QuickMap &fiberExtremity1NeighMeshVertex,
      QuickMap &fiberExtremity2NeighMeshVertex, double connectivityThreshold,
      double distanceThreshold, float connectionLength) {
    Connectivities &conn = *conn_ptr;
    for (QuickMap::const_iterator neigh_2
           = fiberExtremity2NeighMeshVertex.begin();
         neigh_2 != fiberExtremity2NeighMeshVertex.end(); ++neigh_2) {
      for (QuickMap::const_iterator neigh_1
             = fiberExtremity1NeighMeshVertex.begin();
           neigh_1 != fiberExtremity1NeighMeshVertex.end(); ++neigh_1) {
        double e1 = square(neigh_1->second);
        double e2 = square(neigh_2->second);
        double w = exp( -(e1 + e2)
                    / ( 2*distanceThreshold*distanceThreshold ));
        if (w > connectivityThreshold) {
          float connectionLenghtWeight = connectionLength;
          w = connectionLenghtWeight *w;
          conn[neigh_1->first][neigh_2->first] += w;
          conn[neigh_2->first][neigh_1->first] += w;
        }
      }
    }
  }

  //-----------------------------
  //  fillNonSymetricConnMatrix
  //-----------------------------
  /*
  inputs:
    conn_ptr: connectivity matrix to fill, non symmetric in this case:
    shape = (cortexMeshVertex_nb, thalamusMeshVertex_nb) for example
    fiberExtremity1NeighMeshVertex: vector< pair<size_t, double> > of pairs of:
    < i1 = mesh vertex index near the first fiber extremity, distance between
    i1 and first fiber extremity closest vertex in the mesh (or intersection
    between fiber and mesh)>
    fiberExtremity2NeighMeshVertex: idem with the other fiber extremity
    connectivityThreshold: threshold: connectivity contributions below are
    excluded
    distanceThreshold: distance for the smoothing of the connectivity matrix

  output:
    void (*conn_ptr has changed: filled according to input mesh neighborhoods
    (for fiber extrimities))
  */
  void fillNonSymetricConnMatrix(
      Connectivities *conn_ptr, QuickMap &fiberExtremity1NeighMeshVertex_rows,
      QuickMap & fiberExtremity2NeighMeshVertex_cols,
      double connectivityThreshold, double distanceThreshold) {

    Connectivities &conn = *conn_ptr;
    for (QuickMap::const_iterator neigh_2
           = fiberExtremity2NeighMeshVertex_cols.begin();
         neigh_2 != fiberExtremity2NeighMeshVertex_cols.end(); ++neigh_2) {
      for (QuickMap::const_iterator neigh_1
             = fiberExtremity1NeighMeshVertex_rows.begin();
           neigh_1 != fiberExtremity1NeighMeshVertex_rows.end(); ++neigh_1) {
        double e1 = square(neigh_1->second);
        double e2 = square(neigh_2->second);
        double w = exp( -(e1 + e2)
                   / ( 2*distanceThreshold*distanceThreshold ));
        if (w > connectivityThreshold) {
          conn[neigh_1->first][neigh_2->first] += w;
        }
      }
    }
  }

  //-----------------------------------------------
  //  computeIntersectionPointFiberSegmentAndMesh
  //-----------------------------------------------
  template<int D, class T> bool computeIntersectionPointFiberSegmentAndMesh(
      const AimsTimeSurface<D, T> &aimsMesh,
      const vector<set<unsigned> > &polygonsByVertex_Index,
      Point3df fiberPoint1, Point3df fiberPoint2, unsigned meshClosestPoint,
      QuickMap **polygonVerticesDistMap) {
    //fill polygonVerticesDistMap with the D vertex of the intersection polygon
    //(if there is an intersection)
    const vector<Point3df>  &meshVertices  = aimsMesh.vertex();
    const vector<AimsVector<unsigned,D> > &meshPolygons = aimsMesh.polygon();
    set<unsigned> closest_polygons = polygonsByVertex_Index[meshClosestPoint];
    set<unsigned>::const_iterator set_it;
    AimsVector<unsigned,D> currentAims_polygon;
    vector< Point3df > current_polygon(D);
    Point3df * intersection_point = 0;
    bool inter = false;
    for (set_it = closest_polygons.begin();
         set_it != closest_polygons.end(); set_it++) {
      currentAims_polygon = meshPolygons[*set_it];
      for ( unsigned i = 0; i < D; ++i) {
        current_polygon[i] = meshVertices[currentAims_polygon[i]];
      }
      inter = computeIntersectionSegmentPolygon(
          fiberPoint1, fiberPoint2, current_polygon, &intersection_point);
      if (inter) {
        unsigned currentPolygonVerticeIndex;
        Point3df currentPolygonVertice;
        * polygonVerticesDistMap = new QuickMap(D);
          Point3df intersection_point_na(
            (*intersection_point)[0],
            (*intersection_point)[1],
            (*intersection_point)[2]);

        for  ( unsigned j = 0; j < D; ++j) {
          currentPolygonVerticeIndex = currentAims_polygon[j];
          currentPolygonVertice = meshVertices[currentPolygonVerticeIndex];
          Point3df currentPolygonVertice_na(
              currentPolygonVertice[0],
              currentPolygonVertice[1],
              currentPolygonVertice[2]);
          pair<size_t, double > pair_verticeIndex_distanceToIntersection
            = make_pair(currentAims_polygon[j],
                        sqrt(dist2(currentPolygonVertice_na,
                                   intersection_point_na)));
          (**polygonVerticesDistMap)[j]
            = pair_verticeIndex_distanceToIntersection;
        }
        break;
      }
    }
    delete intersection_point;
    return inter;
  }

  //------------------------------------------------
  //  computeIntersectionPointFiberSegmentAndMesh2
  //------------------------------------------------
  template<int D, class T> bool computeIntersectionPointFiberSegmentAndMesh2(
      const AimsTimeSurface<D, T> &aimsMesh,
      const vector<set<unsigned> > &polygonsByVertex_Index,
      Point3df fiberPoint1, Point3df fiberPoint2, unsigned meshClosestPoint,
      QuickMap **closestPointsTointersectionWeightMap) {
    // fill closestPointsTointersectionWeightMap with the polygon vertice
    // closest to the intersection point, and the weight value = 1.0
    const vector<Point3df>  & meshVertices  = aimsMesh.vertex();
    const vector<AimsVector<unsigned, D> > &meshPolygons = aimsMesh.polygon();
    set<unsigned> closest_polygons = polygonsByVertex_Index[meshClosestPoint];
    set<unsigned>::const_iterator set_it;
    AimsVector<unsigned, D> currentAims_polygon;
    vector<Point3df> current_polygon(D);
    Point3df *intersection_point = 0;
    bool inter = false;
    for (set_it = closest_polygons.begin(); set_it != closest_polygons.end();
         set_it++) {
      currentAims_polygon = meshPolygons[*set_it];
      for (unsigned i = 0; i < D; ++i) {
        current_polygon[i] = meshVertices[currentAims_polygon[i]];
      }
      inter = computeIntersectionSegmentPolygon(fiberPoint1,
                                                fiberPoint2,
                                                current_polygon,
                                                &intersection_point);
      if (inter) {
        unsigned currentPolygonVerticeIndex;
        Point3df currentPolygonVertice;
        * closestPointsTointersectionWeightMap = new QuickMap(1);
        Point3df intersection_point_na(
            (*intersection_point)[0],
            (*intersection_point)[1],
            (*intersection_point)[2]);

        vector <float> polygonVerticesIntersection_distances(D);
        float polygonVerticesIntersection_minDistance = 1000;  //mm
        //unsigned closestVerticeNb_On_D;
        size_t closestVerticeIndex_On_mesh;
        for (unsigned j = 0; j < D; ++j) {
          currentPolygonVerticeIndex = currentAims_polygon[j];
          currentPolygonVertice = meshVertices[currentPolygonVerticeIndex];
          Point3df currentPolygonVertice_na(
              currentPolygonVertice[0],
              currentPolygonVertice[1],
              currentPolygonVertice[2]);
          polygonVerticesIntersection_distances[j] = dist2(
              currentPolygonVertice_na,
              intersection_point_na);
          if (polygonVerticesIntersection_distances[j]
              < polygonVerticesIntersection_minDistance) {
            polygonVerticesIntersection_minDistance
              = polygonVerticesIntersection_distances[j];
            //closestVerticeNb_On_D = j;
            closestVerticeIndex_On_mesh = currentPolygonVerticeIndex;
          }
        }

        pair<size_t, double> pair_verticeIndex_distanceToIntersection
          = make_pair(closestVerticeIndex_On_mesh,1.0);
        (**closestPointsTointersectionWeightMap)[0]
          = pair_verticeIndex_distanceToIntersection;
        break;
      }
    }
    delete intersection_point;
    return inter;
  }

  //-----------------------------------------------------------
  //  computeIntersectionPointNeighborhoodFiberSegmentAndMesh
  //-----------------------------------------------------------
  template<int D, class T>
    bool computeIntersectionPointNeighborhoodFiberSegmentAndMesh(
        const AimsTimeSurface<D, T> &aimsMesh,
        const vector<set<unsigned> > &polygonsByVertex_Index,
        Point3df fiberPoint1, Point3df fiberPoint2, unsigned meshClosestPoint,
        vector<QuickMap> &distanceThresholdNeighborhoodByVertex,
        QuickMap **polygonVerticesDistMap_2ptr) {
    // fill polygonVerticesDistMap with the vertex closest to the intersection
    // point and its neighborhood according to distanceThreshold value, if
    // there is an intersection (a changer: fastmarching a partir des sommets
    // du polygone mieux pour construire le voisinage)
    const vector<Point3df>  & meshVertices  = aimsMesh.vertex();
    const vector< AimsVector<unsigned, D> > &meshPolygons= aimsMesh.polygon();
    set<unsigned> closest_polygons = polygonsByVertex_Index[meshClosestPoint];
    set<unsigned>::const_iterator set_it;
    AimsVector<unsigned,D> currentAims_polygon;
    vector< Point3df > current_polygon(D);
    Point3df * intersection_point = 0;
    bool inter = false;
    for (set_it = closest_polygons.begin(); set_it != closest_polygons.end();
         set_it++) {
      currentAims_polygon = meshPolygons[*set_it];
      for (unsigned i = 0; i < D; ++i) {
        current_polygon[i] = meshVertices[currentAims_polygon[i]];
      }
      inter = computeIntersectionSegmentPolygon(fiberPoint1,
                                                fiberPoint2,
                                                current_polygon,
                                                &intersection_point);
      if (inter) {
        unsigned currentPolygonVerticeIndex;
        Point3df currentPolygonVertice;
        * polygonVerticesDistMap_2ptr = new QuickMap(D);
        Point3df intersection_point_na(
            (*intersection_point)[0],
            (*intersection_point)[1],
            (*intersection_point)[2]);

        constel::QuickMap & polygonVerticesDistMap
          = **polygonVerticesDistMap_2ptr;
        vector <float> polygonVerticesIntersection_distances(D);
        float polygonVerticesIntersection_minDistance = 1000;  //mm
        unsigned closestVerticeNb_On_D = 0;
        for  (unsigned j = 0; j < D; ++j) {
          currentPolygonVerticeIndex = currentAims_polygon[j];
          currentPolygonVertice = meshVertices[currentPolygonVerticeIndex];
          Point3df currentPolygonVertice_na(
              currentPolygonVertice[0],
              currentPolygonVertice[1],
              currentPolygonVertice[2]);
          polygonVerticesIntersection_distances[j]= dist2(
              currentPolygonVertice_na,
              intersection_point_na);
          if (polygonVerticesIntersection_distances[j]
              < polygonVerticesIntersection_minDistance) {
            polygonVerticesIntersection_minDistance
              = polygonVerticesIntersection_distances[j];
            closestVerticeNb_On_D = j;
          }
        }
        polygonVerticesDistMap = distanceThresholdNeighborhoodByVertex[
            currentAims_polygon[closestVerticeNb_On_D]];
        break;
      }
    }
    delete intersection_point;
    return inter;
  }

  //---------------------------------
  //  fillconnMatrixWithConnections
  //---------------------------------
  void fillconnMatrixWithConnections(
      Connectivities *conn_ptr, const BundleConnections &connections,
      double connectivityThreshold, double distanceThreshold) {
    size_t connectionsCount = 0;
    if (distanceThreshold > 0.0) {
      for (BundleConnections::const_iterator iConnection = connections.begin();
           iConnection != connections.end();
           ++iConnection, ++connectionsCount) {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr,
                       ConnectionMapFront,
                       ConnectionMapBack,
                       connectivityThreshold,
                       distanceThreshold);
      }
    } else {
      for (BundleConnections::const_iterator iConnection = connections.begin();
           iConnection != connections.end();
           ++iConnection, ++connectionsCount) {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrixNoSmoothing(conn_ptr,
                                  ConnectionMapFront,
                                  ConnectionMapBack);
      }
    }
  }

  //---------------------------------
  //  fillconnMatrixWithConnections
  //---------------------------------
  void fillconnMatrixWithConnections(
      aims::SparseMatrix *conn_ptr, const BundleConnections &connections,
      double connectivityThreshold, double distanceThreshold,
      size_t rowIndex_min, size_t rowIndex_max,
      aims::SparseMatrix *conn_ptr2) {
    size_t connectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin();
         iConnection != connections.end(); ++iConnection, ++connectionsCount) {
      constel::QuickMap ConnectionMapFront = iConnection->front();
      constel::QuickMap ConnectionMapBack = iConnection->back();
      if (not conn_ptr2) {
        fillconnMatrix(conn_ptr,
                       ConnectionMapFront,
                       ConnectionMapBack,
                       connectivityThreshold,
                       distanceThreshold);
      } else {
        fillconnMatrix(conn_ptr,
                       conn_ptr2,
                       ConnectionMapFront,
                       ConnectionMapBack,
                       connectivityThreshold,
                       distanceThreshold,
                       rowIndex_min,
                       rowIndex_max);
      }
    }
  }


  //-------------------------------------------
  //  fillconnMatrixWithConnectionsPlusLength
  //-------------------------------------------
  void fillconnMatrixWithConnectionsPlusLength(
      Connectivities *conn_ptr, const BundleConnections &connections,
      double connectivityThreshold, double distanceThreshold,
      unsigned length_min, unsigned length_max,
      ConnectionsLength &connectionsLength) {
    size_t connectionsCount = 0;
    size_t chosenConnectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin();
         iConnection != connections.end(); ++iConnection, ++connectionsCount) {
      unsigned connection_length = connectionsLength[connectionsCount];
      if (connection_length < length_max
          and connection_length >= length_min) {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr,
                       ConnectionMapFront,
                       ConnectionMapBack,
                       connectivityThreshold,
                       distanceThreshold);
        ++chosenConnectionsCount;
      }
    }
  }

  //-------------------------------------------------
  //  fillconnMatrixWithConnectionsPlusLengthWeight
  //-------------------------------------------------
  void fillconnMatrixWithConnectionsPlusLengthWeight(
      Connectivities *conn_ptr, const BundleConnections &connections,
      double connectivityThreshold, double distanceThreshold,
      unsigned length_min, unsigned length_max,
      ConnectionsLength &connectionsLength) {
    size_t connectionsCount = 0;
    size_t chosenConnectionsCount = 0;

    for (BundleConnections::const_iterator iConnection = connections.begin();
         iConnection != connections.end(); ++iConnection, ++connectionsCount) {
      unsigned connection_length = connectionsLength[connectionsCount];

      if (connection_length < length_max
          and connection_length >= length_min) {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr,
                       ConnectionMapFront,
                       ConnectionMapBack,
                       connectivityThreshold,
                       distanceThreshold,
                       connection_length);
        ++chosenConnectionsCount;
      }
    }
  }

  //------------------------------------------------------
  //  fillconnMatrixWithConnectionsPlusFloatLengthWeight
  //------------------------------------------------------
  void fillconnMatrixWithConnectionsPlusFloatLengthWeight(
      Connectivities *conn_ptr, const BundleConnections &connections,
      double connectivityThreshold, double distanceThreshold,
      float length_min, float length_max,
      ConnectionsFloatLength &connectionsLength) {
    size_t connectionsCount = 0;
    size_t chosenConnectionsCount = 0;

    for (BundleConnections::const_iterator iConnection = connections.begin();
         iConnection != connections.end(); ++iConnection, ++connectionsCount) {
      float connection_length = connectionsLength[connectionsCount];

      if (connection_length < length_max
          and connection_length >= length_min) {
        constel::QuickMap ConnectionMapFront = iConnection->front();
        constel::QuickMap ConnectionMapBack = iConnection->back();
        fillconnMatrix(conn_ptr,
                       ConnectionMapFront,
                       ConnectionMapBack,
                       connectivityThreshold,
                       distanceThreshold,
                       connection_length);
        ++chosenConnectionsCount;
      }
    }
  }

  //--------------------------------------------
  //  fillNonSymetricConnMatrixWithConnections
  //--------------------------------------------
  void fillNonSymetricConnMatrixWithConnections(
      Connectivities *conn_ptr, const BundleConnections &connections,
      double connectivityThreshold, double distanceThreshold) {
    size_t connectionsCount = 0;
    for (BundleConnections::const_iterator iConnection = connections.begin();
         iConnection != connections.end(); ++iConnection, ++connectionsCount) {
      constel::QuickMap ConnectionMapFront = iConnection->front();
      constel::QuickMap ConnectionMapBack = iConnection->back();
      fillNonSymetricConnMatrix(conn_ptr,
                                ConnectionMapFront,
                                ConnectionMapBack,
                                connectivityThreshold,
                                distanceThreshold);
    }
  }

  //-----------------------------------------------
  //  computeIntersectionPointFiberSegmentAndMesh
  //-----------------------------------------------
  template bool computeIntersectionPointFiberSegmentAndMesh(
      const AimsSurfaceTriangle &aimsMesh,
      const vector<set<unsigned> > &polygonsByVertex_Index,
      Point3df fiberPoint1, Point3df fiberPoint2, unsigned meshClosestPoint,
      QuickMap **polygonVerticesDistMap);

  //------------------------------------------------
  //  computeIntersectionPointFiberSegmentAndMesh2
  //------------------------------------------------
  template bool computeIntersectionPointFiberSegmentAndMesh2(
      const AimsSurfaceTriangle &aimsMesh,
      const vector<set<unsigned> > &polygonsByVertex_Index,
      Point3df fiberPoint1, Point3df fiberPoint2, unsigned meshClosestPoint,
      QuickMap **polygonVerticesDistMap);

  //-----------------------------------------------------------
  //  computeIntersectionPointNeighborhoodFiberSegmentAndMesh
  //-----------------------------------------------------------
  template bool computeIntersectionPointNeighborhoodFiberSegmentAndMesh(
      const AimsSurfaceTriangle &aimsMesh,
      const vector<set<unsigned> > &polygonsByVertex_Index,
      Point3df fiberPoint1, Point3df fiberPoint2, unsigned meshClosestPoint,
      vector<QuickMap> &distanceThresholdNeighborhoodByVertex,
      QuickMap **polygonVerticesDistMap_2ptr);

} // namespace constel
