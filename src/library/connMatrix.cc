#include <constellation/connMatrix.h>
// includes from CATHIER
#include <cathier/triangle_mesh_geodesic_map.h>
#include <aims/sparsematrix/sparseordensematrix.h>


using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;
using namespace boost;

namespace constel
{

  Connectivities * connMatrix( const Fibers & fibers,
                               const AimsSurfaceTriangle & inAimsMesh,
                               float distthresh, float wthresh, Motion motion,
                               bool verbose )
  {
    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity contributions
    
    //Convert inAimsMesh to Mesh mesh (Pascal Cathier Format)
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);
    
    // Comuting geomap : neighborhood map
  
    if (verbose) std::cout << "Computing geomap..." << std::flush;
    if (verbose) std::cout << "step 1..." << std::flush;
    std::vector<QuickMap> res(getVertices(mesh).size());
    if (distthresh!=0)
    {
      til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
      boost::shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
      til::Triangle_mesh_geodesic_map<Mesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
          geomap(getVertices(mesh), *pneighc, stopGhost);
      std::vector<std::size_t> startPoints(1);
      std::vector<double> dist(1, 0.0);
//       std::vector<QuickMap> res(getVertices(mesh).size());
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
    }
    if (verbose) std::cout << "OK" << std::endl;
  
    if (verbose) std::cout << "Number of fibers: " << fibers.size() << std::endl;

    // Generating kdtree
//     std::cout << "Generating kdtree" << std::endl;
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
//     std::cout << "getVertices(mesh):" << getVertices(mesh)[0] << ", " << getVertices(mesh)[1]  << std::endl;
    // Computing connectivity matrix
    if (verbose) std::cout << "Computing connectivity matrix" << std::endl;
    Connectivities * conn_ptr = new Connectivities(getVertices(mesh).size(), til::SparseVector<double>(getVertices(mesh).size()));
    Connectivities & conn  = *conn_ptr;
    //Connectivities conn(getVertices(mesh).size(), til::SparseVector<double>(getVertices(mesh).size()));
    int countRemoved = 0;
    std::size_t fiberCount = 0, nFibers = fibers.size();
    std::size_t five_count = nFibers / 20.0;
    double total_conn = 0;
    Point3df p1, p2;
    for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
//     if (fiberCount % int(nFibers/10.0) == 0)
//     {
//        std::cout << 10 * int((fiberCount / int(nFibers/10.0))) << "%..." << std::flush;
//     }
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          if (verbose) std::cout << 5 * int(fiberCount / five_count) << "%..." << std::flush;
        }
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
      if(distthresh==0)
      {
        conn[A][B] += 1.;
        conn[B][A] += 1.;
      }
      else
      {
        double wtotal = 0;
        list<pair<pair<size_t, size_t>, double> > weights;
        for (QuickMap::const_iterator B_neigh = res[B].begin(); B_neigh != res[B].end(); ++B_neigh)
        {
          for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
          {
            double square_sigma = distthresh*distthresh;
            double e1 = til::square(A_neigh->second);
            double e2 = til::square(B_neigh->second);
            double w1 = std::exp( -e1 / ( 2*square_sigma ))/(square_sigma * two_pi);//"smoothing coefficient"
            double w2 = std::exp( -e2 / ( 2*square_sigma ))/(square_sigma * two_pi);//"smoothing coefficient"
            if( w1 < G_THRESH )
              w1 = 0;
            if( w2 < G_THRESH )
              w2 = 0;
            double w = w1 + w2;
            if( w > 0 )
            {
              wtotal += w;
              weights.push_back( make_pair( make_pair( A_neigh->first, B_neigh->first ), w ) );
            }
            // double w = std::exp( -(e1 + e2)  / ( 2*distthresh*distthresh ))/(square_sigma * two_pi);//"smoothing coefficient"
          }
        }
        
        list<pair<pair<size_t, size_t>, double> >::iterator iw, ew = weights.end();
        for ( iw=weights.begin(); iw!=ew; ++iw )
        {
          double w = iw->second / wtotal;
          conn[iw->first.second][iw->first.first] += w;
          conn[iw->first.first][iw->first.second] += w;
          
          total_conn += 2 * w;
        }
      }
    }
    if (verbose) std::cout << "...100%..." << std::flush;
    if (verbose) std::cout << "Done." << std::endl;
    if (verbose) std::cout << "NB: " << countRemoved << " fibers out of " << fibers.size() << " have been discarded" << std::endl;
  
    // To compress a little bit, remove points that are below some threshold
    double nonzero_count=0;
    if (verbose) std::cout << "Removing weak connections: " << std::flush;
    {
      int count1 = 0;
      int count2 = 0;
      for (Connectivities::iterator i = conn.begin(); i != conn.end(); ++i)
      {
        for (Connectivity::Map::iterator j = i->getMap().begin(); j != i->getMap().end(); )
        {
          ++count1;
          if (j->second <= wthresh)
          {
            ++count2;
            i->getMap().erase(j++);
          }
          else
          {
            ++j;
            ++nonzero_count;
          }
        }
      }
      if (verbose) std::cout << "Removed " << count2 << " elements out of " << count1 << std::endl;
    }
//     double mean_conn = total_conn / nonzero_count;
//     std::cout << "total_conn : " << total_conn << " on "<< nonzero_count << " non zero elements" << endl;
//     std::cout << "mean connectivity for non zero elements : " << mean_conn << endl;
    return ( conn_ptr);
  }// connMatrix


  Connectivities * connMatrixSeedMeshToTargetMesh(
    const Fibers & fibers, const AimsSurfaceTriangle & aimsSeedMesh,
    const AimsSurfaceTriangle & aimsTargetMesh, float distthresh,
    float wthresh, Motion motion, bool verbose )
  {
    /*
    input:
                fibers
                aimsSeedMesh
                aimsTargetMesh
                distthresh: distance for the smoothing of the matrix
                wthresh: weight threshold for clean the connectivity matrix (values below are removed)
                motion: transformation from t2 (fibers) to anat(meshes)

    output:
                conn: sparse matrix of shape = [seedMeshVertex_nb,targetMeshVertex_nb]
    */
    const double G_THRESH = 0.001; //threshold for connectivity contributions
    
    //Convert inAimsMesh to Mesh mesh (Pascal Cathier Format)
    Mesh seedMesh;
    til::Mesh1 seedMesh0;
    til::convert(seedMesh0, aimsSeedMesh);
    seedMesh = addNeighborsToMesh(seedMesh0);
    Mesh targetMesh;
    til::Mesh1 targetMesh0;
    til::convert(targetMesh0, aimsTargetMesh);
    targetMesh = addNeighborsToMesh(targetMesh0);
    // Comuting geomap : neighborhood map
    
    //For seedMesh:
    if (verbose) std::cout << "Computing geomap seedMesh..." << std::flush;
    if (verbose) std::cout << "step 1..." << std::flush;
    std::vector<QuickMap> res_seedMesh(getVertices(seedMesh).size());
    
    if (distthresh!=0)
    {
      til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
      boost::shared_ptr<CNeighborhoods> pneighc_seedMesh = til::circular_neighborhoods(getVertices(seedMesh), getFaceIndices(seedMesh));
      til::Triangle_mesh_geodesic_map<Mesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
          geomap_seedMesh(getVertices(seedMesh), *pneighc_seedMesh, stopGhost);
      std::vector<std::size_t> startPoints(1);
      std::vector<double> dist(1, 0.0);
//       std::vector<QuickMap> res(getVertices(mesh).size());
      std::vector<std::size_t> nneigh_seedMesh(til::size(getVertices(seedMesh)));
      
      {
        for (std::size_t i = 0; i < til::size(getVertices(seedMesh)); ++i)
        {
          startPoints[0] = i;
          geomap_seedMesh.init(startPoints, dist);
          geomap_seedMesh.process();
          boost::shared_ptr<til::sparse_vector<double> > tmp = geomap_seedMesh.distanceMap();
          res_seedMesh[i].resize(tmp->getMap().size());
          {
            using namespace til::expr;
            til::detail::loop_xx(castTo(*_1, *_2), res_seedMesh[i], tmp->getMap());
          }
  
          nneigh_seedMesh[i] = res_seedMesh[i].size();
        }
      }
    }
    
    //For targetMesh:
    if (verbose) std::cout << "Computing geomap targetMesh..." << std::flush;
    if (verbose) std::cout << "step 1..." << std::flush;
    std::vector<QuickMap> res_targetMesh(getVertices(targetMesh).size());
    
    if (distthresh!=0)
    {
      til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
      boost::shared_ptr<CNeighborhoods> pneighc_targetMesh = til::circular_neighborhoods(getVertices(targetMesh), getFaceIndices(targetMesh));
      til::Triangle_mesh_geodesic_map<Mesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
          geomap_targetMesh(getVertices(targetMesh), *pneighc_targetMesh, stopGhost);
      std::vector<std::size_t> startPoints(1);
      std::vector<double> dist(1, 0.0);
//       std::vector<QuickMap> res(getVertices(mesh).size());
      std::vector<std::size_t> nneigh_targetMesh(til::size(getVertices(targetMesh)));
      
      {
        for (std::size_t i = 0; i < til::size(getVertices(targetMesh)); ++i)
        {
          startPoints[0] = i;
          geomap_targetMesh.init(startPoints, dist);
          geomap_targetMesh.process();
          boost::shared_ptr<til::sparse_vector<double> > tmp = geomap_targetMesh.distanceMap();
          res_targetMesh[i].resize(tmp->getMap().size());
          {
            using namespace til::expr;
            til::detail::loop_xx(castTo(*_1, *_2), res_targetMesh[i], tmp->getMap());
          }
  
          nneigh_targetMesh[i] = res_targetMesh[i].size();
        }
      }
    }
    
    if (verbose) std::cout << "OK" << std::endl;
  
    if (verbose) std::cout << "Number of fibers: " << fibers.size() << std::endl;

    // Generating kdtree
    //For seedMesh:
    if (verbose) std::cout << "Generating kdtrees" << std::endl;
    KDTree kdt_seedMesh(getVertices(seedMesh));
    makeKDTree(getVertices(seedMesh), kdt_seedMesh);
    //For targetMesh:
    KDTree kdt_targetMesh(getVertices(targetMesh));
    makeKDTree(getVertices(targetMesh), kdt_targetMesh);
    
    // Computing connectivity matrix
    if (verbose) std::cout << "Computing connectivity matrix" << std::endl;
    Connectivities * conn_ptr = new Connectivities(getVertices(seedMesh).size(), til::SparseVector<double>(getVertices(targetMesh).size()));
    Connectivities & conn  = *conn_ptr;
    //Connectivities conn(getVertices(mesh).size(), til::SparseVector<double>(getVertices(mesh).size()));
    int countRemoved = 0;
    std::size_t fiberCount = 0, nFibers = fibers.size();
    std::size_t five_count = nFibers / 20.0;
    double total_conn = 0;
    Point3df p1, p2;
    for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
//     if (fiberCount % int(nFibers/10.0) == 0)
//     {
//        std::cout << 10 * int((fiberCount / int(nFibers/10.0))) << "%..." << std::flush;
//     }
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          if (verbose) std::cout << 5 * int(fiberCount / five_count) << "%..." << std::flush;
        }
      }
      //Looking for closest points
      til::Find_closest< double, KDTree > fc_seedMesh(kdt_seedMesh);
      til::Find_closest< double, KDTree > fc_targetMesh(kdt_targetMesh);
      // Transform iFiber->front() and iFiber->back() from t2 to anat space (with motion)
      p1 = motion.transform(iFiber->front()[0],iFiber->front()[1],iFiber->front()[2]);
      til::numeric_array<float, 3> p1na(p1[0], p1[1], p1[2]);
      p2 = motion.transform(iFiber->back()[0],iFiber->back()[1],iFiber->back()[2]);
      til::numeric_array<float, 3> p2na(p2[0], p2[1], p2[2]);
      std::size_t A_seedMesh = fc_seedMesh(p1na);
      std::size_t B_seedMesh = fc_seedMesh(p2na);
      std::size_t A_targetMesh = fc_targetMesh(p1na);
      std::size_t B_targetMesh = fc_targetMesh(p2na);
      
      // Filling the connectivity matrix
      /*Two retained cases : first case: p1na "touches" seedMesh and p2na "touches" targetMesh
                              second case: the opposite (p1na "touches" targetMesh and p2na "touches" seedMesh
      */
      if (til::dist2(p1na, getVertices(seedMesh)[A_seedMesh], til::prec<float>()) <= 25.0) //p1na "touches" seedMesh
      {
        if (til::dist2(p2na, getVertices(targetMesh)[B_targetMesh], til::prec<float>()) <= 25.0) //first case (p2na "touches" targetMesh)
        {
          if(distthresh==0)
          {
            conn[A_seedMesh][B_targetMesh] += 1.;
          }
          else
          {
            for (QuickMap::const_iterator A_seedMesh_neigh = res_seedMesh[A_seedMesh].begin(); A_seedMesh_neigh != res_seedMesh[A_seedMesh].end(); ++A_seedMesh_neigh)
            {
              for(QuickMap::const_iterator B_targetMesh_neigh = res_targetMesh[B_targetMesh].begin(); B_targetMesh_neigh != res_targetMesh[B_targetMesh].end(); ++B_targetMesh_neigh)
              {
                double e1 = til::square(A_seedMesh_neigh->second);
                double e2 = til::square(B_targetMesh_neigh->second);
                double w = std::exp( -(e1 + e2)  / ( 2*distthresh*distthresh ));//"smoothing coefficient"
                if (w > G_THRESH)
                {
                  conn[A_seedMesh_neigh->first][B_targetMesh_neigh->first] += w;
                  total_conn +=  w;
                }
              }
            }
          }
        }
        else // p2na doesn t touch targetMesh
        {
          ++countRemoved;
          continue;
        }
      }
      else
      {
        if (til::dist2(p1na, getVertices(targetMesh)[A_targetMesh], til::prec<float>()) <= 25.0) //p1na "touches" targetMesh
        {
          if (til::dist2(p2na, getVertices(seedMesh)[B_seedMesh], til::prec<float>()) <= 25.0) //second case (p2na "touches" seedMesh)
          {
            if(distthresh==0)
            {
              conn[B_seedMesh][A_targetMesh] += 1.;
            }
            else
            {
              for (QuickMap::const_iterator A_targetMesh_neigh = res_targetMesh[A_targetMesh].begin(); A_targetMesh_neigh != res_targetMesh[A_targetMesh].end(); ++A_targetMesh_neigh)
              {
                for(QuickMap::const_iterator B_seedMesh_neigh = res_seedMesh[B_seedMesh].begin(); B_seedMesh_neigh != res_seedMesh[B_seedMesh].end(); ++B_seedMesh_neigh)
                {
                  double e1 = til::square(A_targetMesh_neigh->second);
                  double e2 = til::square(B_seedMesh_neigh->second);
                  double w = std::exp( -(e1 + e2)  / ( 2*distthresh*distthresh ));//"smoothing coefficient"
                  if (w > G_THRESH)
                  {
                    conn[B_seedMesh_neigh->first][A_targetMesh_neigh->first] += w;
                    total_conn +=  w;
                  }
                }
              }
            }
          }
          else // p2na doesn t touch seedMesh
          {
            ++countRemoved;
            continue;
          }
        }
        else // p1na touches nor seedMesh neither targetMesh
        {
          ++countRemoved;
          continue;
        }
      }
    }
    
    if (verbose) std::cout << "...100%..." << std::flush;
    
    if (verbose) std::cout << "Done." << std::endl;
    if (verbose) std::cout << "NB: " << countRemoved << " fibers out of " << fibers.size() << " have been discarded" << std::endl;
  
    // To compress a little bit, remove points that are below some threshold
    double nonzero_count=0;
    if (verbose) std::cout << "Removing weak connections: " << std::flush;
    {
      int count1 = 0;
      int count2 = 0;
      for (Connectivities::iterator i = conn.begin(); i != conn.end(); ++i)
      {
        for (Connectivity::Map::iterator j = i->getMap().begin(); j != i->getMap().end(); )
        {
          ++count1;
          if (j->second <= wthresh)
          {
            ++count2;
            i->getMap().erase(j++);
          }
          else
          {
            ++j;
            ++nonzero_count;
          }
        }
      }
      if (verbose) std::cout << "Removed " << count2 << " elements out of " << count1 << std::endl;
    }
    return (conn_ptr);
  }//connMatrixMesh1ToMesh2


  Connectivity * connMatrixToRois(
    const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh,
    float distthresh, float wthresh, Motion motion, bool verbose )
  {
    if (verbose) std::cout << "Computing connectivity Matrix to Rois" << std::endl;
    /*
    Coumputing connectivity matrix between the mesh vertex and ROIs (subcortical structures for example):
    for the moment: on suppose que les fibres ont ete filtres en fonction d une ROI et pour chaque element du maillage, on compte seulement le nombre de fibres qu'il y a...
    */
    const double G_THRESH = 0.001; //threshold for connectivity contributions
    
    //Convert inAimsMesh to Mesh mesh (Pascal Cathier Format)
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);
    
    // Comuting geomap : neighborhood map
    if (verbose) std::cout << "Computing geomap..." << std::flush;
    til::ghost::GMapStop_AboveThreshold<double> stopGhost(distthresh);
    boost::shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(mesh), getFaceIndices(mesh));
    til::Triangle_mesh_geodesic_map<Mesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
        geomap(getVertices(mesh), *pneighc, stopGhost);
    std::vector<std::size_t> startPoints(1);
    std::vector<double> dist(1, 0.0);
    std::vector<QuickMap> res(getVertices(mesh).size());
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
    if (verbose) std::cout << "OK" << std::endl;
    
    if (verbose) std::cout << "Number of fibers: " << fibers.size() << std::endl;

    // Generating kdtree
//     if (verbose) std::cout << "Generating kdtree" << std::endl;
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    
    // Computing connectivity matrix
    if (verbose) std::cout << "Computing connectivity matrix" << std::endl;
    Connectivity * conn_ptr = new Connectivity(getVertices(mesh).size());
    Connectivity & conn  = *conn_ptr;
    int countRemoved = 0;
    std::size_t fiberCount = 0, nFibers = fibers.size();
    std::size_t five_count = nFibers / 20.0;
    double total_conn = 0;
    Point3df p1, p2;
    for (Fibers::const_iterator iFiber = fibers.begin(); iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
//     if (fiberCount % int(nFibers/10.0) == 0)
//     {
//        std::cout << 10 * int((fiberCount / int(nFibers/10.0))) << "%..." << std::flush;
//     }
      if (fiberCount % five_count == 0)
      {
        std::cout << 5 * int(fiberCount / five_count) << "%..." << std::flush;
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

      if (til::dist2(p1na, getVertices(mesh)[A], til::prec<float>()) > 25.0 and til::dist2(p2na, getVertices(mesh)[B], til::prec<float>()) > 25.0)
      {
        ++countRemoved;
        continue;
      }
      // Filling the connectivity matrix
      if (til::dist2(p1na, getVertices(mesh)[A], til::prec<float>()) <= 25.0 )
      {
        if(distthresh==0)
        {
          conn[A] += 1.;
        }
        else
        {
          for(QuickMap::const_iterator A_neigh = res[A].begin(); A_neigh != res[A].end(); ++A_neigh)
          {
            double e = til::square(A_neigh->second);
            double w = std::exp( -e  / ( 2*distthresh*distthresh ));//"smoothing coefficient"
            if (w > G_THRESH)
            {
              conn[A_neigh->first] += w;
            }
          }
        }
      }
      if (til::dist2(p2na, getVertices(mesh)[B], til::prec<float>()) <= 25.0 )
      {
        if(distthresh==0)
        {
          conn[A] += 1.;
        }
        else
        {
          for(QuickMap::const_iterator B_neigh = res[B].begin(); B_neigh != res[B].end(); ++B_neigh)
          {
            double e = til::square(B_neigh->second);
            double w = std::exp( -e  / ( 2*distthresh*distthresh ));//"smoothing coefficient"
            if (w > G_THRESH)
            {
              conn[B_neigh->first] += w;
            }
          }
        }
      }
    }
    if (verbose) std::cout << std::endl << "NB: " << countRemoved << " fibers out of " << fibers.size() << " have been discarded" << std::endl;
    return (conn_ptr);
  }//connMatrixToRois


  TimeTexture<float> * connMatrixRow_TO_TimeTexture_FLOAT(
    Connectivity * conn_ptr )
  {
    Connectivity & conn = *conn_ptr;
    TimeTexture<float> * connTex_ptr = new TimeTexture<float>;
    TimeTexture<float> & connTex = *connTex_ptr;
    std::size_t conn_size = conn.size();
    connTex.reserve(conn_size);
    for (std::size_t i = 0; i < conn_size; i++)
    {
      connTex[0].push_back(conn[i]);
    }
    return connTex_ptr;
  }
  
  Connectivity * connMatrixSumRows(Connectivities * matrix_ptr, bool verbose)
  {
    /*
    Compute the sum of the matrix rows
    input: Connectivities matrix_ptr, shape = (r,l)
    output: Connectivity matrixSumRows, shape = (1,l)
    
    */
    
    Connectivities & matrix = *matrix_ptr;
    std::size_t colNb = matrix[0].size();
    std::size_t rowsNb = matrix.size();
    if (verbose) std::cout << "("<<rowsNb << "," << colNb << ") : " << std::flush;
    Connectivity sparseVectorIdentity = Connectivity(colNb);
    Connectivity * matrixSumRows_ptr = new Connectivity(colNb);//default values = zero
    Connectivity & matrixSumRows  = *matrixSumRows_ptr;
    for(std::size_t j = 0; j < colNb; ++j)
    {
      sparseVectorIdentity[j]=1.;
      matrixSumRows[j]=0;
    }
//     int five_count = int(rowsNb / 20.0);
    for (std::size_t i = 0; i < rowsNb; ++i)
    {
//       if (i % five_count == 0)
//       {
//         std::cout << 5 * int(i / five_count) << "%..." << std::flush;
//       }
      Connectivity line_i = matrix[i];
      if (line_i.is_null()!=1)
      {
        matrixSumRows += line_i;
      }
    }
    return ( matrixSumRows_ptr );
  }//connMatrixSumRows


  vector<double> * connMatrixSumRows( const SparseOrDenseMatrix & matrix, bool verbose )
  {
    /*
    Compute the sum of the matrix rows
    input: Connectivities matrix_ptr, shape = (r,l)
    output: Connectivity matrixSumRows, shape = (1,l)

    */

    std::size_t colNb = matrix.getSize2();
    std::size_t rowsNb = matrix.getSize1();
    if (verbose) std::cout << "("<<rowsNb << "," << colNb << ") : " << std::flush;
    vector<double> * matrixSumRows_ptr = new vector<double>(colNb);//default values = zero
    vector<double> & matrixSumRows  = *matrixSumRows_ptr;
    for(std::size_t j = 0; j < colNb; ++j)
    {
      matrixSumRows[j]=0;
    }
//     int five_count = int(rowsNb / 20.0);
    for (std::size_t i = 0; i < rowsNb; ++i)
    {
//       if (i % five_count == 0)
//       {
//         std::cout << 5 * int(i / five_count) << "%..." << std::flush;
//       }
      vector<double> line_i = matrix.getRow(i);
      for( int j=0; j<colNb; ++j )
      {
        matrixSumRows[j] += line_i[j];
      }
    }
    return ( matrixSumRows_ptr );
  }//connMatrixSumRows


  Connectivities * connMatrixTargetsRegroup(
    Connectivities * connMatrixToAllMesh_ptr,
    const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
    bool verbose )
  {
    /*  Computing connectivity matrix regrouped according to targetRegions given in targetRegionsTexture (labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to targetRegionsNb = target regions
    input: Connectivities * connMatrixToAllMesh_ptr, shape = (meshVertex_nb,meshVertexNb)
    output: Connectivities * regroupConn_ptr, shape = (targetRegions Nb, targetRegionsNb)
    */
    if (verbose) std::cout << "Start connMatrixTargetsRegroup" << std::endl;
    Connectivities & connMatrixToAllMesh  = *connMatrixToAllMesh_ptr;
    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    std::size_t meshVertexNb = connMatrixToAllMesh[0].size();//number of columns
    std::size_t seedRegionsNb = connMatrixToAllMesh.size();//number of rows
    
    if (meshVertexNb!=seedRegionsNb or targetTex0.nItem()!=meshVertexNb)
    {
      if (verbose) std::cout << "error in matrix dimensions" << std::endl;
    }
    Connectivities * regroupConn_ptr = new Connectivities(targetRegionsNb, til::SparseVector<double>(targetRegionsNb));//default values = zero
    Connectivities & regroupConnMatrix  = *regroupConn_ptr;
    
    if (verbose) std::cout << "Writing connectivity matrix between the target regions: (" << targetRegionsNb << ","  << targetRegionsNb << ") matrix."<< std::endl;
    int ten_count = int(meshVertexNb / 10.0);
    
    for (std::size_t i = 0; i < meshVertexNb; i++)
    {
      int label = targetTex0[i];//equivalent to targetRegionsTex[0].item(i)
      if (1 <= label && label <= targetRegionsNb)
      {
        if (i % ten_count == 0)
        {
          if (verbose) std::cout << 10 * int(i / ten_count) << "%..." << std::flush;
        }
  
        for (std::size_t j = 0; j < i; ++j)
        {
          int label2 = targetTex0[j];  // target region (label) is connected to target region (label2)
          if (1 <= label2 && label2 <= targetRegionsNb && i != j)
          {
            regroupConnMatrix[label-1][label2-1] += connMatrixToAllMesh[i][j];
            regroupConnMatrix[label2-1][label-1] += connMatrixToAllMesh[i][j];
//             float connMatrix_value = connMatrixToAllMesh[i][j];
//             if (connMatrix_value!=0)
//             {
//               std::cout << "(i,j): , connMatrix_value:" << connMatrix_value << std::endl;
//             }
          }
        }
      }
    }
    if (verbose) std::cout << "Done."<<std::endl;
    return regroupConn_ptr;
  }//connMatrixTargetsRegroup



  Connectivities * connMatrixRegionExtract(
    const SparseOrDenseMatrix & AllMeshconnMatrixToAllMesh,
    const TimeTexture<short> & seedRegionsTex, int seedRegionLabel,
    std::size_t seedRegionLabelVertexNb,
    std::vector<std::size_t> ** seedVertexIndex, bool verbose )
  {
    if (verbose) cout << "Start connMatrixRegionExtract" << endl;
    /* Write connectivity matrix of the vertex of a given seed region (label = seedRegionLabel)
    input: Connectivities * connMatrixToAllMesh_ptr, shape = (meshVertexNb, other_meshVertex_nb)
    output: Connectivities * extractConn_ptr, shape = (seedRegionVertexNb, other_meshVertexNb)
    */
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    size_t meshVertexNb = AllMeshconnMatrixToAllMesh.getSize2();//number of columns
    size_t seedRegionsNb = AllMeshconnMatrixToAllMesh.getSize1();//number of rows
    size_t seedRegionsMeshVertexNb  = seedTex0.nItem();
    if (verbose) cout << "seedRegionsMeshVertexNb(rows):" << seedRegionsMeshVertexNb
      << ", meshVertexNb(cols):" << meshVertexNb << endl;
    
    if (verbose) std::cout << "seedRegionLabelVertexNb:"
      << seedRegionLabelVertexNb << std::endl;
    Connectivities * extractConn_ptr = new Connectivities(seedRegionLabelVertexNb, til::SparseVector<double>(meshVertexNb));//default values = zero
    Connectivities & extractConnMatrix  = *extractConn_ptr;
    * seedVertexIndex = new std::vector<size_t>;
    (** seedVertexIndex).resize(seedRegionLabelVertexNb);
    size_t vertexCount = 0;

    int five_count = int(seedRegionsMeshVertexNb / 20.0);
    for( size_t i = 0; i < seedRegionsMeshVertexNb; ++i)//iteration on the rows
    {
//       std::cout << "i:" << i << std::endl;
      if( i % five_count == 0 )
      {
        if( verbose ) std::cout << 5 * int(i / five_count) << "%..." << std::flush;
      }
      if( seedTex0[i]==seedRegionLabel )
      {
        if( vertexCount >= seedRegionLabelVertexNb )
        {
          if(verbose) std::cout <<"vertexCount >= seedRegionLabelVertexNb : "
            << vertexCount << ">= " << seedRegionLabelVertexNb << std::endl;
          throw runtime_error( "error in seedRegionLabelVertexNb " );
        }
//         extractConnMatrix[vertexCount] = AllMeshconnMatrixToAllMesh.getRow(i);
//         (** seedVertexIndex)[vertexCount]=i;
//         vertexCount++;
      }
    }
    if (verbose) std::cout << "Done."<<std::endl;
    return extractConn_ptr;
  }//connMatrixRegionExtract


  Connectivities * connMatrixReducedFromRegion(
    Connectivities * AllMeshconnMatrixToAllMesh_ptr,
    const TimeTexture<short> & seedRegionsTex, int seedRegionLabel,
    std::size_t seedRegionLabelVertexNb,
    std::vector<std::size_t> ** seedVertexIndex, bool verbose )
/*
input: AllMeshconnMatrixToAllMesh_ptr: connectivity matrix of shape (meshVertexNb, other_meshVertex_nb)
output: extractConn_ptr: a connectivity matrix of shape (seedRegionVertexNb, other_meshVertexNb) i.e. the connectivity matrix of the vertex of a given seed region 
!!!: erase the AllMeshconnMatrixToAllMesh_ptr : set all its elements to zero!!!
*/
//   Connectivities * connMatrixReducedFromRegion(Connectivities * AllMeshconnMatrixToAllMesh_ptr, const TimeTexture<short> & seedRegionsTex, int seedRegionLabel, std::size_t seedRegionLabelVertexNb)
  {
    if (verbose) std::cout << "Start connMatrixReducedFromRegion" << std::endl;
    Connectivities & AllMeshconnMatrixToAllMesh  = *AllMeshconnMatrixToAllMesh_ptr;
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    std::size_t meshVertexNb = AllMeshconnMatrixToAllMesh[0].size();//number of columns
    std::size_t seedRegionsNb = AllMeshconnMatrixToAllMesh.size();//number of rows
    std::size_t seedRegionsMeshVertexNb  = seedTex0.nItem();
    if (verbose) std::cout << "seedRegionsMeshVertexNb(rows):" << seedRegionsMeshVertexNb << ", meshVertexNb(cols):" << meshVertexNb << std::endl;
    
//     if (meshVertexNb!=seedRegionsNb or seedTex0.nItem()!=meshVertexNb)
//       throw runtime_error( "error in matrix dimensions" );
    if (verbose) std::cout << "seedRegionLabelVertexNb:" << seedRegionLabelVertexNb << std::endl;
    Connectivities * extractConn_ptr = new Connectivities(seedRegionLabelVertexNb, til::SparseVector<double>(meshVertexNb));//default values = zero
    Connectivities & extractConnMatrix  = *extractConn_ptr;
    * seedVertexIndex = new std::vector<size_t>;
    (** seedVertexIndex).resize(seedRegionLabelVertexNb);
    std::size_t vertexCount = 0;
    if (verbose) std::cout << "coucou1" << std::endl;
    int five_count = int(seedRegionsMeshVertexNb / 20.0);
//     for (std::size_t i = 0; i < seedRegionsMeshVertexNb; ++i)//iteration on the rows
    std::size_t rows_count = 0;
    for (Connectivities::iterator i = AllMeshconnMatrixToAllMesh.begin(); i != AllMeshconnMatrixToAllMesh.end(); ++i)
    {
//       std::cout << "i:" << i << std::endl;
      if (rows_count % five_count == 0)
      {
        if (verbose) std::cout << 5 * int(rows_count / five_count) << "%..." << std::flush;
      }
      if (seedTex0[rows_count]==seedRegionLabel)
      {
        if (vertexCount >= seedRegionLabelVertexNb)
        {
          if (verbose) std::cout <<"vertexCount >= seedRegionLabelVertexNb : "<< vertexCount << ">= " << seedRegionLabelVertexNb << std::endl;
          throw runtime_error( "error in seedRegionLabelVertexNb " );
        }
        extractConnMatrix[vertexCount] = AllMeshconnMatrixToAllMesh[rows_count];
        (** seedVertexIndex)[vertexCount]=rows_count;
        vertexCount++;
      }
      for (Connectivity::Map::iterator j = i->getMap().begin(); j != i->getMap().end(); )
      {
        i->getMap().erase(j++);
      }
    rows_count++;
    }
    if (verbose) std::cout << "Done."<<std::endl;
    return extractConn_ptr;
  }//connMatrixReducedFromRegion


  Connectivities * connMatrixRegionExtractTargetsRegroup(
    Connectivities * allMeshConnMatrixToAllMesh_ptr,
    const TimeTexture<short> & seedRegionsTex, int seedRegionLabel,
    const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
    std::size_t seedRegionLabelVertexNb,
    std::vector<std::size_t> ** seedVertexIndex, bool verbose )
  {
    if (verbose) std::cout << "Start connMatrixRegionExtractRegroup" << std::endl;
    /* Write connectivity matrix of the vertex of a given seed region (label = seedRegionLabel)with a grouping according to targetRegions given in targetRegionsTexture (labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to targetRegionsNb = target regions
    input: Connectivities * connMatrixToAllMesh_ptr, shape = (meshVertexNb, meshVertex_nb)
    output: Connectivities * extractConn_ptr, shape = (seedRegionVertexNb, targetRegionsNb)
    */
    Connectivities & allMeshConnMatrixToAllMesh  = *allMeshConnMatrixToAllMesh_ptr;
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    std::size_t meshVertexNb = allMeshConnMatrixToAllMesh[0].size();//number of columns
    std::size_t seedRegionsNb = allMeshConnMatrixToAllMesh.size();//number of rows
    if (verbose) std::cout << "meshVertexNb:" << meshVertexNb << std::endl;
    if (meshVertexNb!=seedRegionsNb or seedTex0.nItem()!=meshVertexNb or targetTex0.nItem()!=meshVertexNb)
      throw runtime_error( "error in matrix dimensions" );
    Connectivities * extractConn_ptr = new Connectivities(seedRegionLabelVertexNb, til::SparseVector<double>(targetRegionsNb));//default values = zero
    Connectivities & extractConnMatrix  = *extractConn_ptr;
    * seedVertexIndex = new std::vector<size_t>;
    (** seedVertexIndex).resize(seedRegionLabelVertexNb);
    std::size_t vertexCount = 0;
    int targetLabel;
    int five_count = int(meshVertexNb / 20.0);
    for (std::size_t i = 0; i < meshVertexNb; ++i)//iteration on the rows
    {
      if (i % five_count == 0)
      {
        if (verbose) std::cout << 5 * int(i / five_count) << "%..." << std::flush;
      }
      if (seedTex0[i]==seedRegionLabel)
      {
        for (std::size_t j = 0; j < meshVertexNb; ++j)
        {
          targetLabel = targetTex0[j];
          if (targetLabel>0 && targetLabel<= targetRegionsNb)
          {
            extractConnMatrix[vertexCount][targetLabel-1] += allMeshConnMatrixToAllMesh[i][j];
          }
        }
        (** seedVertexIndex)[vertexCount]=i;
        vertexCount++;
      }
    }
    if (verbose) std::cout << "Done." << std::endl;
    return extractConn_ptr;
  }//connMatrixRegionExtractTargetsRegroup


  Connectivities * connMatrixSeedMesh_to_targetMeshTargets_regroup(
    const SparseOrDenseMatrix & connMatrixToAllMesh,
    const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
    bool verbose )
  {
    if (verbose)
      cout << "Start connMatrixRegionExtractRegroup" << endl;
    /* Write connectivity matrix of the vertex of the seedMesh with a grouping according to targetRegions (in targetMesh) given in targetRegionsTexture (labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to targetRegionsNb = target regions
    input: Connectivities * , shape = (seedMeshVertexNb, targetMeshVertexNb), generated by connMatrixSeedMeshToTargetMesh Method
    output: Connectivities * connMatrixSeedMeshToTargetMeshTargets_ptr, shape = (seedMeshVertexNb, targetMeshTargetsRegionsNb)
    */
    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    size_t targetMeshVertexNb = connMatrixToAllMesh.getSize2();//number of columns
    size_t seedMeshVertexNb = connMatrixToAllMesh.getSize1();//number of rows
    if (verbose) cout << "("<< seedMeshVertexNb << "," << targetMeshVertexNb << ") : " << flush;
    if (targetTex0.nItem()!=targetMeshVertexNb)
      throw runtime_error( "error in matrix dimensions" );
    Connectivities * connMatrixSeedMeshToTargetMeshTargets_ptr = new Connectivities(seedMeshVertexNb, til::SparseVector<double>(targetRegionsNb));//default values = zero
    Connectivities & connMatrixSeedMeshToTargetMeshTargets  = *connMatrixSeedMeshToTargetMeshTargets_ptr;
    int targetLabel;
    int five_count = int(seedMeshVertexNb / 20.0);
    for( size_t i = 0; i < seedMeshVertexNb; ++i )//iteration on the rows
    {
      for(size_t j = 0; j < targetMeshVertexNb; ++j)
      {
        targetLabel = targetTex0[j];
        if (targetLabel>0 && targetLabel<= targetRegionsNb)
        {
          connMatrixSeedMeshToTargetMeshTargets[i][targetLabel-1] += connMatrixToAllMesh(i,j);
        }
      }
    }
    if (verbose)
      cout << "Done." << endl;
    return connMatrixSeedMeshToTargetMeshTargets_ptr;
  }//connMatrixSeedMesh_to_targetMeshTargets_regroup


  Connectivities * connMatrixSeedMeshRegions_to_targetMeshTargets_regroup(
    Connectivities * connMatrixSeedMeshToTargetMesh_ptr,
    const TimeTexture<short> & seedRegionsTex,
    const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
    int seedRegionsNb, bool verbose )
  {
    /* Write connectivity matrix of regions of the seedMesh with a grouping according to targetRegions (in targetMesh) given in targetRegionsTexture (labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to targetRegionsNb = target regions) (and a groupin according to seedRegions (in seedMesh) given in seedRegionsTexture)
    input: Connectivities * connMatrixSeedMeshToTargetMesh_ptr, shape = (seedMeshVertexNb, targetMeshVertexNb), generated by connMatrixSeedMeshToTargetMesh Method
    output: Connectivities * connMatrixSeedMeshRegionsToTargetMeshTargets_ptr, shape = (seedMeshVertexNb, targetMeshTargetsRegionsNb)
    */
    Connectivities & connMatrixSeedMeshToTargetMesh  = *connMatrixSeedMeshToTargetMesh_ptr;
    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    std::size_t targetMeshVertexNb = connMatrixSeedMeshToTargetMesh[0].size();//number of columns
    std::size_t seedMeshVertexNb = connMatrixSeedMeshToTargetMesh.size();//number of rows
    if (verbose) std::cout << "seedMeshVertexNb:" << seedMeshVertexNb << std::endl;
    if (targetTex0.nItem()!=targetMeshVertexNb or seedTex0.nItem()!=seedMeshVertexNb)
      throw runtime_error( "error in matrix dimensions" );
    Connectivities * connMatrixSeedMeshRegionsToTargetMeshTargets_ptr = new Connectivities(seedRegionsNb, til::SparseVector<double>(targetRegionsNb));//default values = zero
    Connectivities & connMatrixSeedMeshRegionsToTargetMeshTargets  = *connMatrixSeedMeshRegionsToTargetMeshTargets_ptr;
    int targetLabel;
    int seedRegionLabel;
    int five_count = int(seedMeshVertexNb / 20.0);
    for (std::size_t i = 0; i < seedMeshVertexNb; ++i)//iteration on the rows
    {
      if (i % five_count == 0)
      {
        if (verbose) std::cout << 5 * int(i / five_count) << "%..." << std::flush;
      }
      seedRegionLabel = seedTex0[i];
      if (seedRegionLabel>0 && seedRegionLabel<= seedRegionsNb)
      {
        for (std::size_t j = 0; j < targetMeshVertexNb; ++j)
        {
          targetLabel = targetTex0[j];
          if (targetLabel>0 && targetLabel<= targetRegionsNb)
          {
            connMatrixSeedMeshRegionsToTargetMeshTargets[seedRegionLabel-1][targetLabel-1] += connMatrixSeedMeshToTargetMesh[i][j];
          }
        }
      }
    }
    if (verbose) std::cout << "Done." << std::endl;
    return connMatrixSeedMeshRegionsToTargetMeshTargets_ptr;
  }// connMatrixSeedMeshRegions_to_targetMeshTargets_regroup
  
  void writeAimsFmtConnMatrix(Connectivities * connMatrix_ptr,std::string file_name, bool verbose)
  {
    //Write connectivity matrix (Aims format)
    if (verbose) std::cout << "Writing connectivity matrix ";
    typedef AimsData< float > Matrix;

    Connectivities & connMatrix  = *connMatrix_ptr;
    std::size_t colNb = connMatrix[0].size();
    std::size_t rowsNb = connMatrix.size();
    if (verbose) std::cout << "("<< rowsNb << ", " << colNb <<") in Aims Format...";
    Matrix matrix(rowsNb, colNb, 1);
    matrix.setSizeXYZT(100.0/rowsNb , 80.0/colNb);
    for (std::size_t i = 0; i < rowsNb; ++i)
    {
      for (std::size_t j = 0; j < colNb; ++j)
      {
        matrix(i, j, 0) = connMatrix[i][j];
      }
    }
    Writer<Matrix> w( file_name );
    w.write( matrix );
    if (verbose) std::cout << "done." << std::endl;
  }//writeAimsFmtConnMatrix


  void connMatrixNormalize( SparseOrDenseMatrix & inConnMatrix,
                            bool verbose )
  {
    if (verbose) std::cout << "Normalize connectivity matrix " << std::flush;
    std::size_t colNb = inConnMatrix.getSize2();
    std::size_t rowsNb = inConnMatrix.getSize1();
//     Connectivity sparceVectorIdentity = Connectivity(colNb);
//     for(std::size_t j = 0; j < colNb; ++j) sparceVectorIdentity[j]=1.;
    if (verbose) std::cout << "("<<rowsNb << "," << colNb << ") : " << std::flush;
//    int five_count = int(rowsNb / 20.0);
    for (std::size_t i = 0; i < rowsNb; ++i)
    {
//       if (i % five_count == 0)
//       {
//         if (five_count != 0)
//         {
//           std::cout << 5 * int(i / five_count) << "%..." << std::flush;
//         }
//       }
      vector<double> line_i = inConnMatrix.getSparseRow<vector<double> >(i);
      double norm_i = 0;
      for( std::size_t j = 0; j < colNb; ++j )
        norm_i += line_i[j] * line_i[j];
      norm_i = sqrt( norm_i );
//         std::cout << ", norm_i:" << norm_i <<", sum_i:" << sum_i <<std::endl;
      if(norm_i!=0)
      {
        for(std::size_t j = 0; j < colNb; ++j)
        {
          float connval = inConnMatrix( i, j );
          if( connval != 0 )
            inConnMatrix.set_element( i, j, connval / norm_i );
        }
      }
    }
    if (verbose) std::cout << "Done." << std::endl;
  }//connMatrixNormalize
  
  TimeTexture<float> densityTexture(Connectivities * connMatrixToAllMesh_ptr, std::vector<std::size_t> vertexIndex, bool verbose)
  {
    /*
    input:    connMatrixToAllMesh : connectivity matrix of shape (p, n), p = seed region vertex nb, n = mesh vertex nb
              vertexIndex : index of the seed region vertex (matrix rows) in the main mesh ( vertexIndex[i]=index corresponding to the ith row of the connectivity matrix in the mesh)

    output:   outputDensityTex : texture of connection density : for each vertex of the seed region : sum of its connections to all the other mesh vertex
    */
    if (verbose) std::cout << "Computing density texture" << std::flush;
    
    Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
    std::size_t colNb = connMatrixToAllMesh[0].size();
    std::size_t rowsNb = connMatrixToAllMesh.size();
    if (verbose) std::cout << "("<<rowsNb << "," << colNb << ") : " << std::flush;
    if(rowsNb != vertexIndex.size())
      throw runtime_error( "Not the same dimensions as vertexIndex vector" );
    Connectivity sparseVectorIdentity = Connectivity(colNb);
    for(std::size_t j = 0; j < colNb; ++j) sparseVectorIdentity[j]=1.;
    int five_count = int(rowsNb / 20.0);
    float currentSum;
    int currentVertexIndex;
    TimeTexture<float> outputDensityTex; //(1, til::size(res));
    int countVertex = 0;
    outputDensityTex[0].reserve(colNb);
    for (std::size_t v = 0; v < colNb; ++v) outputDensityTex[0].push_back(-1);
    for (std::size_t i = 0; i < rowsNb; ++i)
    {
      currentVertexIndex = vertexIndex[i];
      if (i % five_count == 0)
      {
        if (verbose) std::cout << 5 * int(i / five_count) << "%..." << std::flush;
      }
      outputDensityTex[0][currentVertexIndex]=0;
      Connectivity line_i = connMatrixToAllMesh[i];
      if (line_i.is_null()!=1)
      {
        currentSum = til::dot<float>(line_i, sparseVectorIdentity);
        if (currentSum==-1 and verbose) std::cout<<"-1,"<<std::flush;
        if (isnan(currentSum) and verbose) std::cout << "currentSum is nan!!" <<std::flush;
        if (currentSum>0)
        {
          outputDensityTex[0][currentVertexIndex]=currentSum;
        }
      }
      countVertex++;
    }
    if (verbose) std::cout << "countVertex++:" << countVertex++ << ", Done." << std::endl;
    return outputDensityTex;
  }


  TimeTexture<float> densityTexture( const SparseOrDenseMatrix & connMatrixToAllMesh, std::vector<std::size_t> vertexIndex, bool verbose)
  {
    /*
    input:    connMatrixToAllMesh : connectivity matrix of shape (p, n), p = seed region vertex nb, n = mesh vertex nb
              vertexIndex : index of the seed region vertex (matrix rows) in the main mesh ( vertexIndex[i]=index corresponding to the ith row of the connectivity matrix in the mesh)

    output:   outputDensityTex : texture of connection density : for each vertex of the seed region : sum of its connections to all the other mesh vertex
    */
    if (verbose) std::cout << "Computing density texture" << std::flush;

    std::size_t colNb = connMatrixToAllMesh.getSize2();
    std::size_t rowsNb = connMatrixToAllMesh.getSize1();
    if (verbose) std::cout << "("<<rowsNb << "," << colNb << ") : " << std::flush;
    if(rowsNb != vertexIndex.size())
      throw runtime_error( "Not the same dimensions as vertexIndex vector" );
    int five_count = int(rowsNb / 20.0);
    float currentSum;
    int currentVertexIndex;
    TimeTexture<float> outputDensityTex; //(1, til::size(res));
    int countVertex = 0;
    outputDensityTex[0].reserve(colNb);
    for (std::size_t v = 0; v < colNb; ++v)
      outputDensityTex[0].push_back(-1);
    for (std::size_t i = 0; i < rowsNb; ++i)
    {
      currentVertexIndex = vertexIndex[i];
      if (i % five_count == 0)
      {
        if (verbose) std::cout << 5 * int(i / five_count) << "%..." << std::flush;
      }
      outputDensityTex[0][currentVertexIndex]=0;
      vector<double> line_i = connMatrixToAllMesh.getRow(i);
      for( int j=0; j<colNb; ++j )
        currentSum += line_i[j];
      if (currentSum==-1 and verbose) std::cout<<"-1,"<<std::flush;
      if (isnan(currentSum) and verbose) std::cout << "currentSum is nan!!" <<std::flush;
      if (currentSum>0)
      {
        outputDensityTex[0][currentVertexIndex] = (float) currentSum;
      }
      countVertex++;
    }
    if (verbose) std::cout << "countVertex++:" << countVertex++ << ", Done." << std::endl;
    return outputDensityTex;
  }


  TimeTexture<float> meshDensityTexture(
    Connectivities * connMatrixToAllMesh_ptr, bool verbose )
  {
  /*
    input:    connMatrixToAllMesh : connectivity matrix of shape (p, n), p = seed region vertex nb, n = mesh vertex nb

    output:   outputDensityTex : texture of connection density : for each vertex of the mesh : sum of its connections to all seed region vertex
  */
    if (verbose) std::cout << "Computing Mesh density texture" << std::flush;
  
    Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
    std::size_t colNb = connMatrixToAllMesh[0].size();
    std::size_t rowsNb = connMatrixToAllMesh.size();
    if (verbose) std::cout << "("<<rowsNb << "," << colNb << ") : " << std::flush;
    TimeTexture<float> outputDensityTex; //(1, til::size(res));
    outputDensityTex[0].reserve(colNb);
    Connectivity * totalConnectionDensity_ptr = connMatrixSumRows(connMatrixToAllMesh_ptr);
    Connectivity & totalConnectionDensity = *totalConnectionDensity_ptr;
    for (std::size_t v = 0; v < colNb; ++v) outputDensityTex[0].push_back(-1);
    for (std::size_t j = 0; j < colNb; ++j)
    {
      outputDensityTex[0][j]= totalConnectionDensity[j];
//       if (outputDensityTex[0][j]!= 0)
//       {
//         std::cout << "outputDensityTex[0][j]:, j:"<< j << ", " << outputDensityTex[0][j] << std::endl;
//       }
    }
    if (verbose) std::cout << "Done." << std::endl;
    delete totalConnectionDensity_ptr;
    return outputDensityTex;
  }//meshDensityTexture


  TimeTexture<float> meshDensityTexture(
    const SparseOrDenseMatrix & connMatrixToAllMesh, bool verbose )
  {
  /*
    input:    connMatrixToAllMesh : connectivity matrix of shape (p, n), p = seed region vertex nb, n = mesh vertex nb

    output:   outputDensityTex : texture of connection density : for each vertex of the mesh : sum of its connections to all seed region vertex
  */
    if (verbose) std::cout << "Computing Mesh density texture" << std::flush;

    std::size_t colNb = connMatrixToAllMesh.getSize2();
    std::size_t rowsNb = connMatrixToAllMesh.getSize1();
    if (verbose) std::cout << "("<<rowsNb << "," << colNb << ") : " << std::flush;
    TimeTexture<float> outputDensityTex; //(1, til::size(res));
    outputDensityTex[0].reserve(colNb);
    vector<double> * totalConnectionDensity_ptr = connMatrixSumRows(connMatrixToAllMesh);
    vector<double> & totalConnectionDensity = *totalConnectionDensity_ptr;
    for (std::size_t v = 0; v < colNb; ++v)
      outputDensityTex[0].push_back(-1);
    for (std::size_t j = 0; j < colNb; ++j)
    {
      outputDensityTex[0][j]= totalConnectionDensity[j];
//       if (outputDensityTex[0][j]!= 0)
//       {
//         std::cout << "outputDensityTex[0][j]:, j:"<< j << ", " << outputDensityTex[0][j] << std::endl;
//       }
    }
    if (verbose) std::cout << "Done." << std::endl;
    delete totalConnectionDensity_ptr;
    return outputDensityTex;
  }//meshDensityTexture


  template <typename T>
  TimeTexture<float> * oneTargetDensityTargetsRegroupTexture_T(
    const T * lineMatrixToTargetRegions_ptr,
    const TimeTexture<short> & targetRegionsTex, int timestep )
  {
    /*
    Compute density texture (associated to a mesh) of a targetRegion towards all the target regions
    inputs:     lineMatrixToTargetsRegions_ptr : connectivity of one target to all the target regions, shape = (1,targetRegionsNb)
                targetRegionsTex : labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to targetRegionsNb = target regions
    output:     TimeTexture<float> * outputDensityTex : pointer on output density texture
    */
    TimeTexture<short>::const_iterator it = targetRegionsTex.find( timestep );
    if( it == targetRegionsTex.end() )
      it = targetRegionsTex.begin();
    const Texture<short> & targetTex0 = it->second;
    std::size_t meshVertexNb = targetTex0.nItem();
    const T & lineMatrixToTargetRegions = *lineMatrixToTargetRegions_ptr;
    std::size_t targetRegionsNb = lineMatrixToTargetRegions.size();
    TimeTexture<float> * outputDensityTex_ptr = new TimeTexture<float>;
    TimeTexture<float> & outputDensityTex = *outputDensityTex_ptr;
    outputDensityTex[0].reserve(meshVertexNb);
    for (std::size_t v = 0; v < meshVertexNb; ++v) outputDensityTex[0].push_back(-1);
    
    for (std::size_t i = 0; i < meshVertexNb; i++)
    {
      int label = targetTex0[i];//equivalent to targetRegionsTex[0].item(i)
      if (1 <= label && label <= targetRegionsNb)
      {
        outputDensityTex[0][i] = (float) lineMatrixToTargetRegions[label-1];
      }
      else outputDensityTex[0][i] = -1;
    }
    return outputDensityTex_ptr;
  }//oneTargetDensityTargetsRegroupTexture


  TimeTexture<float> * oneTargetDensityTargetsRegroupTexture(
    const Connectivity * lineMatrixToTargetRegions_ptr,
    const TimeTexture<short> & targetRegionsTex, int timestep )
  {
    return oneTargetDensityTargetsRegroupTexture_T(
      lineMatrixToTargetRegions_ptr, targetRegionsTex, timestep );
  }


  TimeTexture<float> * oneTargetDensityTargetsRegroupTexture(
    const vector<double> * lineMatrixToTargetRegions_ptr,
    const TimeTexture<short> & targetRegionsTex, int timestep )
  {
    return oneTargetDensityTargetsRegroupTexture_T(
      lineMatrixToTargetRegions_ptr, targetRegionsTex, timestep );
  }


  SparseOrDenseMatrix* connectivitiesToSparseOrDenseMatrix(
    const Connectivities & conn )
  {
    unsigned ncol = 0;
    if( !conn.empty() )
      ncol = conn.begin()->size();

    SparseOrDenseMatrix* mat = new SparseOrDenseMatrix( conn.size(), ncol );

    unsigned line, nline = conn.size();
    unsigned long count = 0;
    for( line=0; line<nline; ++line )
      count += conn[line].getMap().size();

    if( count >= mat->optimalShapeThreshold() )
      mat->muteToDense();

    Connectivity::sparse_const_iterator ic, ec;

    for( line=0; line<nline; ++line )
    {
      const Connectivity & cline = conn[line];
      if( !cline.empty() )
      {
        for( ic=cline.sparse_begin(), ec=cline.sparse_end(); ic!=ec; ++ic )
          mat->set_element( line, ic->first, ic->second );
      }
    }

    return mat;
  }


} // namespace constel
