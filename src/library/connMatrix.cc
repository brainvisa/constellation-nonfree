#include <constellation/connMatrix.h>
#include <aims/sparsematrix/sparseordensematrix.h>
#include <aims/io/writer.h>
#include <aims/distancemap/meshdistance.h>
#include <aims/utility/stl_conversion.h>


using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;
using namespace boost;

namespace constel {

  //---------------------
  //  weightedConnMatrix
  //---------------------
  Connectivities * weightedConnMatrix(
      const WeightedFibers &fibers, const AimsSurfaceTriangle &inAimsMesh,
      Motion motion, bool verbose)
  {

    size_t nv = inAimsMesh.vertex().size();
    if (verbose) cout << "Number of fibers: " << fibers.size() << endl;

    // Generating kdtree
    KDTreeVertices m = kdt_vertices( inAimsMesh );
    KDTree kdt( m.begin(), m.end() );

    // Computing connectivity matrix
    if (verbose) cout << "Computing connectivity matrix" << endl;
    Connectivities *conn_ptr = new Connectivities( nv, Connectivity( nv ) );
    Connectivities &conn  = *conn_ptr;
    size_t fiberCount = 0;
    size_t nFibers = fibers.size();
    Point3df p1, p2;

    for (WeightedFibers::const_iterator iFiber = fibers.begin();
         iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      // Get fiber informations
      Fiber fiber = iFiber->first;
      double weight = iFiber->second;

      // Transform first point and last point of fiber from t2 to anat space
      // (with motion)

      p1 = motion.transform(fiber.front()[0],
                            fiber.front()[1],
                            fiber.front()[2]);
      p2 = motion.transform(fiber.back()[0],
                            fiber.back()[1],
                            fiber.back()[2]);

      size_t A = kdt.find_nearest( make_pair( 0U, p1 ) ).first->first;
      size_t B = kdt.find_nearest( make_pair( 0U, p2 ) ).first->first;

      // Filling the connectivity matrix with the fiber's weight
      conn[A][B] += weight;
      conn[B][A] += weight;
    }

    return (conn_ptr);
  }

  //--------------
  //  connMatrix
  //--------------
  Connectivities * connMatrix(
      const Fibers &fibers, const AimsSurfaceTriangle &inAimsMesh,
      float distthresh, float wthresh, Motion motion, bool verbose)
  {

    double two_pi = 2*3.1415926535897931;
    const double G_THRESH = 0.001; //threshold for connectivity contributions

    // Comuting geomap : neighborhood map
    if (verbose) cout << "Computing geomap..." << flush;
    if (verbose) cout << "step 1..." << flush;
    size_t nv = inAimsMesh.vertex().size();

    std::vector< std::map<size_t, float> > res;

    if (distthresh != 0)
      meshdistance::pairwiseDistanceMaps( inAimsMesh.begin()->second, res,
                                          distthresh );

    if (verbose) cout << "Number of fibers: " << fibers.size() << endl;

    // Generating kdtree
    KDTreeVertices m = kdt_vertices( inAimsMesh );
    KDTree kdt( m.begin(), m.end() );

    // Computing connectivity matrix
    if (verbose) cout << "Computing connectivity matrix" << endl;
    Connectivities *conn_ptr = new Connectivities( nv, Connectivity( nv ) );
    Connectivities &conn  = *conn_ptr;
    int countRemoved = 0;
    size_t fiberCount = 0;
    size_t nFibers = fibers.size();
    size_t five_count = nFibers / 20.0;
    double total_conn = 0;
    Point3df p1, p2;
    for (Fibers::const_iterator iFiber = fibers.begin();
         iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (five_count != 0)
      {
        if (fiberCount % five_count == 0)
        {
          if (verbose)
            cout << 5 * int(fiberCount / five_count) << "%..." << flush;
        }
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
      if (distthresh == 0)
      {
        conn[A][B] += 1.;
        conn[B][A] += 1.;
      }
      else
      {
        double wtotal = 0;
        list<pair<pair<size_t, size_t>, double> > weights;
        for (map<size_t, float>::const_iterator B_neigh = res[B].begin();
             B_neigh != res[B].end(); ++B_neigh)
        {
          for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
               A_neigh != res[A].end(); ++A_neigh)
          {
            double square_sigma = distthresh*distthresh;
            double e1 = square(A_neigh->second);
            double e2 = square(B_neigh->second);
            // "smoothing coefficient"
            double w1 = exp(-e1 / (2*square_sigma))/(square_sigma * two_pi);
            double w2 = exp(-e2 / (2*square_sigma))/(square_sigma * two_pi);
            if (w1 < G_THRESH) w1 = 0;
            if (w2 < G_THRESH) w2 = 0;
            double w = w1 + w2;
            if (w > 0)
            {
              wtotal += w;
              weights.push_back(
                  make_pair(make_pair(A_neigh->first, B_neigh->first), w));
            }
          }
        }

        list<pair<pair<size_t, size_t>, double> >::iterator iw, ew
          = weights.end();
        for (iw=weights.begin(); iw!=ew; ++iw)
        {
          double w = iw->second / wtotal;
          conn[iw->first.second][iw->first.first] += w;
          conn[iw->first.first][iw->first.second] += w;
          total_conn += 2 * w;
        }
      }
    }
    if (verbose)
    {
      cout << "...100%..." << flush;
      cout << "NB: " << countRemoved << " fibers out of " << fibers.size()
        << " have been discarded" << endl;
    }

    // To compress a little bit, remove points that are below some threshold
    // Don't do this is wthresh is 0
    if ( wthresh > 0 )
    {
      double nonzero_count = 0;
      if (verbose) cout << "Removing weak connections: " << flush;
      int count1 = 0;
      int count2 = 0;
      for (Connectivities::iterator i = conn.begin(); i != conn.end(); ++i)
      {
        for (Connectivity::iterator j = i->begin();
             j != i->end();)
        {
          ++count1;
          if (j->second <= wthresh)
          {
            ++count2;
            i->erase(j++);
          }
          else
          {
            ++j;
            ++nonzero_count;
          }
        }
      }
      if (verbose)
        cout << "Removed " << count2 << " elements out of " << count1 << endl;
    }
    return (conn_ptr);
  }

  //----------------------------------
  //  connMatrixSeedMeshToTargetMesh
  //----------------------------------
  /*
    input:
      distthresh: distance for the smoothing of the matrix
      wthresh: weight threshold for clean the connectivity matrix
               (values below are removed)
      motion: transformation from t2 (fibers) to anat(meshes)

   output:
      conn: sparse matrix of shape = [seedMeshVertex_nb,targetMeshVertex_nb]
    */
  Connectivities * connMatrixSeedMeshToTargetMesh(
    const Fibers & fibers, const AimsSurfaceTriangle & aimsSeedMesh,
    const AimsSurfaceTriangle & aimsTargetMesh, float distthresh,
    float wthresh, Motion motion, bool verbose)
  {

    const double G_THRESH = 0.001; //threshold for connectivity contributions

    // Computing geomap : neighborhood map

    // For seedMesh
    if (verbose) cout << "Computing geomap seedMesh..." << flush;
    if (verbose) cout << "step 1..." << flush;

    size_t nv = aimsSeedMesh.vertex().size();

    std::vector< std::map<size_t, float> > res_seedMesh;

    if (distthresh != 0)
      meshdistance::pairwiseDistanceMaps( aimsSeedMesh.begin()->second,
                                          res_seedMesh, distthresh );

    //For targetMesh
    if (verbose) cout << "Computing geomap targetMesh..." << flush;
    if (verbose) cout << "step 1..." << flush;

    size_t nvt = aimsTargetMesh.vertex().size();

    std::vector< std::map<size_t, float> > res_targetMesh;

    if (distthresh != 0)
      meshdistance::pairwiseDistanceMaps( aimsTargetMesh.begin()->second,
                                          res_targetMesh, distthresh );

    if (verbose) cout << "Number of fibers: " << fibers.size() << endl;

    // Generating kdtree
    // For seedMesh
    if (verbose) cout << "Generating kdtrees" << endl;
    KDTreeVertices m = kdt_vertices( aimsSeedMesh );
    KDTree kdt_seedMesh( m.begin(), m.end() );
    // For targetMesh
    KDTreeVertices m2 = kdt_vertices( aimsTargetMesh );
    KDTree kdt_targetMesh( m2.begin(), m2.end() );

    // Computing connectivity matrix
    if (verbose) cout << "Computing connectivity matrix" << endl;
    Connectivities * conn_ptr = new Connectivities(
        aimsSeedMesh.vertex().size(),
        Connectivity(aimsTargetMesh.vertex().size()));
    Connectivities & conn  = *conn_ptr;
    int countRemoved = 0;
    size_t fiberCount = 0;
    size_t nFibers = fibers.size();
    size_t five_count = nFibers / 20.0;
    double total_conn = 0;
    Point3df p1, p2;
    for (Fibers::const_iterator iFiber = fibers.begin();
         iFiber != fibers.end(); ++iFiber, ++fiberCount) {
      if (five_count != 0) {
        if (fiberCount % five_count == 0) {
          if (verbose)
            cout << 5 * int(fiberCount / five_count) << "%..." << flush;
        }
      }
      //Looking for closest points
      // Transform iFiber->front() and iFiber->back()
      // from t2 to anat space (with motion)
      p1 = motion.transform(iFiber->front()[0],
                            iFiber->front()[1],
                            iFiber->front()[2]);
      p2 = motion.transform(iFiber->back()[0],
                            iFiber->back()[1],
                            iFiber->back()[2]);
      size_t A_seedMesh = kdt_seedMesh.find_nearest(
        make_pair( 0U, p1 ) ).first->first;
      size_t B_seedMesh = kdt_seedMesh.find_nearest(
        make_pair( 0U, p2 ) ).first->first;
      size_t A_targetMesh = kdt_targetMesh.find_nearest(
        make_pair( 0U, p1 ) ).first->first;
      size_t B_targetMesh = kdt_targetMesh.find_nearest(
        make_pair( 0U, p2 ) ).first->first;

      // Filling the connectivity matrix
      /*
      Two retained cases:
        first case: p1 "touches" seedMesh and p2 "touches" targetMesh
        second case: the opposite (p1 "touches" targetMesh and p2 "touches"
                     seedMesh)
      */
      if( dist2( p1, aimsSeedMesh.vertex()[A_seedMesh] ) <= 25.0)
      { //p1 "touches" seedMesh
        if( dist2( p2, aimsTargetMesh.vertex()[B_targetMesh] ) <= 25.0)
        {
          //first case (p2 "touches" targetMesh)
          if (distthresh == 0)
          {
            conn[A_seedMesh][B_targetMesh] += 1.;
          }
          else
          {
            for (map<size_t, float>::const_iterator A_seedMesh_neigh
                 = res_seedMesh[A_seedMesh].begin();
                 A_seedMesh_neigh != res_seedMesh[A_seedMesh].end();
                 ++A_seedMesh_neigh)
            {
              for (map<size_t, float>::const_iterator B_targetMesh_neigh
                  = res_targetMesh[B_targetMesh].begin();
                  B_targetMesh_neigh != res_targetMesh[B_targetMesh].end();
                  ++B_targetMesh_neigh)
              {
                double e1 = square(A_seedMesh_neigh->second);
                double e2 = square(B_targetMesh_neigh->second);
                // "smoothing coefficient"
                double w = exp( -(e1 + e2)  / ( 2*distthresh*distthresh ));
                if (w > G_THRESH)
                {
                  conn[A_seedMesh_neigh->first][B_targetMesh_neigh->first]
                    += w;
                  total_conn += w;
                }
              }
            }
          }
        }
        else
        { // p2na doesn t touch targetMesh
          ++countRemoved;
          continue;
        }
      }
      else
      {
        if( dist2( p1, aimsTargetMesh.vertex()[A_targetMesh] ) <= 25.0)
        {//p1 "touches" targetMesh
          if( dist2( p2, aimsSeedMesh.vertex()[B_seedMesh] ) <= 25.0)
          {
            //second case (p2 "touches" seedMesh)
            if (distthresh == 0)
            {
              conn[B_seedMesh][A_targetMesh] += 1.;
            }
            else
            {
              for (map<size_t, float>::const_iterator A_targetMesh_neigh
                   = res_targetMesh[A_targetMesh].begin();
                   A_targetMesh_neigh != res_targetMesh[A_targetMesh].end();
                   ++A_targetMesh_neigh)
              {
                for (map<size_t, float>::const_iterator B_seedMesh_neigh
                     = res_seedMesh[B_seedMesh].begin();
                     B_seedMesh_neigh != res_seedMesh[B_seedMesh].end();
                     ++B_seedMesh_neigh)
                {
                  double e1 = square(A_targetMesh_neigh->second);
                  double e2 = square(B_seedMesh_neigh->second);
                  // "smoothing coefficient"
                  double w = exp( -(e1 + e2)  / ( 2*distthresh*distthresh ));
                  if (w > G_THRESH)
                  {
                    conn[B_seedMesh_neigh->first][A_targetMesh_neigh->first]
                      += w;
                    total_conn += w;
                  }
                }
              }
            }
          } else { // p2na doesn t touch seedMesh
            ++countRemoved;
            continue;
          }
        } else {// p1 touches nor seedMesh neither targetMesh
          ++countRemoved;
          continue;
        }
      }
    }

    if (verbose) {
      cout << "...100%..." << flush;
      cout << "Done." << endl;
      cout << "NB: " << countRemoved << " fibers out of " << fibers.size()
        << " have been discarded" << endl;
    }
    // To compress a little bit, remove points that are below some threshold
    double nonzero_count = 0;
    if (verbose) cout << "Removing weak connections: " << flush;
    int count1 = 0;
    int count2 = 0;
    for (Connectivities::iterator i = conn.begin(); i != conn.end(); ++i) {
      for (Connectivity::iterator j = i->begin();
           j != i->end(); )
      {
        ++count1;
        if (j->second <= wthresh)
        {
          ++count2;
          i->erase(j++);
        }
        else
        {
          ++j;
          ++nonzero_count;
        }
      }
    }
    if (verbose)
      cout << "Removed " << count2 << " elements out of " << count1 << endl;
    return (conn_ptr);
  }


  //--------------------
  //  connMatrixToRois
  //--------------------
  Connectivity * connMatrixToRois(
    const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh,
    float distthresh, float /* wthresh */, Motion motion, bool verbose )
  {
    if (verbose) cout << "Computing connectivity Matrix to Rois" << endl;
    /*
    Coumputing connectivity matrix between the mesh vertex and ROIs
    (subcortical structures for example): for the moment: on suppose que les
    fibres ont ete filtres en fonction d une ROI et pour chaque element
    du maillage, on compte seulement le nombre de fibres qu'il y a...
    */
    const double G_THRESH = 0.001; //threshold for connectivity contributions

    // Comuting geomap : neighborhood map
    if (verbose) cout << "Computing geomap..." << flush;

    size_t nv = inAimsMesh.vertex().size();
    std::vector< std::map<size_t, float> > res;

    if (distthresh != 0)
      meshdistance::pairwiseDistanceMaps( inAimsMesh.begin()->second, res,
                                          distthresh );

    if (verbose) cout << "Number of fibers: " << fibers.size() << endl;

    // Generating kdtree
    KDTreeVertices m = kdt_vertices( inAimsMesh );
    KDTree kdt( m.begin(), m.end() );

    // Computing connectivity matrix
    if (verbose) cout << "Computing connectivity matrix" << endl;
    Connectivity * conn_ptr = new Connectivity( inAimsMesh.vertex().size() );
    Connectivity & conn  = *conn_ptr;
    int countRemoved = 0;
    size_t fiberCount = 0, nFibers = fibers.size();
    size_t five_count = nFibers / 20.0;
    Point3df p1, p2;
    for (Fibers::const_iterator iFiber = fibers.begin();
         iFiber != fibers.end(); ++iFiber, ++fiberCount)
    {
      if (fiberCount % five_count == 0)
      {
        cout << 5 * int(fiberCount / five_count) << "%..." << flush;
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

      if( dist2( p1, inAimsMesh.vertex()[A] ) > 25.0 &&
          dist2( p2, inAimsMesh.vertex()[B] ) > 25.0)
      {
        ++countRemoved;
        continue;
      }
      // Filling the connectivity matrix
      if( dist2( p1, inAimsMesh.vertex()[A] ) <= 25.0)
      {
        if (distthresh == 0)
        {
          conn[A] += 1.;
        }
        else
        {
          for (map<size_t, float>::const_iterator A_neigh = res[A].begin();
               A_neigh != res[A].end(); ++A_neigh)
          {
            double e = square(A_neigh->second);
            // "smoothing coefficient"
            double w = exp(-e  / (2*distthresh*distthresh));
            if (w > G_THRESH)
            {
              conn[A_neigh->first] += w;
            }
          }
        }
      }
      if( dist2( p2, inAimsMesh.vertex()[B] ) <= 25.0)
      {
        if(distthresh == 0)
        {
          conn[A] += 1.;
        }
        else
        {
          for (map<size_t, float>::const_iterator B_neigh = res[B].begin();
               B_neigh != res[B].end(); ++B_neigh)
          {
            double e = square(B_neigh->second);
            // smoothing coefficient
            double w = exp(-e  / (2*distthresh*distthresh));
            if (w > G_THRESH)
            {
              conn[B_neigh->first] += w;
            }
          }
        }
      }
    }
    if (verbose)
      cout << endl << "NB: " << countRemoved << " fibers out of "
        << fibers.size() << " have been discarded" << endl;
    return (conn_ptr);
  }

  //--------------------------------------
  //  connMatrixRow_TO_TimeTexture_FLOAT
  //--------------------------------------
  TimeTexture<float> *connMatrixRow_TO_TimeTexture_FLOAT(
      Connectivity *conn_ptr)
  {
    Connectivity &conn = *conn_ptr;
    TimeTexture<float> *connTex_ptr = new TimeTexture<float>;
    TimeTexture<float> &connTex = *connTex_ptr;
    size_t conn_size = conn.size();
    connTex.reserve(conn_size);
    for (size_t i = 0; i < conn_size; i++)
    {
      connTex[0].push_back(conn[i]);
    }
    return connTex_ptr;
  }

  //---------------------
  //  connMatrixSumRows
  //---------------------
  /*
  Compute the sum of the matrix rows
  input: Connectivities matrix_ptr, shape = (r,l)
  output: Connectivity matrixSumRows, shape = (1,l)
  */
  Connectivity *connMatrixSumRows(Connectivities *matrix_ptr, bool verbose)
  {
    Connectivities &matrix = *matrix_ptr;
    size_t colNb = matrix[0].size();
    size_t rowsNb = matrix.size();
    if (verbose) cout << "("<<rowsNb << "," << colNb << ") : " << flush;
    Connectivity sparseVectorIdentity = Connectivity(colNb);
    Connectivity * matrixSumRows_ptr = new Connectivity(colNb);
    Connectivity & matrixSumRows  = *matrixSumRows_ptr;
    for (size_t j = 0; j < colNb; ++j)
    {
      sparseVectorIdentity[j]=1.;
      matrixSumRows[j]=0;
    }
    for (size_t i = 0; i < rowsNb; ++i)
    {
      Connectivity &line_i = matrix[i];
      if( !line_i.empty() )
      {
        matrixSumRows += line_i;
      }
    }
    return (matrixSumRows_ptr);
  }

  //---------------------
  //  connMatrixSumRows
  //---------------------
  /*
  Compute the sum of the matrix rows
  input: Connectivities matrix_ptr, shape = (r,l)
  output: Connectivity matrixSumRows, shape = (1,l)
  */
  vector<double> *connMatrixSumRows(
      const SparseOrDenseMatrix & matrix, bool verbose)
  {
    size_t colNb = matrix.getSize2();
    size_t rowsNb = matrix.getSize1();
    if (verbose) cout << "(" <<rowsNb << "," << colNb << ") : " << flush;
    vector<double> * matrixSumRows_ptr = new vector<double>(colNb);
    vector<double> & matrixSumRows  = *matrixSumRows_ptr;
    for(size_t j = 0; j < colNb; ++j)
    {
      matrixSumRows[j]=0;
    }
    for (size_t i = 0; i < rowsNb; ++i)
    {
      vector<double> line_i = matrix.getRow(i);
      for( size_t j=0; j<colNb; ++j )
      {
        matrixSumRows[j] += line_i[j];
      }
    }
    return (matrixSumRows_ptr);
  }

  //----------------------------
  //  connMatrixTargetsRegroup
  //----------------------------
  /*
  Computing connectivity matrix regrouped according to targetRegions given in
  targetRegionsTexture (labeled texture, between 0 and targetRegionsNb,
  0 = background, 1 to targetRegionsNb = target regions
  input: Connectivities * connMatrixToAllMesh_ptr,
         shape = (meshVertex_nb,meshVertexNb)
  output: Connectivities * regroupConn_ptr,
          shape = (targetRegions Nb, targetRegionsNb)
  */
  Connectivities * connMatrixTargetsRegroup(
      Connectivities * connMatrixToAllMesh_ptr,
      const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
      bool verbose) {
    if (verbose) cout << "Start connMatrixTargetsRegroup" << endl;
    Connectivities & connMatrixToAllMesh  = *connMatrixToAllMesh_ptr;
    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    size_t meshVertexNb = connMatrixToAllMesh[0].size(); //number of columns
    size_t seedRegionsNb = connMatrixToAllMesh.size(); //number of rows

    if (meshVertexNb != seedRegionsNb or targetTex0.nItem() != meshVertexNb) {
      if (verbose) cout << "error in matrix dimensions" << endl;
    }
    Connectivities *regroupConn_ptr = new Connectivities(
        targetRegionsNb, Connectivity(targetRegionsNb));
    Connectivities &regroupConnMatrix  = *regroupConn_ptr;

    if (verbose)
      cout << "Writing connectivity matrix between the target regions: ("
        << targetRegionsNb << ","  << targetRegionsNb << ") matrix."<< endl;
    int ten_count = int(meshVertexNb / 10.0);

    for (size_t i = 0; i < meshVertexNb; i++) {
      int label = targetTex0[i];//equivalent to targetRegionsTex[0].item(i)
      if (1 <= label && label <= targetRegionsNb) {
        if (i % ten_count == 0) {
          if (verbose) cout << 10 * int(i / ten_count) << "%..." << flush;
        }

        for (size_t j = 0; j < i; ++j) {
          int label2 = targetTex0[j];
          if (1 <= label2 && label2 <= targetRegionsNb && i != j) {
            regroupConnMatrix[label-1][label2-1] += connMatrixToAllMesh[i][j];
            regroupConnMatrix[label2-1][label-1] += connMatrixToAllMesh[i][j];
          }
        }
      }
    }
    return regroupConn_ptr;
  }

  //---------------------------
  //  connMatrixRegionExtract
  //---------------------------
  /*
  Write connectivity matrix of the vertex of a given seed region
  (label = seedRegionLabel)
  input: Connectivities * connMatrixToAllMesh_ptr,
         shape = (meshVertexNb, other_meshVertex_nb)
  output: Connectivities * extractConn_ptr,
          shape = (seedRegionVertexNb, other_meshVertexNb)
  */
  Connectivities * connMatrixRegionExtract(
    const SparseOrDenseMatrix & AllMeshconnMatrixToAllMesh,
    const TimeTexture<short> & seedRegionsTex, int seedRegionLabel,
    size_t seedRegionLabelVertexNb,
    vector<size_t> ** seedVertexIndex, bool verbose) {
    if (verbose) cout << "Start connMatrixRegionExtract" << endl;
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    size_t meshVertexNb = AllMeshconnMatrixToAllMesh.getSize2();//nb columns
    //size_t seedRegionsNb = AllMeshconnMatrixToAllMesh.getSize1();//nb of rows
    size_t seedRegionsMeshVertexNb  = seedTex0.nItem();
    if (verbose)
      cout << "seedRegionsMeshVertexNb(rows):" << seedRegionsMeshVertexNb
        << ", meshVertexNb(cols):" << meshVertexNb << endl;

    if (verbose) cout << "seedRegionLabelVertexNb:"
      << seedRegionLabelVertexNb << endl;
    Connectivities * extractConn_ptr = new Connectivities(
        seedRegionLabelVertexNb, Connectivity(meshVertexNb));
    * seedVertexIndex = new vector<size_t>;
    (** seedVertexIndex).resize(seedRegionLabelVertexNb);
    size_t vertexCount = 0;

    int five_count = int(seedRegionsMeshVertexNb / 20.0);
    for (size_t i = 0; i < seedRegionsMeshVertexNb; ++i) {
      if (i % five_count == 0) {
        if (verbose) cout << 5 * int(i / five_count) << "%..." << flush;
      }
      if (seedTex0[i]==seedRegionLabel) {
        if (vertexCount >= seedRegionLabelVertexNb) {
          if(verbose) cout <<"vertexCount >= seedRegionLabelVertexNb : "
            << vertexCount << ">= " << seedRegionLabelVertexNb << endl;
          throw runtime_error( "error in seedRegionLabelVertexNb ");
        }
      }
    }
    return extractConn_ptr;
  }

  //-------------------------------
  //  connMatrixReducedFromRegion
  //-------------------------------
  /*
  input: AllMeshconnMatrixToAllMesh_ptr: connectivity matrix of shape
         (meshVertexNb, other_meshVertex_nb)
  output: extractConn_ptr: a connectivity matrix of shape
          (seedRegionVertexNb, other_meshVertexNb) i.e. the connectivity matrix
          of the vertex of a given seed region
  !: erase the AllMeshconnMatrixToAllMesh_ptr : set all its elements to zero!
  */
  Connectivities * connMatrixReducedFromRegion(
    Connectivities * AllMeshconnMatrixToAllMesh_ptr,
    const TimeTexture<short> & seedRegionsTex, int seedRegionLabel,
    int seedRegionLabelVertexNb,
    vector<size_t> ** seedVertexIndex, bool verbose)
  {
    if (verbose) cout << "Start connMatrixReducedFromRegion" << endl;
    Connectivities &AllMeshconnMatrixToAllMesh
      = *AllMeshconnMatrixToAllMesh_ptr;
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    size_t meshVertexNb = AllMeshconnMatrixToAllMesh[0].size();//nb of columns
    //size_t seedRegionsNb = AllMeshconnMatrixToAllMesh.size();//nb of rows
    size_t seedRegionsMeshVertexNb  = seedTex0.nItem();
    if (verbose)
      cout << "seedRegionsMeshVertexNb(rows):" << seedRegionsMeshVertexNb
        << ", meshVertexNb(cols):" << meshVertexNb << endl;
    if (verbose)
      cout << "seedRegionLabelVertexNb:" << seedRegionLabelVertexNb << endl;
    Connectivities * extractConn_ptr = new Connectivities(
        seedRegionLabelVertexNb, Connectivity(meshVertexNb));
    Connectivities & extractConnMatrix  = *extractConn_ptr;
    * seedVertexIndex = new vector<size_t>;
    (** seedVertexIndex).resize(seedRegionLabelVertexNb);
    size_t vertexCount = 0;
    int five_count = int(seedRegionsMeshVertexNb / 20.0);
    size_t rows_count = 0;
    for (Connectivities::iterator i = AllMeshconnMatrixToAllMesh.begin();
         i != AllMeshconnMatrixToAllMesh.end(); ++i) {
      if (rows_count % five_count == 0) {
        if (verbose)
          cout << 5 * int(rows_count / five_count) << "%..." << flush;
      }
      if (seedTex0[rows_count]==seedRegionLabel) {
        if (int(vertexCount) >= seedRegionLabelVertexNb) {
          if (verbose)
            cout <<"vertexCount >= seedRegionLabelVertexNb : "
              << vertexCount << ">= " << seedRegionLabelVertexNb << endl;
          throw runtime_error("error in seedRegionLabelVertexNb ");
        }
        extractConnMatrix[vertexCount]
          = AllMeshconnMatrixToAllMesh[rows_count];
        (** seedVertexIndex)[vertexCount]=rows_count;
        vertexCount++;
      }
      for (Connectivity::iterator j = i->begin();
           j != i->end();)
      {
        i->erase(j++);
      }
    rows_count++;
    }
    return extractConn_ptr;
  }

  //-----------------------------------------
  //  connMatrixRegionExtractTargetsRegroup
  //-----------------------------------------
  /*
  Write connectivity matrix of the vertex of a given seed region
  (label = seedRegionLabel)with a grouping according to targetRegions given in
  targetRegionsTexture (labeled texture, between 0 and targetRegionsNb,
  0 = background, 1 to targetRegionsNb = target regions
  input: Connectivities *connMatrixToAllMesh_ptr, shape
           = (meshVertexNb, meshVertex_nb)
  output: Connectivities *extractConn_ptr, shape
           = (seedRegionVertexNb, targetRegionsNb)
  */
  Connectivities * connMatrixRegionExtractTargetsRegroup(
      Connectivities * allMeshConnMatrixToAllMesh_ptr,
      const TimeTexture<short> & seedRegionsTex, int seedRegionLabel,
      const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
      size_t seedRegionLabelVertexNb,
      vector<size_t> ** seedVertexIndex, bool verbose) {
    if (verbose) cout << "Start connMatrixRegionExtractRegroup" << endl;

    Connectivities & allMeshConnMatrixToAllMesh
      = *allMeshConnMatrixToAllMesh_ptr;
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    size_t meshVertexNb = allMeshConnMatrixToAllMesh[0].size();//nb of columns
    size_t seedRegionsNb = allMeshConnMatrixToAllMesh.size();//nb of rows
    if (verbose) cout << "meshVertexNb:" << meshVertexNb << endl;
    if (meshVertexNb != seedRegionsNb or seedTex0.nItem() != meshVertexNb
        or targetTex0.nItem() != meshVertexNb) {
      throw runtime_error("error in matrix dimensions");
    }
    Connectivities *extractConn_ptr = new Connectivities(
        seedRegionLabelVertexNb, Connectivity(targetRegionsNb));
    Connectivities & extractConnMatrix  = *extractConn_ptr;
    * seedVertexIndex = new vector<size_t>;
    (** seedVertexIndex).resize(seedRegionLabelVertexNb);
    size_t vertexCount = 0;
    int targetLabel;
    int five_count = int(meshVertexNb / 20.0);
    for (size_t i = 0; i < meshVertexNb; ++i) {
      if (i % five_count == 0) {
        if (verbose) cout << 5 * int(i / five_count) << "%..." << flush;
      }
      if (seedTex0[i]==seedRegionLabel) {
        for (size_t j = 0; j < meshVertexNb; ++j) {
          targetLabel = targetTex0[j];
          if (targetLabel>0 && targetLabel<= targetRegionsNb) {
            extractConnMatrix[vertexCount][targetLabel-1]
              += allMeshConnMatrixToAllMesh[i][j];
          }
        }
        (** seedVertexIndex)[vertexCount]=i;
        vertexCount++;
      }
    }
    return extractConn_ptr;
  }

  //---------------------------------------------------
  //  connMatrixSeedMesh_to_targetMeshTargets_regroup
  //---------------------------------------------------
  /*
  Write connectivity matrix of the vertex of the seedMesh with a grouping
  according to targetRegions (in targetMesh) given in targetRegionsTexture
  (labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to
  targetRegionsNb = target regions
  input: Connectivities * , shape = (seedMeshVertexNb, targetMeshVertexNb),
         generated by connMatrixSeedMeshToTargetMesh Method
  output: Connectivities * connMatrixSeedMeshToTargetMeshTargets_ptr,
          shape = (seedMeshVertexNb, targetMeshTargetsRegionsNb)
  */
  Connectivities *connMatrixSeedMesh_to_targetMeshTargets_regroup(
      const SparseOrDenseMatrix &connMatrixToAllMesh,
      const TimeTexture<short> &targetRegionsTex, int targetRegionsNb,
      bool verbose) {
    if (verbose) cout << "Start connMatrixRegionExtractRegroup" << endl;

    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    size_t targetMeshVertexNb = connMatrixToAllMesh.getSize2(); //nb columns
    size_t seedMeshVertexNb = connMatrixToAllMesh.getSize1(); //nb rows
    if (verbose)
      cout << "("<< seedMeshVertexNb << "," << targetMeshVertexNb << ") : "
        << flush;
    if (targetTex0.nItem()!=targetMeshVertexNb)
      throw runtime_error( "error in matrix dimensions" );
    Connectivities *connMatrixSeedMeshToTargetMeshTargets_ptr
      = new Connectivities(seedMeshVertexNb,
                           Connectivity(targetRegionsNb));
    Connectivities &connMatrixSeedMeshToTargetMeshTargets
      = *connMatrixSeedMeshToTargetMeshTargets_ptr;
    int targetLabel;
    //int five_count = int(seedMeshVertexNb / 20.0);
    for( size_t i = 0; i < seedMeshVertexNb; ++i ) {
      for(size_t j = 0; j < targetMeshVertexNb; ++j) {
        targetLabel = targetTex0[j];
        if (targetLabel>0 && targetLabel<= targetRegionsNb) {
          connMatrixSeedMeshToTargetMeshTargets[i][targetLabel-1]
            += connMatrixToAllMesh(i,j);
        }
      }
    }
    return connMatrixSeedMeshToTargetMeshTargets_ptr;
  }

  //----------------------------------------------------------
  //  connMatrixSeedMeshRegions_to_targetMeshTargets_regroup
  //----------------------------------------------------------
  /*
  Write connectivity matrix of regions of the seedMesh with a grouping
  according to targetRegions (in targetMesh) given in targetRegionsTexture
  (labeled texture, between 0 and targetRegionsNb, 0 = background, 1 to
  targetRegionsNb = target regions) (and a groupin according to seedRegions
  (in seedMesh) given in seedRegionsTexture)
  input: Connectivities * connMatrixSeedMeshToTargetMesh_ptr,
         shape = (seedMeshVertexNb, targetMeshVertexNb),
         generated by connMatrixSeedMeshToTargetMesh Method
  output: Connectivities * connMatrixSeedMeshRegionsToTargetMeshTargets_ptr,
          shape = (seedMeshVertexNb, targetMeshTargetsRegionsNb)
  */
  Connectivities * connMatrixSeedMeshRegions_to_targetMeshTargets_regroup(
      Connectivities * connMatrixSeedMeshToTargetMesh_ptr,
      const TimeTexture<short> & seedRegionsTex,
      const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
      int seedRegionsNb, bool verbose) {

    Connectivities & connMatrixSeedMeshToTargetMesh
      = *connMatrixSeedMeshToTargetMesh_ptr;
    const Texture<short> & targetTex0 = targetRegionsTex.begin()->second;
    const Texture<short> & seedTex0 = seedRegionsTex.begin()->second;
    size_t targetMeshVertexNb = connMatrixSeedMeshToTargetMesh[0].size();
    size_t seedMeshVertexNb = connMatrixSeedMeshToTargetMesh.size();
    if (verbose) cout << "seedMeshVertexNb:" << seedMeshVertexNb << endl;
    if (targetTex0.nItem()!=targetMeshVertexNb
        or seedTex0.nItem()!=seedMeshVertexNb)
      throw runtime_error("error in matrix dimensions");
    Connectivities *connMatrixSeedMeshRegionsToTargetMeshTargets_ptr
      = new Connectivities(seedRegionsNb,
                           Connectivity(targetRegionsNb));
    Connectivities &connMatrixSeedMeshRegionsToTargetMeshTargets
      = *connMatrixSeedMeshRegionsToTargetMeshTargets_ptr;
    int targetLabel;
    int seedRegionLabel;
    int five_count = int(seedMeshVertexNb / 20.0);
    for (size_t i = 0; i < seedMeshVertexNb; ++i) {
      if (i % five_count == 0) {
        if (verbose) cout << 5 * int(i / five_count) << "%..." << flush;
      }
      seedRegionLabel = seedTex0[i];
      if (seedRegionLabel > 0 && seedRegionLabel <= seedRegionsNb) {
        for (size_t j = 0; j < targetMeshVertexNb; ++j) {
          targetLabel = targetTex0[j];
          if (targetLabel > 0 && targetLabel <= targetRegionsNb) {
            connMatrixSeedMeshRegionsToTargetMeshTargets[
              seedRegionLabel-1][targetLabel-1]
                += connMatrixSeedMeshToTargetMesh[i][j];
          }
        }
      }
    }
    return connMatrixSeedMeshRegionsToTargetMeshTargets_ptr;
  }

  //--------------------------
  //  writeAimsFmtConnMatrix
  //--------------------------
  void writeAimsFmtConnMatrix(
      Connectivities *connMatrix_ptr, string file_name, bool verbose) {
    if (verbose) cout << "Writing connectivity matrix ";
    typedef Volume< float > Matrix;
    Connectivities &connMatrix  = *connMatrix_ptr;
    size_t colNb = connMatrix[0].size();
    size_t rowsNb = connMatrix.size();
    if (verbose) cout << "("<< rowsNb << ", " << colNb <<") in Aims Format...";
    Matrix matrix(rowsNb, colNb, 1);
    matrix.setVoxelSize(100.0/rowsNb , 80.0/colNb);
    for (size_t i = 0; i < rowsNb; ++i)
    {
      for (size_t j = 0; j < colNb; ++j)
      {
        matrix(i, j, 0) = connMatrix[i][j];
      }
    }
    Writer<Matrix> w( file_name );
    w.write( matrix );
    if (verbose) cout << "done." << endl;
  }

  //----------------------
  // connMatrixNormalize
  //----------------------
  void connMatrixNormalize(SparseOrDenseMatrix &inConnMatrix, bool verbose) {
    if (verbose) cout << "Normalize connectivity matrix " << flush;
    size_t colNb = inConnMatrix.getSize2();
    size_t rowsNb = inConnMatrix.getSize1();
    if (verbose) cout << "("<<rowsNb << "," << colNb << ") : " << flush;
    for (size_t i = 0; i < rowsNb; ++i) {
      vector<double> line_i = inConnMatrix.getSparseRow<vector<double> >(i);
      double norm_i = 0;
      for(size_t j = 0; j < colNb; ++j) norm_i += line_i[j] * line_i[j];
      norm_i = sqrt(norm_i);
      if (norm_i != 0) {
        for (size_t j = 0; j < colNb; ++j) {
          float connval = inConnMatrix(i, j);
          if (connval != 0) inConnMatrix.set_element(i, j, connval / norm_i);
        }
      }
    }
    if (verbose) cout << "Done." << endl;
  }

  //------------------
  //  densityTexture
  //------------------
  /*
  input:  connMatrixToAllMesh: connectivity matrix of shape (p, n),
                               p = seed region vertex nb, n = mesh vertex nb
          vertexIndex: index of the seed region vertex (matrix rows) in the
                       main mesh (vertexIndex[i]=index corresponding to the
                       ith row of the connectivity matrix in the mesh)
  output: outputDensityTex: texture of connection density: for each vertex of
                            the seed region : sum of its connections to all the
                            other mesh vertex
    */
  TimeTexture<float> densityTexture(
      Connectivities *connMatrixToAllMesh_ptr, vector<size_t> vertexIndex,
      bool verbose) {
    if (verbose) cout << "Computing density texture" << flush;

    Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
    size_t colNb = connMatrixToAllMesh[0].size();
    size_t rowsNb = connMatrixToAllMesh.size();
    if (verbose) cout << "("<<rowsNb << "," << colNb << ") : " << flush;
    if (rowsNb != vertexIndex.size())
      throw runtime_error("Not the same dimensions as vertexIndex vector");
    Connectivity sparseVectorIdentity = Connectivity(colNb);
    for(size_t j = 0; j < colNb; ++j) sparseVectorIdentity[j]=1.;
    int five_count = int(rowsNb / 20.0);
    float currentSum;
    int currentVertexIndex;
    TimeTexture<float> outputDensityTex;
    int countVertex = 0;
    outputDensityTex[0].reserve(colNb);
    for (size_t v = 0; v < colNb; ++v) outputDensityTex[0].push_back(-1);
    for (size_t i = 0; i < rowsNb; ++i) {
      currentVertexIndex = vertexIndex[i];
      if (i % five_count == 0) {
        if (verbose) cout << 5 * int(i / five_count) << "%..." << flush;
      }
      outputDensityTex[0][currentVertexIndex] = 0;
      Connectivity & line_i = connMatrixToAllMesh[i];
      if( !line_i.empty() )
      {
        currentSum = line_i.dot( sparseVectorIdentity );
        if (currentSum == -1 and verbose) cout << "-1," << flush;
        if (isnan(currentSum) and verbose)
          cout << "currentSum is nan!!" << flush;
        if (currentSum > 0) {
          outputDensityTex[0][currentVertexIndex] = currentSum;
        }
      }
      countVertex++;
    }
    if (verbose)
      cout << "countVertex++:" << countVertex++ << ", Done." << endl;
    return outputDensityTex;
  }

  //------------------
  //  densityTexture
  //------------------
  /*
  input:  connMatrixToAllMesh : connectivity matrix of shape (p, n),
                                p = seed region vertex nb, n = mesh vertex nb
          vertexIndex : index of the seed region vertex (matrix rows) in the
                        main mesh (vertexIndex[i]=index corresponding to the
                        ith row of the connectivity matrix in the mesh)
  output: outputDensityTex : texture of connection density : for each vertex of
                             the seed region : sum of its connections to all
                             the other mesh vertex
    */
  TimeTexture<float> densityTexture(
      const SparseOrDenseMatrix & connMatrixToAllMesh,
      vector<size_t> vertexIndex, bool verbose) {
    if (verbose) cout << "Computing density texture" << flush;
    size_t colNb = connMatrixToAllMesh.getSize2();
    size_t rowsNb = connMatrixToAllMesh.getSize1();
    if (verbose) cout << "("<<rowsNb << "," << colNb << ") : " << flush;
    if (rowsNb != vertexIndex.size())
      throw runtime_error("Not the same dimensions as vertexIndex vector");
    int five_count = int(rowsNb / 20.0);
    float currentSum = 0.0f;
    int currentVertexIndex;
    TimeTexture<float> outputDensityTex;
    int countVertex = 0;
    outputDensityTex[0].reserve(colNb);
    for (size_t v = 0; v < colNb; ++v) outputDensityTex[0].push_back(-1);
    for (size_t i = 0; i < rowsNb; ++i) {
      currentVertexIndex = vertexIndex[i];
      if (i % five_count == 0) {
        if (verbose) cout << 5 * int(i / five_count) << "%..." << flush;
      }
      outputDensityTex[0][currentVertexIndex]=0;
      vector<double> line_i = connMatrixToAllMesh.getRow(i);
      for (size_t j=0; j<colNb; ++j) currentSum += line_i[j];
      if ((currentSum == -1) && verbose) cout <<"-1," << flush;
      if (isnan(currentSum) && verbose)
        cout << "currentSum is nan!!" << flush;
      if (currentSum > 0) {
        outputDensityTex[0][currentVertexIndex] = (float) currentSum;
      }
      countVertex++;
    }
    if (verbose)
      cout << "countVertex++:" << countVertex++ << ", Done." << endl;
    return outputDensityTex;
  }

  //----------------------
  //  meshDensityTexture
  //----------------------
  /*
  input:  connMatrixToAllMesh: connectivity matrix of shape (p, n),
          p = seed region vertex nb, n = mesh vertex nb
  output: outputDensityTex: texture of connection density : for each vertex of
          the mesh : sum of its connections to all seed region vertex
  */
  TimeTexture<float> meshDensityTexture(
      Connectivities * connMatrixToAllMesh_ptr, bool verbose) {
    if (verbose) cout << "Computing Mesh density texture" << flush;
    Connectivities &connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
    size_t colNb = connMatrixToAllMesh[0].size();
    size_t rowsNb = connMatrixToAllMesh.size();
    if (verbose) cout << "("<<rowsNb << "," << colNb << ") : " << flush;
    TimeTexture<float> outputDensityTex;
    outputDensityTex[0].reserve(colNb);
    Connectivity *totalConnectionDensity_ptr = connMatrixSumRows(
       connMatrixToAllMesh_ptr);
    Connectivity &totalConnectionDensity = *totalConnectionDensity_ptr;
    for (size_t v = 0; v < colNb; ++v) outputDensityTex[0].push_back(-1);
    for (size_t j = 0; j < colNb; ++j) {
      outputDensityTex[0][j]= totalConnectionDensity[j];
    }
    if (verbose) cout << "Done." << endl;
    delete totalConnectionDensity_ptr;
    return outputDensityTex;
  }


  //----------------------
  //  meshDensityTexture
  //----------------------
  /*
  input: connMatrixToAllMesh: connectivity matrix of shape (p, n),
                              p = seed region vertex nb, n = mesh vertex nb
  output: outputDensityTex: texture of connection density : for each vertex of
          the mesh: sum of its connections to all seed region vertex
  */
  TimeTexture<float> meshDensityTexture(
      const SparseOrDenseMatrix &connMatrixToAllMesh, bool verbose)
  {
    if (verbose) cout << "Computing Mesh density texture" << flush;

    size_t colNb = connMatrixToAllMesh.getSize2();
    size_t rowsNb = connMatrixToAllMesh.getSize1();
    if (verbose) cout << "("<<rowsNb << "," << colNb << ") : " << flush;
    TimeTexture<float> outputDensityTex;
    outputDensityTex[0].reserve(colNb);
    vector<double> *totalConnectionDensity_ptr = connMatrixSumRows(
        connMatrixToAllMesh);
    vector<double> &totalConnectionDensity = *totalConnectionDensity_ptr;
    for (size_t v = 0; v < colNb; ++v) outputDensityTex[0].push_back(-1);
    for (size_t j = 0; j < colNb; ++j) {
      outputDensityTex[0][j]= totalConnectionDensity[j];
    }
    if (verbose) cout << "Done." << endl;
    delete totalConnectionDensity_ptr;
    return outputDensityTex;
  }


  template <typename T>
  double get_item( const T& container, size_t item )
  {
    return container[item];
  }

  template <>
  double get_item( const Connectivity & container, size_t item )
  {
    return container.get( item );
  }

  //-------------------------------------------
  //  oneTargetDensityTargetsRegroupTexture_T
  //-------------------------------------------
  /*
  Compute density texture (associated to a mesh) of a targetRegion towards
  all the target regions
  inputs: lineMatrixToTargetsRegions_ptr: connectivity of one target to all
          the target regions, shape = (1,targetRegionsNb)
          targetRegionsTex: labeled texture, between 0 and targetRegionsNb,
          0 = background, 1 to targetRegionsNb = target regions
  output: TimeTexture<float> *outputDensityTex: pointer on output density
          texture
  */
  template <typename T>
  TimeTexture<float> *oneTargetDensityTargetsRegroupTexture_T(
      const T *lineMatrixToTargetRegions_ptr,
      const TimeTexture<short> &targetRegionsTex, int timestep) {
    TimeTexture<short>::const_iterator it = targetRegionsTex.find(timestep);
    if (it == targetRegionsTex.end()) it = targetRegionsTex.begin();
    const Texture<short> & targetTex0 = it->second;
    size_t meshVertexNb = targetTex0.nItem();
    const T &lineMatrixToTargetRegions = *lineMatrixToTargetRegions_ptr;
    size_t targetRegionsNb = lineMatrixToTargetRegions.size();
    TimeTexture<float> *outputDensityTex_ptr = new TimeTexture<float>;
    TimeTexture<float> &outputDensityTex = *outputDensityTex_ptr;
    outputDensityTex[0].reserve(meshVertexNb);
    for (size_t v = 0; v < meshVertexNb; ++v)
      outputDensityTex[0].push_back(-1);

    for (size_t i = 0; i < meshVertexNb; i++)
    {
      int label = targetTex0[i];//equivalent to targetRegionsTex[0].item(i)
      if (( 1 <= label ) && ( label <= int(targetRegionsNb) ))
      {
        outputDensityTex[0][i]
          = (float) get_item( lineMatrixToTargetRegions, label-1 );
      }
      else
      {
        outputDensityTex[0][i] = -1;
      }
    }
    return outputDensityTex_ptr;
  }

  //-----------------------------------------
  //  oneTargetDensityTargetsRegroupTexture
  //-----------------------------------------
  TimeTexture<float> *oneTargetDensityTargetsRegroupTexture(
      const Connectivity *lineMatrixToTargetRegions_ptr,
      const TimeTexture<short> &targetRegionsTex, int timestep)
  {
    return oneTargetDensityTargetsRegroupTexture_T(
      lineMatrixToTargetRegions_ptr, targetRegionsTex, timestep);
  }

  //-----------------------------------------
  //  oneTargetDensityTargetsRegroupTexture
  //-----------------------------------------
  TimeTexture<float> *oneTargetDensityTargetsRegroupTexture(
      const vector<double> *lineMatrixToTargetRegions_ptr,
      const TimeTexture<short> &targetRegionsTex, int timestep)
  {
    return oneTargetDensityTargetsRegroupTexture_T(
      lineMatrixToTargetRegions_ptr, targetRegionsTex, timestep);
  }


  //---------------------------------------
  //  connectivitiesToSparseOrDenseMatrix
  //---------------------------------------
  SparseOrDenseMatrix *connectivitiesToSparseOrDenseMatrix(
      const Connectivities &conn) {
    unsigned ncol = 0;
    unsigned line, nline = conn.size();
    unsigned long count = 0;

    if (!conn.empty()) ncol = conn.begin()->size();
    SparseOrDenseMatrix *mat = new SparseOrDenseMatrix(conn.size(), ncol);
    for (line = 0; line < nline; ++line) count += conn[line].size();
    if (count >= mat->optimalShapeThreshold()) mat->muteToDense();

    Connectivity::const_iterator ic, ec;

    for (line = 0; line < nline; ++line)
    {
      const Connectivity & cline = conn[line];
      if (!cline.empty())
      {
        for (ic = cline.begin(), ec = cline.end();
             ic != ec; ++ic)
        {
          mat->set_element(line, ic->first, ic->second);
        }
      }
    }
    return mat;
  }


} // namespace constel
