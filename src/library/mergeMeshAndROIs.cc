#include <constellation/mergeMeshAndROIs.h>
#include <connectomist/fibertracking/bundleRegroup.h>
#include <constellation/selectFiberListenerFromMesh.h>
#include <constellation/bundleTools.h>
#include <connectomist/fibertracking/bundleSampler.h>
#include <cartobase/config/verbose.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;
using namespace comist;
using namespace constel;


namespace constel {
  
  //---------------------------------
  //  fusionMeshROIs_to_Texture_S16
  //---------------------------------
  TimeTexture<short> *fusionMeshROIs_to_Texture_S16(
      const AimsSurfaceTriangle &inAimsMesh, const Graph &ROI_graph,
      float distthresh, string vertex_property) {
    /*
    input:    inAimsMesh : .mesh
              Roi_buckets : .arg (ROI graph)
    output:   texture from the merging of mesh and Roi_buckets
    */
    
    cout << "Computing Fusion Texture" << endl;
    //Converting AimsMesh to Cathier format
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);
    
    // Generating kdtree
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    
    vector<float> voxel_size(3);
    
    if (ROI_graph.getProperty("voxel_size", voxel_size)) {
      //Computing texture
      TimeTexture<short> *tex_ptr = new TimeTexture<short>;
      TimeTexture<short> &tex = *tex_ptr;
      size_t meshVertexNb = getVertices(mesh).size();
      tex[0].reserve(meshVertexNb);
      TimeTexture<float> distROIToMeshTex;
      distROIToMeshTex[0].reserve(meshVertexNb);
      for (size_t i = 0; i < meshVertexNb; i++) {
        tex[0].push_back(0);
        distROIToMeshTex[0].push_back(distthresh);
      }
      cout << "voxel_size:" << voxel_size[0] << ", " << voxel_size[1] << ", "
        << voxel_size[2] << endl;
      const set<Vertex*> vset = ROI_graph.vertices();
      int roi_nb = vset.size();
      cout << "roi nb:" << roi_nb << endl;
      set<Vertex*>::const_iterator roi_it, roi_it_end = vset.end();
      rc_ptr<BucketMap<Void> > roi_vertex_buckmap;
      int tex_roiLabel;
      AimsVector<short int, 3> p;
      til::numeric_array<float, 3> pFloat;
      string roi_name;
      int roi_count = 1;
      cout << "Iteration on: " << flush;
      for(roi_it = vset.begin(); roi_it != roi_it_end; roi_it++) {
        const Vertex* roi_vertex = *roi_it;
        if (!(roi_vertex->getProperty("roi_label",tex_roiLabel))) {
          tex_roiLabel = roi_count;
        }
        roi_vertex->getProperty("name", roi_name);
        cout << roi_name << ", roiLabel:" << tex_roiLabel << ", " << flush;
        if (roi_vertex->getProperty(vertex_property, roi_vertex_buckmap)) {
          BucketMap<Void>::Bucket bck = roi_vertex_buckmap->begin()->second;
          BucketMap<Void>::Bucket::iterator itbck, itbck_end = bck.end();
//       Looking for closest points
          til::Find_closest< double, KDTree > fc(kdt);
          for (itbck = bck.begin(); itbck != itbck_end; itbck++) {
            p = itbck->first;
            til::numeric_array<float, 3> pFloat(p[0]*voxel_size[0],
                                                p[1]*voxel_size[1],
                                                p[2]*voxel_size[2]);
            size_t pIndex_Mesh = fc(pFloat);
            float distROIToMeshVertex = til::dist2(
                pFloat, getVertices(mesh)[pIndex_Mesh], til::prec<float>());
            if (distROIToMeshVertex <= distROIToMeshTex[0][pIndex_Mesh]) {
              tex[0][pIndex_Mesh] = tex_roiLabel;
              distROIToMeshTex[0][pIndex_Mesh] = distROIToMeshVertex;
            }
          }
        } else {
          cout << roi_name << " empty" << endl;
        }
        roi_count++;
      }
      return tex_ptr;
    } else {
      cout << "error on ROI arg voxel_size" << endl;
      return 0;
    }
  }

  //----------------------------------
  //  fusionMeshFoldGraph_to_Texture
  //----------------------------------
  TimeTexture<short> *fusionMeshFoldGraph_to_Texture(
      const AimsSurfaceTriangle &inAimsMesh, const Graph &ROI_graph,
      float distthresh, string vertex_property, string label_name) {
    /*
    input:    inAimsMesh : .mesh
              Roi_buckets : .arg (fold graph) with labeled sulci
    output:   texture from the merging of mesh and Roi_buckets
    */
    
    cout << "Computing Fusion Texture" << endl;
    //Converting AimsMesh to Cathier format
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);
    
    // Generating kdtree
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    
    vector<float> voxel_size(3);

    if (ROI_graph.getProperty("voxel_size", voxel_size)) {
      //Computing texture
      TimeTexture<short> * tex_ptr = new TimeTexture<short>;
      TimeTexture<short> & tex = *tex_ptr;
      size_t meshVertexNb = getVertices(mesh).size();
      tex[0].reserve(meshVertexNb);
      TimeTexture<float> distROIToMeshTex;
      distROIToMeshTex[0].reserve(meshVertexNb);
      for (size_t i = 0; i < meshVertexNb; i++) {
        tex[0].push_back(0);
        distROIToMeshTex[0].push_back(distthresh);
      }
      cout << "voxel_size:" << voxel_size[0] << ", " << voxel_size[1] << ", "
        << voxel_size[2] << endl;
      const set<Vertex*> vset = ROI_graph.vertices();
      int roi_nb = vset.size();
      cout << "roi nb:" << roi_nb << endl;
      set<Vertex*>::const_iterator roi_it, roi_it_end = vset.end();
      rc_ptr<BucketMap<Void> > roi_vertex_buckmap;
      int tex_roiLabel;
      AimsVector<short int, 3> p;
      til::numeric_array<float, 3> pFloat;
      string roi_name;
      int roi_count = 1;
      cout << "Iteration on: " << flush;
      for (roi_it = vset.begin(); roi_it != roi_it_end; roi_it++) {
        const Vertex* roi_vertex = *roi_it;
        if (!(roi_vertex->getProperty(label_name,tex_roiLabel))) {
          tex_roiLabel = roi_count;
        }
        roi_vertex->getProperty("name", roi_name);
        cout << roi_name << ", roiLabel name :" << label_name << ":"
          << tex_roiLabel << ", " << flush;
        if (roi_vertex->getProperty(vertex_property, roi_vertex_buckmap)) {
          BucketMap<Void>::Bucket bck = roi_vertex_buckmap->begin()->second;
          BucketMap<Void>::Bucket::iterator itbck, itbck_end = bck.end();
//       Looking for closest points
          til::Find_closest< double, KDTree > fc(kdt);
          for (itbck = bck.begin(); itbck != itbck_end; itbck++) {
            p = itbck->first;
            til::numeric_array<float, 3> pFloat(p[0]*voxel_size[0],
                                                p[1]*voxel_size[1],
                                                p[2]*voxel_size[2]);
            size_t pIndex_Mesh = fc(pFloat);
            float distROIToMeshVertex = til::dist2(
                pFloat, getVertices(mesh)[pIndex_Mesh], til::prec<float>());
            if (distROIToMeshVertex <= distROIToMeshTex[0][pIndex_Mesh]) {
              tex[0][pIndex_Mesh] = tex_roiLabel;
              distROIToMeshTex[0][pIndex_Mesh] = distROIToMeshVertex;
            }
          }
        } else {
          cout << roi_name << " empty" << endl;
        }
        roi_count++;
      }
      return tex_ptr;
    } else {
      cout << "error on ROI arg voxel_size" << endl;
      return 0;
    }
  }

  //-------------------------------------
  //  texMeshAndBundles_to_BundlesGraph
  //-------------------------------------
  Graph *texMeshAndBundles_to_BundlesGraph(
      const AimsSurfaceTriangle &inAimsMesh, rc_ptr<TimeTexture<short> > tex,
      string namesMode, string bundlesFile_name,
      Motion motion, const string &BundlesNamesfile_name,
      float filter_proportion, int texture_time_step ) {
    //Converting AimsMesh to Cathier format (neighbouring points)
    rc_ptr<Mesh> mesh( new Mesh );
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    *mesh = addNeighborsToMesh(mesh0);
    stringstream osBundlesNamesfile_name;
    SelectFiberListenerFromMesh fiberBundle( mesh, tex, namesMode, 0, motion,
      BundlesNamesfile_name, texture_time_step );
    if (BundlesNamesfile_name.empty())
      fiberBundle.setStream(osBundlesNamesfile_name);
    BundleReader bundle(bundlesFile_name);
    bundle.addBundleListener(fiberBundle);
    bundle.read();

    BundleReader bundle2(bundlesFile_name);
    rc_ptr< BundleRegroup > bundleRegroup;
    bundleRegroup.reset(new BundleRegroup(BundlesNamesfile_name));
    if (BundlesNamesfile_name.empty())
      bundleRegroup->setStream(osBundlesNamesfile_name);
    bundle2.addBundleListener(*bundleRegroup);

    rc_ptr< BundleSampler > bundleSampler;
    if (filter_proportion != 0.) {
      bundleSampler.reset(new BundleSampler(filter_proportion * 100, "toto",
                                             "tutu", 0));
      bundleRegroup->addBundleListener(*bundleSampler);
    }
    Graph *result = new Graph("RoiArg");
    BundleToGraph bundleToGraph(*result);
    if (filter_proportion != 0.)
      bundleSampler->addBundleListener(bundleToGraph);
    else
      bundleRegroup->addBundleListener(bundleToGraph);
    bundle2.read();

    return result;
  }

  //---------------------------------------------------------------
  //  texMeshAndBundles_to_BundlesGraph_WithIntersectionComputing
  //---------------------------------------------------------------
  Graph * texMeshAndBundles_to_BundlesGraph_WithIntersectionComputing(
      const AimsSurfaceTriangle &inAimsMesh, const TimeTexture<short> &tex,
      string bundlesFile_name, Motion motion, string BundlesNamesfile_name) {
    double meshDistanceThreshold = 0.;
    double meshClosestPointMaxDistance = 5.;  //mm

    //Computing fiber and mesh intersections
    //BundleReader creation:
    constel::BundleInteractionReader bundleInteractionReader(bundlesFile_name);
    BundleProducer *finalProducer = &bundleInteractionReader;
    
    //Applying transformation from t2 to anat to the fibers
    rc_ptr< BundleMotion > bundleMotion;
    bundleMotion.reset(new BundleMotion(motion));
    finalProducer->addBundleListener(*bundleMotion);
    finalProducer = bundleMotion.get();
    
    // Initializing BundleListener which memorizes the anterior fiber point
    // during the fiber's reading:
    constel::MemAntBundleListener memAntBundleListener(
        bundleInteractionReader);

    // BundleListener which computes the curvilinear abscissa of current fiber
    // point during the fiber's reading:
    constel::CurvilinearAbscissaBundleListener
    curvilinearAbscissaBundleListener(bundleInteractionReader);

    // BundleListener which gives a name to the fiber, according to its
    // intersections with the input meshes:
    constel::FiberNameAccordingToMeshTextureIntersectionBundleListener
    fiberNameAccordingToMeshTextureIntersectionBundleListener(
        bundleInteractionReader, tex, BundlesNamesfile_name);

    // Initializing for each input mesh the BundleListener for computing
    // fiber's intersection with the mesh
    constel::MeshIntersectionBundleListener meshIntersectionBundleListener(
        inAimsMesh, bundleInteractionReader, meshDistanceThreshold,
        meshClosestPointMaxDistance);
    meshIntersectionBundleListener.setMeshIdentity(0);

    //Additions to the BundleReader:
    finalProducer->addBundleListener(curvilinearAbscissaBundleListener);
    finalProducer->addBundleListener(meshIntersectionBundleListener);
    finalProducer->addBundleListener(memAntBundleListener);
    finalProducer->addBundleListener(
        fiberNameAccordingToMeshTextureIntersectionBundleListener);
    bundleInteractionReader.read();
    
    BundleReader bundle2(bundlesFile_name);
    rc_ptr< BundleRegroup > bundleRegroup;
    bundleRegroup.reset(new BundleRegroup(BundlesNamesfile_name));
    bundle2.addBundleListener(*bundleRegroup);

    Graph *result = new Graph("RoiArg");
    BundleToGraph bundleToGraph(*result);
    bundleRegroup->addBundleListener(bundleToGraph);
    bundle2.read();
    
    return result;
  }

} // namespace constel

