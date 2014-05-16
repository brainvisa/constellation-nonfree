#include <roca/mergeMeshAndROIs.h>
#include <connectomist/fibertracking/bundleRegroup.h>
#include <constellation/selectFiberListenerFromMesh.h>
#include <constellation/bundleTools.h>
#include <cartobase/config/verbose.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;
using namespace comist;
using namespace constel;


namespace constel {
  
  TimeTexture<short> * fusionMeshROIs_to_Texture_S16(const AimsSurfaceTriangle & inAimsMesh, const Graph & ROI_graph, float distthresh, std::string vertex_property)
  {
    /*
    input:    inAimsMesh : .mesh
              Roi_buckets : .arg (ROI graph)
    output:   texture from the merging of mesh and Roi_buckets
    */
    
    std::cout << "Computing Fusion Texture" << std::endl;
    //Converting AimsMesh to Cathier format
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);
    
    // Generating kdtree
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    
    std::vector<float> voxel_size(3);
    
    if(ROI_graph.getProperty("voxel_size", voxel_size))
    {
      //Computing texture
      TimeTexture<short> * tex_ptr = new TimeTexture<short>;
      TimeTexture<short> & tex = *tex_ptr;
      std::size_t meshVertexNb = getVertices(mesh).size();
      tex[0].reserve(meshVertexNb);
      TimeTexture<float> distROIToMeshTex;
      distROIToMeshTex[0].reserve(meshVertexNb);
      for (std::size_t i = 0; i < meshVertexNb; i++)
      {
        tex[0].push_back(0);
        distROIToMeshTex[0].push_back(distthresh);
      }
      std::cout << "voxel_size:" << voxel_size[0] << ", " << voxel_size[1] << ", " << voxel_size[2] << std::endl;
      const set<Vertex*> vset = ROI_graph.vertices();
      int roi_nb = vset.size();
      std::cout << "roi nb:" << roi_nb << std::endl;
      set<Vertex*>::const_iterator roi_it, roi_it_end = vset.end();
      rc_ptr<BucketMap<Void> > roi_vertex_buckmap;
      int tex_roiLabel;
      AimsVector<short int, 3> p;
      til::numeric_array<float, 3> pFloat;
      std::string roi_name;
      int roi_count = 1;
      std::cout << "Iteration on: " <<std::flush;
      for(roi_it = vset.begin(); roi_it != roi_it_end; roi_it++)
      {
        
        const Vertex* roi_vertex = *roi_it;
        if (!(roi_vertex->getProperty("roi_label",tex_roiLabel)))
        {
          tex_roiLabel = roi_count;
        }
        roi_vertex->getProperty("name", roi_name);
        std::cout << roi_name << ", roiLabel:" << tex_roiLabel << ", " << std::flush;
        if (roi_vertex->getProperty(vertex_property, roi_vertex_buckmap))
        {
          BucketMap<Void>::Bucket bck = roi_vertex_buckmap->begin()->second;
          BucketMap<Void>::Bucket::iterator itbck, itbck_end = bck.end();
//       Looking for closest points
          til::Find_closest< double, KDTree > fc(kdt);
          for (itbck = bck.begin(); itbck != itbck_end; itbck++)
          {
            p = itbck->first;
            til::numeric_array<float, 3> pFloat(p[0]*voxel_size[0], p[1]*voxel_size[1], p[2]*voxel_size[2]);
//             std::cout << pFloat << ", " << std::flush;
            std::size_t pIndex_Mesh = fc(pFloat);
            float distROIToMeshVertex = til::dist2(pFloat, getVertices(mesh)[pIndex_Mesh], til::prec<float>());
            if (distROIToMeshVertex <= distROIToMeshTex[0][pIndex_Mesh])
            {
              tex[0][pIndex_Mesh] = tex_roiLabel;
              distROIToMeshTex[0][pIndex_Mesh] = distROIToMeshVertex;
            }
          }
        }
        else
        {
          std::cout << roi_name << " empty" << std::endl;
        }
        roi_count++;
      }
      std::cout << "..., OK." << std::endl;
      return tex_ptr;
    }
    else
    {
      std::cout << "error on ROI arg voxel_size" << std::endl;
      return 0;
    }
    
  }//fusionMeshRois_to_Texture_S16
  
  TimeTexture<short> * fusionMeshFoldGraph_to_Texture(const AimsSurfaceTriangle & inAimsMesh, const Graph & ROI_graph, float distthresh, std::string vertex_property, std::string label_name)
  {
    /*
    input:    inAimsMesh : .mesh
              Roi_buckets : .arg (fold graph) with labeled sulci
    output:   texture from the merging of mesh and Roi_buckets
    */
    
    std::cout << "Computing Fusion Texture" << std::endl;
    //Converting AimsMesh to Cathier format
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);
    
    // Generating kdtree
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    
    std::vector<float> voxel_size(3);

    if(ROI_graph.getProperty("voxel_size", voxel_size))
    {
      //Computing texture
      TimeTexture<short> * tex_ptr = new TimeTexture<short>;
      TimeTexture<short> & tex = *tex_ptr;
      std::size_t meshVertexNb = getVertices(mesh).size();
      tex[0].reserve(meshVertexNb);
      TimeTexture<float> distROIToMeshTex;
      distROIToMeshTex[0].reserve(meshVertexNb);
      for (std::size_t i = 0; i < meshVertexNb; i++)
      {
        tex[0].push_back(0);
        distROIToMeshTex[0].push_back(distthresh);
      }
      std::cout << "voxel_size:" << voxel_size[0] << ", " << voxel_size[1] << ", " << voxel_size[2] << std::endl;
      const set<Vertex*> vset = ROI_graph.vertices();
      int roi_nb = vset.size();
      std::cout << "roi nb:" << roi_nb << std::endl;
      set<Vertex*>::const_iterator roi_it, roi_it_end = vset.end();
      rc_ptr<BucketMap<Void> > roi_vertex_buckmap;
      int tex_roiLabel;
      AimsVector<short int, 3> p;
      til::numeric_array<float, 3> pFloat;
      std::string roi_name;
      int roi_count = 1;
      std::cout << "Iteration on: " <<std::flush;
      for(roi_it = vset.begin(); roi_it != roi_it_end; roi_it++)
      {
        
        const Vertex* roi_vertex = *roi_it;
        if (!(roi_vertex->getProperty(label_name,tex_roiLabel)))
        {
          tex_roiLabel = roi_count;
        }
        roi_vertex->getProperty("name", roi_name);
        std::cout << roi_name << ", roiLabel name :" << label_name << ":" << tex_roiLabel << ", " << std::flush;
        if (roi_vertex->getProperty(vertex_property, roi_vertex_buckmap))
        {
          BucketMap<Void>::Bucket bck = roi_vertex_buckmap->begin()->second;
          BucketMap<Void>::Bucket::iterator itbck, itbck_end = bck.end();
//       Looking for closest points
          til::Find_closest< double, KDTree > fc(kdt);
          for (itbck = bck.begin(); itbck != itbck_end; itbck++)
          {
            p = itbck->first;
            til::numeric_array<float, 3> pFloat(p[0]*voxel_size[0], p[1]*voxel_size[1], p[2]*voxel_size[2]);
//             std::cout << pFloat << ", " << std::flush;
            std::size_t pIndex_Mesh = fc(pFloat);
            float distROIToMeshVertex = til::dist2(pFloat, getVertices(mesh)[pIndex_Mesh], til::prec<float>());
            if (distROIToMeshVertex <= distROIToMeshTex[0][pIndex_Mesh])
            {
              tex[0][pIndex_Mesh] = tex_roiLabel;
              distROIToMeshTex[0][pIndex_Mesh] = distROIToMeshVertex;
            }
          }
        }
        else
        {
          std::cout << roi_name << " empty" << std::endl;
        }
        roi_count++;
      }
      std::cout << "..., OK." << std::endl;
      return tex_ptr;
    }
    else
    {
      std::cout << "error on ROI arg voxel_size" << std::endl;
      return 0;
    }
    
  }//fusionMeshFoldGraph_to_Texture

  Graph * texMeshAndBundles_to_BundlesGraph( const AimsSurfaceTriangle & inAimsMesh, rc_ptr<TimeTexture<short> > tex, string namesMode, string bundlesFile_name, Motion motion, const string & BundlesNamesfile_name)
  {
    //Converting AimsMesh to Cathier format (neighbouring points)
    rc_ptr<Mesh> mesh( new Mesh );
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    *mesh = addNeighborsToMesh(mesh0);
    //int addInt = 0; don't use
    stringstream osBundlesNamesfile_name;
    SelectFiberListenerFromMesh fiberBundle(mesh,tex,namesMode,0,motion,BundlesNamesfile_name);
    if( BundlesNamesfile_name.empty() )
      fiberBundle.setStream( osBundlesNamesfile_name );
    BundleReader bundle(bundlesFile_name);
    bundle.addBundleListener(fiberBundle);
    bundle.read();

    BundleReader bundle2(bundlesFile_name);
    rc_ptr< BundleRegroup > bundleRegroup;
    bundleRegroup.reset( new BundleRegroup(BundlesNamesfile_name) );
    if( BundlesNamesfile_name.empty() )
      bundleRegroup->setStream( osBundlesNamesfile_name );
    bundle2.addBundleListener(*bundleRegroup);

    Graph *result = new Graph("RoiArg");
    BundleToGraph bundleToGraph(*result);
    bundleRegroup->addBundleListener( bundleToGraph );
    bundle2.read();
    
    return result;
  }//texMeshAndROIs_to_BundlesGraph

  Graph * texMeshAndBundles_to_BundlesGraph_WithIntersectionComputing( const AimsSurfaceTriangle & inAimsMesh, const TimeTexture<short> & tex, string bundlesFile_name, Motion motion, string BundlesNamesfile_name)
  {
    double meshDistanceThreshold = 0.;
    double meshClosestPointMaxDistance = 5.;//mm

    //Computing fiber and mesh intersections:
    if(verbose) std::cout << "Loading fibers and computing intersections with the input mesh:" << std::endl;
    //BundleReader creation:
    constel::BundleInteractionReader bundleInteractionReader(bundlesFile_name);
    BundleProducer *finalProducer = &bundleInteractionReader;
    
    //Applying transformation from t2 to anat to the fibers:
    rc_ptr< BundleMotion > bundleMotion;
    if ( verbose ) cout << "Creating motion filter: " << flush;
    bundleMotion.reset( new BundleMotion( motion ) );
    finalProducer->addBundleListener( *bundleMotion );
    finalProducer = bundleMotion.get();
    if ( verbose ) cout << "done" << endl;
  
    //Adding the different BundleListener to the BundleReader:
    if (verbose)  std::cout << "Adding the different BundleListener to the BundleReader..." << std::flush;
    
    //Initializing BundleListener which memorizes the anterior fiber point during the fiber's reading:
    constel::MemAntBundleListener memAntBundleListener(bundleInteractionReader);

    //BundleListener which computes the curvilinear abscissa of current fiber point during the fiber's reading:
    constel::CurvilinearAbscissaBundleListener curvilinearAbscissaBundleListener(bundleInteractionReader);

    //BundleListener which gives a name to the fiber, according to its intersections with the input meshes:
    constel::FiberNameAccordingToMeshTextureIntersectionBundleListener fiberNameAccordingToMeshTextureIntersectionBundleListener(bundleInteractionReader, tex, BundlesNamesfile_name);

    //Initializing for each input mesh the BundleListener for computing fiber's intersection with the mesh
    constel::MeshIntersectionBundleListener meshIntersectionBundleListener(inAimsMesh, bundleInteractionReader, meshDistanceThreshold, meshClosestPointMaxDistance);
    meshIntersectionBundleListener.setMeshIdentity(0);

    //Additions to the BundleReader:
    finalProducer->addBundleListener(curvilinearAbscissaBundleListener);
    finalProducer->addBundleListener(meshIntersectionBundleListener);
    finalProducer->addBundleListener(memAntBundleListener);
    finalProducer->addBundleListener(fiberNameAccordingToMeshTextureIntersectionBundleListener);
    if (verbose)  std::cout << "Done." << std::endl;
    if (verbose)  std::cout << "Reading bundleInteractionReader: intersections computing and fiber\'s names filename saving in "<< BundlesNamesfile_name << "..." << std::flush;
    bundleInteractionReader.read();
    if (verbose)  std::cout << "Done." << std::endl;
    

    BundleReader bundle2(bundlesFile_name);
    rc_ptr< BundleRegroup > bundleRegroup;
    bundleRegroup.reset( new BundleRegroup(BundlesNamesfile_name) );
    bundle2.addBundleListener(*bundleRegroup);

    Graph *result = new Graph("RoiArg");
    BundleToGraph bundleToGraph(*result);
    bundleRegroup->addBundleListener( bundleToGraph );
    bundle2.read();
    
    return result;
  }//texMeshAndBundles_to_BundlesGraph_WithIntersectionComputing
/*  
  Graph * fusionMeshROIs_to_Graph(const AimsSurfaceTriangle & inAimsMesh, const Graph & ROI_graph, float distthresh)
  {
    
    //input:    inAimsMesh : .mesh
    //Roi_buckets : .arg (ROI graph)
    //output:   Graph from the "merging" of mesh and Roi_buckets:
      //        graph nodes : -submesh extracted from intersection of the Roi_bucket and the mesh if the distance between both structures is lower than distthresh
        //                    -aims_roi (bucket) otherwise.
    
    std::cout << "Computing Fusion Graph" << std::endl;
    //Converting AimsMesh to Cathier format
    Mesh mesh;
    til::Mesh1 mesh0;
    til::convert(mesh0, inAimsMesh);
    mesh = addNeighborsToMesh(mesh0);
    
    // Generating kdtree
    KDTree kdt(getVertices(mesh));
    makeKDTree(getVertices(mesh), kdt);
    
    std::vector<float> voxel_size(3);
    if(ROI_graph.getProperty("voxel_size", voxel_size))
    {
      //Computing Graph
      Graph * outputGraph_ptr = new Graph;
      Graph & outputGraph = *outputGraph_ptr;
      
      std::cout << "voxel_size:" << voxel_size[0] << ", " << voxel_size[1] << ", " << voxel_size[2] << std::endl;
      const set<Vertex*> vset = ROI_graph.vertices();
      set<Vertex*>::const_iterator roi_it, roi_it_end = vset.end();
      rc_ptr<BucketMap<Void> > roi_vertex_buckmap;
      short tex_roiLabel = 0;
      AimsVector<short int, 3> p;
      til::numeric_array<float, 3> pFloat;
      std::string roi_name;
    
      std::cout << "Iteration on: " <<std::flush;
      for(roi_it = vset.begin(); roi_it != roi_it_end; roi_it++)
      {
        tex_roiLabel += 1;
        const Vertex* roi_vertex = *roi_it;
        roi_vertex->getProperty("name", roi_name);
        std::cout << roi_name << ", roiLabel:" << tex_roiLabel << ", " << std::flush;
        if (roi_vertex->getProperty("aims_roi", roi_vertex_buckmap))
        {
          BucketMap<Void>::Bucket bck = roi_vertex_buckmap->begin()->second;
          BucketMap<Void>::Bucket::iterator itbck, itbck_end = bck.end();
//       Looking for closest points
          til::Find_closest< double, KDTree > fc(kdt);
          for (itbck = bck.begin(); itbck != itbck_end; itbck++)
          {
            p = itbck->first;
            til::numeric_array<float, 3> pFloat(p[0]*voxel_size[0], p[1]*voxel_size[1], p[2]*voxel_size[2]);
            std::cout << pFloat << ", " << std::flush;
            std::size_t pIndex_Mesh = fc(pFloat);
            if (til::dist2(pFloat, getVertices(mesh)[pIndex_Mesh], til::prec<float>()) <= distthresh)
            {
              tex[0][pIndex_Mesh] = tex_roiLabel;
            }
          }
        }
        else
        {
          std::cout << roi_name << " empty" << std::endl;
        }
      }
      std::cout << "..., OK." << std::endl;
      return tex_ptr;
    }
    else
    {
      std::cout << "error on ROI arg voxel_size" << std::endl;
      return 0;
    }
  }//fusionMeshRois_to_Graph*/

} // namespace constel

