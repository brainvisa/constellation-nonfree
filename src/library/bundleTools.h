#ifndef CONSTELLATION_BUNDLETOOLS_H
#define CONSTELLATION_BUNDLETOOLS_H

#include <constellation/bundleSet.h>
#include <constellation/sparseMatrix.h>
#include <constellation/tildefs.h>
#include <constellation/connectivities.h>


  //-------------//
 //  ListenedFiberInfo  //
//-------------//

namespace constel
{

class MeshIntersectionBundleListener;


class ListenedFiberInfo
{
public:
  
  ListenedFiberInfo();
  ListenedFiberInfo( int id );
  int id() const;
  void setCurvilinearAbscissa(float value);
  float getCurvilinearAbscissa() const {return _curvilinearAbscissa;};
  void setAntFiberPoint(comist::FiberPoint antFiberPoint);
  comist::FiberPoint getAntFiberPoint() const {return _antFiberPoint;};
  void pushBackMeshIntersectionNeighbourhood(constel::QuickMap fiberIntersectionNeighDistMap, int meshIntersectionMeshId);
  std::vector<constel::QuickMap > getFiberIntersectionNeighDistMapVector() const {return _fiberIntersectionNeighDistMapVector;};
  std::vector<float> getFiberMeshIntersectionCurvilinearAbscissaVector() const {return _fiberMeshIntersectionCurvilinearAbscissaVector;};
  std::vector<int> getFiberMeshIntersectionMeshIdentityVector() const {return _fiberMeshIntersectionMeshIdentityVector;};
  void clearFiberMeshIntersectionInfo();
  void setAntFiberPointExistingMeshIntersection(bool intersect);
  bool getAntFiberPointExistingMeshIntersection() const { return _antFiberPoint_ExistingMeshIntersection; };
  void setAntFiberPointMeshClosestPointIndex(std::size_t meshVertex_index);
  std::size_t getAntFiberPointMeshClosestPointIndex() const { return _antFiberPointMeshClosestPoint_index; };
//   void savingFiberMeshIntersectionPoints(std::fstream file);

  
protected:

  int _id;
  float _curvilinearAbscissa;
  comist::FiberPoint _antFiberPoint;
  bool _antFiberPoint_ExistingMeshIntersection;
  std::size_t _antFiberPointMeshClosestPoint_index;
  std::vector<constel::QuickMap > _fiberIntersectionNeighDistMapVector;
  std::vector<float> _fiberMeshIntersectionCurvilinearAbscissaVector;
  std::vector<int> _fiberMeshIntersectionMeshIdentityVector;
  std::vector<std::size_t> _fiberMeshIntersectionNbVector;//nb of intersection per mesh: _fiberMeshIntersectionNb[i] = fiber intersection nb with mesh i (ith mesh added to the BundleInteractionReader by addBundleListener
  

};

  //-------------//
 //  BundleInteractionReader  //
//-------------//

class BundleInteractionReader : public comist::BundleReader
{
  public:
    BundleInteractionReader( const std::string &fileName );
//     void addBundleListener( MeshIntersectionBundleListener & meshIntersectionBundleListener);
    virtual ~BundleInteractionReader();
  
  protected:

    friend class MemAntBundleListener;
    friend class AfficheAntFiberPointBundleListener;
    friend class CurvilinearAbscissaBundleListener;
    friend class MeshIntersectionBundleListener;
    friend class SavingMeshIntersectionBundleListener;
    friend class MeshConnectionBundleListener;
    friend class MeshHistoLengthConnectionBundleListener;
    friend class FiberNameAccordingToMeshIntersectionBundleListener;
    friend class SubSamplerFromMeshIntersectionBundleListener;
    friend class FiberNameAccordingToMeshTextureIntersectionBundleListener;
    friend class MeshIntersectionNoSmoothingBundleListener;
    friend class MeshIntersectionNoSmoothingFasterBundleListener;
    friend class MeshIntersectionMatrixBundleListener;
    friend class MeshIntersectionMatrixWithLengthBundleListener;

  protected:
    ListenedFiberInfo _listenedFiberInfo;
    int _meshIntersectionBundleListener_nb;
  
};

  //-------------//
 //  MemAntBundleListener  //
//-------------//
//always at the end of the BundleListenerList (of the associated BundleInteractionReader

class MemAntBundleListener : public comist::BundleListener
{
   public:

      MemAntBundleListener(BundleInteractionReader &bundleInteractionReader );

      virtual void newFiberPoint( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &, const comist::FiberPoint & );

      virtual ~MemAntBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;
    

};


 //-------------//
 //  AfficheAntFiberPointBundleListener  //
//-------------//
class AfficheAntFiberPointBundleListener : public comist::BundleListener
{
   public:

      AfficheAntFiberPointBundleListener(BundleInteractionReader &bundleInteractionReader );

      virtual void newFiberPoint( const comist::BundleProducer &, const comist::BundleInfo &,
                                  const comist::FiberInfo &, const comist::FiberPoint & );

      virtual ~AfficheAntFiberPointBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;
    

};

//-------------//
 //  CurvilinearAbscissaBundleListener  //
//-------------//
class CurvilinearAbscissaBundleListener : public comist::BundleListener
{
   public:

      CurvilinearAbscissaBundleListener(BundleInteractionReader &bundleInteractionReader );

      virtual void fiberStarted( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void newFiberPoint( const comist::BundleProducer &, const comist::BundleInfo &,
                                  const comist::FiberInfo &, const comist::FiberPoint & );

      virtual ~CurvilinearAbscissaBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;

   protected:
    unsigned _fiberPointCount;
    float _fiberLength;
    

};

//-------------//
 //  MeshIntersectionBundleListener  //
//-------------//

class MeshIntersectionBundleListener : public comist::BundleListener
{
    
    public:

      MeshIntersectionBundleListener(const AimsSurfaceTriangle & aimsMesh,  BundleInteractionReader &bundleInteractionReader, double meshDistanceThreshold = 0., double meshClosestPointMaxDistance = 1., bool verbose = false );

//       virtual void bundleTerminated( const comist::BundleProducer &, const comist::BundleInfo & );
      virtual void fiberStarted( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void newFiberPoint( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &, const comist::FiberPoint & );
      virtual void noMoreBundle( const comist::BundleProducer & );
      void setMeshIdentity(int meshIdentity);

      virtual ~MeshIntersectionBundleListener();

    protected:
      BundleInteractionReader * _bundleInteractionReader;
      double _meshDistanceThreshold;
      double _meshClosestPointMaxDistance;
      int _meshIdentity;
      bool _verbose;
      const AimsSurfaceTriangle & _aimsMesh;

      PolygonsByVertexIndex _meshPolygonsByVertex_Index;
      constel::Mesh _mesh;
      //For computing closest point to meshes in a fast computation time:
      constel::KDTree * _mesh_kdt_ptr;
      til::Find_closest< double, constel::KDTree > * _mesh_fc_ptr;
      std::vector<constel::QuickMap> _meshDistanceThresholdNeighborhoodByVertex;

      //For intersection computing:
      bool _antFiberPoint_ExistingMeshIntersection;
      std::size_t _antFiberPointMeshClosestPoint_index;
      float _antFiberPointMeshClosestPoint_dist;
      //variables for count:
      std::size_t _fiberPointCount;
};

//-------------//
 //  MeshIntersectionNoSmoothingBundleListener  //
//-------------//

class MeshIntersectionNoSmoothingBundleListener : public comist::BundleListener
{
    
    public:

      MeshIntersectionNoSmoothingBundleListener(const AimsSurfaceTriangle & aimsMesh,  BundleInteractionReader &bundleInteractionReader, double meshClosestPointMaxDistance = 7., bool verbose = false );

//       virtual void bundleTerminated( const comist::BundleProducer &, const comist::BundleInfo & );
      virtual void fiberStarted( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void newFiberPoint( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &, const comist::FiberPoint & );
      virtual void noMoreBundle( const comist::BundleProducer & );
      void setMeshIdentity(int meshIdentity);

      virtual ~MeshIntersectionNoSmoothingBundleListener();

    protected:
      BundleInteractionReader * _bundleInteractionReader;
      double _meshClosestPointMaxDistance;
      bool _verbose;
      int _meshIdentity;
      const AimsSurfaceTriangle & _aimsMesh;

      PolygonsByVertexIndex _meshPolygonsByVertex_Index;
      constel::Mesh _mesh;
      //For computing closest point to meshes in a fast computation time:
      constel::KDTree * _mesh_kdt_ptr;
      til::Find_closest< double, constel::KDTree > * _mesh_fc_ptr;

      //For intersection computing:
      bool _antFiberPoint_ExistingMeshIntersection;
      std::size_t _antFiberPointMeshClosestPoint_index;
      float _antFiberPointMeshClosestPoint_dist;
      //variables for count:
      std::size_t _fiberPointCount;
      std::size_t _fiberCount;
};//class MeshIntersectionNoSmoothingBundleListener

//-------------//
 //  MeshIntersectionNoSmoothingFasterBundleListener  //
// idem MeshIntersectionNoSmoothingBundleListener but before looking for an intersection, test if the fiberPoints are in the input rois mask (ribbon around the input Mesh, given by _roisMask)
//-------------//

class MeshIntersectionNoSmoothingFasterBundleListener : public MeshIntersectionNoSmoothingBundleListener
{
    
    public:

      MeshIntersectionNoSmoothingFasterBundleListener(const AimsSurfaceTriangle & aimsMesh, AimsData<short> &roisMask, BundleInteractionReader &bundleInteractionReader, double meshClosestPointMaxDistance = 7., bool verbose = false );

//       virtual void bundleTerminated( const comist::BundleProducer &, const comist::BundleInfo & );
      virtual void fiberStarted( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void newFiberPoint( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &, const comist::FiberPoint & );
      virtual void noMoreBundle( const comist::BundleProducer & );
      void setMeshIdentity(int meshIdentity);

      virtual ~MeshIntersectionNoSmoothingFasterBundleListener();

    protected:
      AimsData<short>& _roisMask;
      bool _antFiberPoint_inRoisMask;
      
};//class MeshIntersectionNoSmoothingFasterBundleListener

//-------------//
 //  SavingMeshIntersectionBundleListener //
//-------------//
class SavingMeshIntersectionBundleListener : public comist::BundleListener
{
   public:

      SavingMeshIntersectionBundleListener(BundleInteractionReader &bundleInteractionReader, std::string fileName );

      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );

      virtual ~SavingMeshIntersectionBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;

   protected:
    std::string _fileName;//fiber names according to intersection points
    std::fstream _file;
    

};
//-------------//
 //  MeshIntersectionMatrixBundleListener //
//for one mesh first, of meshId = 0
//-------------//
class MeshIntersectionMatrixBundleListener : public comist::BundleListener
{
   public:

      MeshIntersectionMatrixBundleListener(BundleInteractionReader &bundleInteractionReader, int meshIdentity, unsigned meshVertexNb, std::string file_name = "", std::string matrixRowsBundleNames_file_name = "", int bundlesNb = 1 );

      virtual void bundleStarted( const comist::BundleProducer &, const comist::BundleInfo & );
      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void noMoreBundle( const comist::BundleProducer & bundleProducer);

      virtual ~MeshIntersectionMatrixBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;
    int _meshIdentity;
    unsigned _meshVertexNb;
    aims::SparseMatrix _meshIntersectionsMatrix;
    std::string _file_name;
    std::string _matrixRowsBundleNames_file_name;
    std::fstream _matrixRowsBundleNames_file;
    int _bundlesNb;
    int _bundlesCount;

};//class MeshIntersectionMatrixBundleListener

//-------------//
 //  MeshIntersectionMatrixWithLengthBundleListener //
//for one mesh first, of meshId = 0, and one bundles per region (basal ganglia for exemple)
// idem MeshIntersectionMatrixBundleListener, but each connectivity value is weighted with the fiber tract length
//-------------//
class MeshIntersectionMatrixWithLengthBundleListener : public comist::BundleListener
{
   public:

      MeshIntersectionMatrixWithLengthBundleListener(BundleInteractionReader &bundleInteractionReader, int meshIdentity, unsigned meshVertexNb, std::string file_name = "", std::string matrixRowsBundleNames_file_name = "", int bundlesNb = 1 );

      virtual void bundleStarted( const comist::BundleProducer &, const comist::BundleInfo & );
      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void noMoreBundle( const comist::BundleProducer & bundleProducer);

      virtual ~MeshIntersectionMatrixWithLengthBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;
    int _meshIdentity;
    unsigned _meshVertexNb;
    aims::SparseMatrix _meshIntersectionsMatrix;
    std::string _file_name;
    std::string _matrixRowsBundleNames_file_name;
    std::fstream _matrixRowsBundleNames_file;
    int _bundlesNb;
    int _bundlesCount;

};//class MeshIntersectionMatrixWithLengthBundleListener

//-------------//
 //  MeshConnectionBundleListener //
//for one mesh first, of meshId = 0
//-------------//
class MeshConnectionBundleListener : public comist::BundleListener
{
   public:

      MeshConnectionBundleListener(BundleInteractionReader &bundleInteractionReader, int meshIdentity, bool verbose = false );

      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );

      boost::shared_ptr<constel::BundleConnections> getBundleMeshConnections() const {return _bundleMeshConnections; }
      boost::shared_ptr<constel::ConnectionsLength> getBundleMeshConnectionsLength()  const  { return _bundleMeshConnectionsLength; }

      virtual ~MeshConnectionBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;
    int _meshIdentity;
    bool _verbose;
    boost::shared_ptr<constel::BundleConnections> _bundleMeshConnections;
    boost::shared_ptr<constel::ConnectionsLength> _bundleMeshConnectionsLength;
    unsigned _bundleMeshConnectionsCount;

};

//-------------//
 //  MeshHistoLengthConnectionBundleListener //
//for one mesh first, of meshId = 0
//-------------//
class MeshHistoLengthConnectionBundleListener : public comist::BundleListener
{
   public:

      MeshHistoLengthConnectionBundleListener(BundleInteractionReader &bundleInteractionReader, int meshIdentity, unsigned meshVertexNb, int connectionLengthMin, int connectionLengthMax, double meshDistanceThreshold = 0. );

      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );

      boost::shared_ptr<AimsData< float > > getBundleMeshConnectionsHistoLength() const {return _meshConnectionsHistoLength_ptr; }

      virtual ~MeshHistoLengthConnectionBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;
    int _meshIdentity;
    unsigned _bundleMeshConnectionsCount;
    unsigned _meshVertexNb;

    //For histogram creation:
    boost::shared_ptr<AimsData< float > > _meshConnectionsHistoLength_ptr;
    AimsData< float >  & _meshConnectionsHistoLength;
    int _connectionLengthMin;
    int _connectionLengthMax;//correspond to mm, corresponding to a python histogram, with bins = range(_connectionLengthMin,_connectionLengthMax+1)

    //For smoothing:
    double _meshDistanceThreshold;
    double _two_pi;
    double _square_sigma;

};

//-------------//
 //  FiberNameAccordingToMeshIntersectionBundleListener //
//-------------//
class FiberNameAccordingToMeshIntersectionBundleListener : public comist::BundleListener
{
   public:

      FiberNameAccordingToMeshIntersectionBundleListener(BundleInteractionReader &bundleInteractionReader, std::string fileName );

      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual ~FiberNameAccordingToMeshIntersectionBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;

   protected:
    std::string _fileName;//fiber names according to intersection points
    std::fstream _file;
    

};

//-------------//
 //  FiberNameAccordingToMeshTextureIntersectionBundleListener //
//-------------//
class FiberNameAccordingToMeshTextureIntersectionBundleListener : public comist::BundleListener
{
   public:

      FiberNameAccordingToMeshTextureIntersectionBundleListener(BundleInteractionReader &bundleInteractionReader, const TimeTexture<short> & labeled_tex, std::string fileName );
      virtual void fiberStarted( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void newFiberPoint( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo &, const comist::FiberPoint & );
      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual ~FiberNameAccordingToMeshTextureIntersectionBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;

   protected:
    std::string _fileName;//fiber names according to intersection points
    std::fstream _file;
    const Texture<short> & _labeledTex0;
    std::set<short> _intersectionPointsMeshClosestPointsTexLabels_set;
    

};


//-------------//
 //  SubSamplerFromMeshIntersectionBundleListener //
//-------------//
class SubSamplerFromMeshIntersectionBundleListener : public comist::BundleProducer, public comist::BundleListener
{
   public:

      SubSamplerFromMeshIntersectionBundleListener(BundleInteractionReader &bundleInteractionReader, int step);
      virtual void bundleStarted( const comist::BundleProducer & bundleProducer, const comist::BundleInfo & bundleInfo);
      virtual void bundleTerminated( const comist::BundleProducer &, const comist::BundleInfo & );
      virtual void fiberStarted( const comist::BundleProducer & bundleProducer, const comist::BundleInfo & bundleInfo, const comist::FiberInfo & fiberInfo);
      virtual void newFiberPoint( const comist::BundleProducer & bundleProducer, const comist::BundleInfo & bundleInfo, const comist::FiberInfo & fiberInfo, const comist::FiberPoint &point );
      virtual void fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & );
      virtual void noMoreBundle( const comist::BundleProducer & bundleProducer);
      virtual ~SubSamplerFromMeshIntersectionBundleListener();


   private:
    BundleInteractionReader * _bundleInteractionReader;

   protected:
    Fiber _fiber;
    Fibers _faisceau;
    BundlesSet _bundlesSet;
    std::vector<std::string> _names;
    int _step;
    std::size_t _fiberPointCount;    

};

} //namespace constel
#endif // ifndef CONSTELLATION_BUNDLETOOLS_H

