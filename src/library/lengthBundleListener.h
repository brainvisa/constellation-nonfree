#ifndef CONSTELLATION_LENGTHBUNDLELISTENER_H
#define CONSTELLATION_LENGTHBUNDLELISTENER_H

#include <connectomist/fibertracking/bundles.h>

namespace constel
{

class LengthBundleListener : public comist::BundleListener
{
   public:

      LengthBundleListener(std::string & fibersLengthFileName);

      virtual void bundleStarted( const comist::BundleProducer &,
                                  const comist::BundleInfo & );
      virtual void bundleTerminated( const comist::BundleProducer &,
                                     const comist::BundleInfo & );
      virtual void fiberStarted( const comist::BundleProducer &,
                                 const comist::BundleInfo &,
                                 const comist::FiberInfo & );
      virtual void fiberTerminated( const comist::BundleProducer &,
                                    const comist::BundleInfo &,
                                    const comist::FiberInfo & );
      virtual void newFiberPoint( const comist::BundleProducer &,
                                  const comist::BundleInfo &,
                                  const comist::FiberInfo &, const comist::FiberPoint & );
      virtual void noMoreBundle( const comist::BundleProducer & );

      virtual ~LengthBundleListener();


   protected:
    float _fiberLength;//in mm
    comist::FiberPoint _antFiberPoint;
    uint _fiberPointCount;
    std::string _fileName;
    std::fstream _file;

};


inline void LengthBundleListener::bundleTerminated(
  const comist::BundleProducer &, const comist::BundleInfo & )
{
}

inline void LengthBundleListener::noMoreBundle( const comist::BundleProducer & )
{
}


} // namespace constel

#endif // ifndef CONSTELLATION_LENGTHBUNDLELISTENER_H
