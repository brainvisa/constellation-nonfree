/* Copyright (c) 1995-2012 CEA
 *
 *  This software and supporting documentation were developed by
 *      CEA/DSV/SHFJ
 *      4 place du General Leclerc
 *      91401 Orsay cedex
 *      France
 *
 * This software is governed by the CeCILL license version 2 under
 * French law and abiding by the rules of distribution of free software.
 * You can  use, modify and/or redistribute the software under the
 * terms of the CeCILL license version 2 as circulated by CEA, CNRS
 * and INRIA at the following URL "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license version 2 and that you accept its terms.
 */

#ifndef CONSTELLATION_BUNDLELOADER_H
#define CONSTELLATION_BUNDLELOADER_H

/* Extracted from code from Pascal Cathier
 */

#include <connectomist/fibertracking/bundles.h>

namespace constel
{

  class BundleLoader : public comist::BundleListener
  {
  public: // typedefs

    typedef comist::FiberPoint Point;
    typedef std::vector<Point> Fiber;
    typedef std::vector<Fiber> Bundle;

  public: // constructors & destructor

    BundleLoader() : m_fibers(new Bundle) {}
    virtual ~BundleLoader() {}

  public: // set & get

    carto::rc_ptr<Bundle> getFibers() const { return m_fibers; }

  private: // functions

    void bundleStarted( const comist::BundleProducer &,
                        const comist::BundleInfo & )
    {
    }

    void bundleTerminated( const comist::BundleProducer &,
                           const comist::BundleInfo & )
    {
    }

    void fiberStarted( const comist::BundleProducer &,
                       const comist::BundleInfo &, const comist::FiberInfo & )
    {
      m_fibers->push_back(Fiber());
    }

    void fiberTerminated( const comist::BundleProducer &,
                          const comist::BundleInfo &,
                          const comist::FiberInfo & )
    {
    }

    void newFiberPoint( const comist::BundleProducer &,
                        const comist::BundleInfo &, const comist::FiberInfo &,
                        const comist::FiberPoint & point )
    {
      m_fibers->back().push_back(point);
    }

    void noMoreBundle( const comist::BundleProducer & )
    {
    }

  private: // data

    carto::rc_ptr<Bundle> m_fibers;
  };

}

#endif


