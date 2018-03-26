/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating edge data objects
 *
 ************************************************************************/

#ifndef included_pdat_EdgeDataFactory
#define included_pdat_EdgeDataFactory

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/tbox/Complex.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace pdat {

/**
 * Class EdgeDataFactory is a factory class used to allocate new
 * instances of EdgeData objects.  It is a subclass of the patch
 * data factory class and edge data is a subclass of patch data.  Both
 * the factory and data classes are templated on the type of the contained
 * object (e.g., double or int).
 *
 * @see pdat::EdgeData
 * @see pdat::PatchDataFactory
 */

template<class TYPE>
class EdgeDataFactory:public hier::PatchDataFactory
{
public:
   /**
    * The constructor for the edge data factory class.  The ghost cell width, depth
    * (number of components), and fine boundary representation arguments give the
    * defaults for all edge data objects created with this factory.  See
    * the EdgeVariable<DIM> class header file for more information.
    */
   EdgeDataFactory(
      int depth,
      const hier::IntVector& ghosts,
      bool fine_boundary_represents_var);

   /**
    * Virtual destructor for the edge data factory class.
    */
   virtual ~EdgeDataFactory<TYPE>();

   /**
    * @brief Abstract virtual function to clone a patch data factory.
    *
    * This will return a new instantiation of the abstract factory
    * with the same properties.  The properties of the cloned factory
    * can then be changed without modifying the original.
    *
    * @param ghosts default ghost cell width for concrete classes created from
    * the factory.
    */
   virtual boost::shared_ptr<hier::PatchDataFactory>
   cloneFactory(
      const hier::IntVector& ghosts);

   /**
    * Virtual factory function to allocate a concrete edge data object.
    * The default information about the object (e.g., ghost cell width)
    * is taken from the factory.
    */

   virtual boost::shared_ptr<hier::PatchData>
   allocate(
      const hier::Patch& patch) const;

   /**
    * Allocate the box geometry object associated with the patch data.
    * This information will be used in the computation of intersections
    * and data dependencies between objects.
    */
   virtual boost::shared_ptr<hier::BoxGeometry>
   getBoxGeometry(
      const hier::Box& box) const;

   /**
    * Get the depth (number of components).  This is the depth that
    * will be used in the instantiation of edge data objects.
    */
   int
   getDepth() const;

   /**
    * Calculate the amount of memory needed to store the edge data object,
    * including object data and dynamically allocated data.
    */
   virtual size_t
   getSizeOfMemory(
      const hier::Box& box) const;

   /**
    * Return a boolean value indicating how data for the edge quantity will be
    * treated on coarse-fine interfaces.  This value is passed into the
    * constructor.  See the EdgeVariable<DIM> class header file for more
    * information.
    */
   bool
   fineBoundaryRepresentsVariable() const;

   /**
    * Return true since the edge data index space extends beyond the interior
    * of patches.  That is, edge data lives on patch borders.
    */
   bool
   dataLivesOnPatchBorder() const;

   /**
    * Return whether it is valid to copy this EdgeDataFactory to the
    * supplied destination patch data factory.  It will return true if
    * dst_pdf is EdgeDataFactory or OuteredgeDataFactory, false otherwise.
    */
   bool
   validCopyTo(
      const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const;

private:
   int d_depth;
   bool d_fine_boundary_represents_var;

};

}
}

#include "SAMRAI/pdat/EdgeDataFactory.C"

#endif
