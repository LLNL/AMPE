/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating node data objects
 *
 ************************************************************************/

#ifndef included_pdat_NodeDataFactory_C
#define included_pdat_NodeDataFactory_C

#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/OuternodeDataFactory.h"
#include "SAMRAI/hier/Patch.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the default ghost cell width and depth.
 *
 *************************************************************************
 */

template<class TYPE>
NodeDataFactory<TYPE>::NodeDataFactory(
   int depth,
   const hier::IntVector& ghosts,
   bool fine_boundary_represents_var):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
NodeDataFactory<TYPE>::~NodeDataFactory()
{
}

/*
 *************************************************************************
 *
 * Clone the factory and copy the default parameters to the new factory.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchDataFactory>
NodeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return boost::make_shared<NodeDataFactory<TYPE> >(
      d_depth,
      ghosts,
      d_fine_boundary_represents_var);
}

/*
 *************************************************************************
 *
 * Allocate the concrete node data classes.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchData>
NodeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   return boost::make_shared<NodeData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for node data objects.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::BoxGeometry>
NodeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   return boost::make_shared<NodeGeometry>(box, d_ghosts);
}

template<class TYPE>
int
NodeDataFactory<TYPE>::getDepth() const
{
   return d_depth;
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory needed to allocate the data object.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
NodeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(NodeData<TYPE>));
   const size_t data =
      NodeData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from NodeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
NodeDataFactory<TYPE>::validCopyTo(
   const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are NodeData and OuternodeData.
    */
   if (!valid_copy) {
      boost::shared_ptr<NodeDataFactory<TYPE> > ndf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (ndf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      boost::shared_ptr<OuternodeDataFactory<TYPE> > ondf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (ondf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean value indicating how data for the node quantity will be
 * treated on coarse-fine interfaces.  This value is passed into the
 * constructor.  See the NodeVariable<DIM> class header file for more
 * information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
NodeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const {
   return d_fine_boundary_represents_var;
}

/*
 *************************************************************************
 *
 * Return true since the node data index space extends beyond the interior
 * of patches.  That is, node data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
NodeDataFactory<TYPE>::dataLivesOnPatchBorder() const {
   return true;
}

}
}
#endif
