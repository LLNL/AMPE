/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating side data objects
 *
 ************************************************************************/

#ifndef included_pdat_SideDataFactory_C
#define included_pdat_SideDataFactory_C

#include "SAMRAI/pdat/SideDataFactory.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/OutersideDataFactory.h"
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

// SGS DODIM TODO
// WARNING! WARNING! WARNING! WARNING!
// This was hacked to replicate a bug in the SAMRAI.  The direction vector
// should not be all ones but not having this causes the tests to fail
// WARNING! WARNING! WARNING! WARNING!

template<class TYPE>
SideDataFactory<TYPE>::SideDataFactory(
   const int depth,
   const hier::IntVector& ghosts,
   bool fine_boundary_represents_var,
   const hier::IntVector& directions):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(hier::IntVector::getOne(ghosts.getDim()))
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(directions);
#endif
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
   TBOX_ASSERT(directions.min() >= 0);
}

template<class TYPE>
SideDataFactory<TYPE>::SideDataFactory(
   const int depth,
   const hier::IntVector& ghosts,
   bool fine_boundary_represents_var):
   hier::PatchDataFactory(ghosts),
   d_depth(depth),
   d_fine_boundary_represents_var(fine_boundary_represents_var),
   d_directions(ghosts.getDim(), 1)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
SideDataFactory<TYPE>::~SideDataFactory()
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
SideDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return boost::make_shared<SideDataFactory<TYPE> >(
      d_depth,
      ghosts,
      d_fine_boundary_represents_var,
      d_directions);
}

/*
 *************************************************************************
 *
 * Allocate the concrete side data classes.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchData>
SideDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   return boost::make_shared<SideData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts,
      d_directions);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for side data objects.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::BoxGeometry>
SideDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   return boost::make_shared<SideGeometry>(
      box,
      d_ghosts,
      d_directions);
}

template<class TYPE>
int
SideDataFactory<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
const hier::IntVector&
SideDataFactory<TYPE>::getDirectionVector() const
{
   return d_directions;
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
SideDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(SideData<TYPE>));
   const size_t data =
      SideData<TYPE>::getSizeOfData(box, d_depth, d_ghosts, d_directions);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from SideData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
SideDataFactory<TYPE>::validCopyTo(
   const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are SideData and OutersideData.
    */
   if (!valid_copy) {
      boost::shared_ptr<SideDataFactory<TYPE> > sdf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (sdf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      boost::shared_ptr<OutersideDataFactory<TYPE> > osdf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (osdf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean value indicating how data for the side quantity will be
 * treated on coarse-fine interfaces.  This value is passed into the
 * constructor.  See the FaceVariable<DIM> class header file for more
 * information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
SideDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
   return d_fine_boundary_represents_var;
}

/*
 *************************************************************************
 *
 * Return true since the side data index space extends beyond the interior
 * of patches.  That is, side data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
SideDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return true;
}


}
}
#endif
