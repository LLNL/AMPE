/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outeredge data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuteredgeDataFactory_C
#define included_pdat_OuteredgeDataFactory_C

#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/pdat/OuteredgeDataFactory.h"
#include "SAMRAI/pdat/OuteredgeGeometry.h"
#include "SAMRAI/hier/Patch.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace pdat {

/*
 *************************************************************************
 *
 * The constructor simply caches the depth of the patch data.
 *
 *************************************************************************
 */

template<class TYPE>
OuteredgeDataFactory<TYPE>::OuteredgeDataFactory(
   const tbox::Dimension& dim,
   int depth):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim))
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OuteredgeDataFactory<TYPE>::~OuteredgeDataFactory()
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
OuteredgeDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return boost::make_shared<OuteredgeDataFactory<TYPE> >(
      ghosts.getDim(),
      d_depth);
}

/*
 *************************************************************************
 *
 * Allocate the concrete outeredge data classes.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchData>
OuteredgeDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   return boost::make_shared<OuteredgeData<TYPE> >(patch.getBox(), d_depth);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outeredge data objects.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::BoxGeometry>
OuteredgeDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));

   return boost::make_shared<OuteredgeGeometry>(box, zero_vector);
}

template<class TYPE>
int
OuteredgeDataFactory<TYPE>::getDepth() const
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
OuteredgeDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OuteredgeData<TYPE>));
   const size_t data = OuteredgeData<TYPE>::getSizeOfData(box,
         d_depth);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that fine data for the outeredge
 * quantity will take precedence on coarse-fine interfaces.  See the
 * OuteredgeVariable<DIM> class header file for more information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuteredgeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const {
   return true;
}

/*
 *************************************************************************
 *
 * Return true since the outeredge data index space extends beyond the
 * interior of patches.  That is, outeredge data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuteredgeDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from EdgeData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
OuteredgeDataFactory<TYPE>::validCopyTo(
   const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are EdgeData and OuteredgeData.
    */
   if (!valid_copy) {
      boost::shared_ptr<EdgeDataFactory<TYPE> > edf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (edf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      boost::shared_ptr<OuteredgeDataFactory<TYPE> > oedf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());
      if (oedf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
