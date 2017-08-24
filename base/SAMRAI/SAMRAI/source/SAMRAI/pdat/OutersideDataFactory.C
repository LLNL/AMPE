/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outerside data objects
 *
 ************************************************************************/

#ifndef included_pdat_OutersideDataFactory_C
#define included_pdat_OutersideDataFactory_C

#include "SAMRAI/pdat/OutersideDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/SideDataFactory.h"

#include "boost/make_shared.hpp"

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
OutersideDataFactory<TYPE>::OutersideDataFactory(
   const tbox::Dimension& dim,
   int depth):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim))
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OutersideDataFactory<TYPE>::~OutersideDataFactory()
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
OutersideDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   return boost::make_shared<OutersideDataFactory<TYPE> >(
             ghosts.getDim(),
             d_depth);
}

/*
 *************************************************************************
 *
 * Allocate the concrete outerside data classes.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchData>
OutersideDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   return boost::make_shared<OutersideData<TYPE> >(patch.getBox(), d_depth);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outerside data objects.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::BoxGeometry>
OutersideDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(getDim()));

   return boost::make_shared<OutersideGeometry>(box, zero_vector);
}

template<class TYPE>
int
OutersideDataFactory<TYPE>::getDepth() const
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
OutersideDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OutersideData<TYPE>));
   const size_t data = OutersideData<TYPE>::getSizeOfData(box, d_depth);
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
OutersideDataFactory<TYPE>::validCopyTo(
   const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are SideData and OutersideData.
    */
   if (!valid_copy) {
      boost::shared_ptr<SideDataFactory<TYPE> > sdf(
         boost::dynamic_pointer_cast<SideDataFactory<TYPE>,
                                     hier::PatchDataFactory>(dst_pdf));
      if (sdf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      boost::shared_ptr<OutersideDataFactory<TYPE> > osdf(
         boost::dynamic_pointer_cast<OutersideDataFactory<TYPE>,
                                     hier::PatchDataFactory>(
            dst_pdf));
      if (osdf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that fine data for the outerside
 * quantity will take precedence on coarse-fine interfaces.  See the
 * OutersideVariable<TYPE> class header file for more information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OutersideDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Return true since the outerside data index space extends beyond the
 * interior of patches.  That is, outerside data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OutersideDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return true;
}

}
}
#endif
