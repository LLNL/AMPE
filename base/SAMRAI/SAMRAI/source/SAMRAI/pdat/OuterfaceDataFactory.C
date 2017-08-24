/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating outerface data objects
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceDataFactory_C
#define included_pdat_OuterfaceDataFactory_C

#include "SAMRAI/pdat/OuterfaceDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OuterfaceGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/FaceDataFactory.h"

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
OuterfaceDataFactory<TYPE>::OuterfaceDataFactory(
   const tbox::Dimension& dim,
   int depth):
   hier::PatchDataFactory(hier::IntVector::getZero(dim)),
   d_depth(depth),
   d_no_ghosts(hier::IntVector::getZero(dim))
{
   TBOX_ASSERT(depth > 0);
}

template<class TYPE>
OuterfaceDataFactory<TYPE>::~OuterfaceDataFactory()
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
OuterfaceDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, ghosts);

   return boost::make_shared<OuterfaceDataFactory<TYPE> >(
             ghosts.getDim(),
             d_depth);
}

/*
 *************************************************************************
 *
 * Allocate the concrete outerface data classes.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchData>
OuterfaceDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, patch);

   return boost::make_shared<OuterfaceData<TYPE> >(patch.getBox(), d_depth);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for outerface data objects.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::BoxGeometry>
OuterfaceDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const hier::IntVector zero_vector(hier::IntVector::getZero(getDim()));

   return boost::make_shared<OuterfaceGeometry>(box, zero_vector);
}

template<class TYPE>
int
OuterfaceDataFactory<TYPE>::getDepth() const
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
OuterfaceDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, box);

   const size_t obj = tbox::MemoryUtilities::align(sizeof(OuterfaceData<TYPE>));
   const size_t data = OuterfaceData<TYPE>::getSizeOfData(box, d_depth);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that fine data for the outerface
 * quantity will take precedence on coarse-fine interfaces.  See the
 * OuterfaceVariable<TYPE> class header file for more information.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuterfaceDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Return true since the outerface data index space extends beyond the
 * interior of patches.  That is, outerface data lives on patch borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
OuterfaceDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from OuterfaceData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
OuterfaceDataFactory<TYPE>::validCopyTo(
   const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Valid options are FaceData and OuterfaceData.
    */
   if (!valid_copy) {
      boost::shared_ptr<FaceDataFactory<TYPE> > fdf(
         boost::dynamic_pointer_cast<FaceDataFactory<TYPE>,
                                     hier::PatchDataFactory>(dst_pdf));
      if (fdf) {
         valid_copy = true;
      }
   }

   if (!valid_copy) {
      boost::shared_ptr<OuterfaceDataFactory<TYPE> > ofdf(
         boost::dynamic_pointer_cast<OuterfaceDataFactory<TYPE>,
                                     hier::PatchDataFactory>(
            dst_pdf));
      if (ofdf) {
         valid_copy = true;
      }
   }

   return valid_copy;
}

}
}
#endif
