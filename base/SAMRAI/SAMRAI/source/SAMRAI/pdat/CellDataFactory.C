/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory class for creating cell data objects
 *
 ************************************************************************/

#ifndef included_pdat_CellDataFactory_C
#define included_pdat_CellDataFactory_C

#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/hier/Patch.h"

#include <boost/make_shared.hpp>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

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
CellDataFactory<TYPE>::CellDataFactory(
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchDataFactory(ghosts),
   d_depth(depth)
{
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
}

template<class TYPE>
CellDataFactory<TYPE>::~CellDataFactory()
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
CellDataFactory<TYPE>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);

   return boost::make_shared<CellDataFactory<TYPE> >(d_depth, ghosts);
}

/*
 *************************************************************************
 *
 * Allocate the concrete cell data classes.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::PatchData>
CellDataFactory<TYPE>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);

   return boost::make_shared<CellData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

/*
 *************************************************************************
 *
 * Return the box geometry type for cell data objects.
 *
 *************************************************************************
 */

template<class TYPE>
boost::shared_ptr<hier::BoxGeometry>
CellDataFactory<TYPE>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   return boost::make_shared<CellGeometry>(box, d_ghosts);
}

template<class TYPE>
int
CellDataFactory<TYPE>::getDepth() const
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
CellDataFactory<TYPE>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   const size_t obj =
      tbox::MemoryUtilities::align(sizeof(CellData<TYPE>));
   const size_t data =
      CellData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);
   return obj + data;
}

/*
 *************************************************************************
 *
 * Determine whether this is a valid copy operation to/from CellData
 * between the supplied datatype.
 *
 *************************************************************************
 */

template<class TYPE>
bool
CellDataFactory<TYPE>::validCopyTo(
   const boost::shared_ptr<hier::PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);

   bool valid_copy = false;

   /*
    * Only valid option is CellData.
    */
   boost::shared_ptr<CellDataFactory<TYPE> > cdf(
      dst_pdf,
      boost::detail::dynamic_cast_tag());
   if (cdf) {
      valid_copy = true;
   }
   return valid_copy;
}

/*
 *************************************************************************
 *
 * Return a boolean true value indicating that the cell data quantities will
 * always be treated as though fine values represent them on coarse-fine
 * interfaces.
 *
 *************************************************************************
 */
template<class TYPE>
bool
CellDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
   return true;
}

/*
 *************************************************************************
 *
 * Return false since the cell data index space matches the cell-centered
 * index space for AMR patches.  Thus, cell data does not live on patch
 * borders.
 *
 *************************************************************************
 */
template<class TYPE>
bool
CellDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
   return false;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
