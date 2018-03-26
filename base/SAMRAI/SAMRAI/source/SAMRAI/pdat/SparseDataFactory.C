/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Implementation for SparseDataFactory
 *
 ************************************************************************/
#ifndef included_pdat_SparseDataFactory_C
#define included_pdat_SparseDataFactory_C

#include "SAMRAI/pdat/SparseDataFactory.h"

#ifdef HAVE_BOOST_HEADERS

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/SparseData.h"
#include "SAMRAI/tbox/MemoryUtilities.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace pdat {

/*
 * C'tors and d'tors.
 */
template<typename BOX_GEOMETRY>
SparseDataFactory<BOX_GEOMETRY>::SparseDataFactory(
   const hier::IntVector& ghosts,
   const std::vector<std::string>& dbl_attributes,
   const std::vector<std::string>& int_attributes):
   hier::PatchDataFactory(ghosts),
   d_dbl_attributes(dbl_attributes),
   d_int_attributes(int_attributes)
{
}

template<typename BOX_GEOMETRY>
SparseDataFactory<BOX_GEOMETRY>::~SparseDataFactory()
{
}

/*
 * Implementation of base class pure virtual functions
 */
template<typename BOX_GEOMETRY>
boost::shared_ptr<hier::PatchDataFactory>
SparseDataFactory<BOX_GEOMETRY>::cloneFactory(
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);
   return boost::make_shared<SparseDataFactory<BOX_GEOMETRY> >(
      ghosts,
      d_dbl_attributes,
      d_int_attributes);
}

template<typename BOX_GEOMETRY>
boost::shared_ptr<hier::PatchData>
SparseDataFactory<BOX_GEOMETRY>::allocate(
   const hier::Patch& patch) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, patch);
   return boost::make_shared<SparseData<BOX_GEOMETRY> >(
         patch.getBox(),
         d_ghosts,
         d_dbl_attributes,
         d_int_attributes);
}

template<typename BOX_GEOMETRY>
boost::shared_ptr<hier::BoxGeometry>
SparseDataFactory<BOX_GEOMETRY>::getBoxGeometry(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   return boost::make_shared<BOX_GEOMETRY>(box, d_ghosts);
}

template<typename BOX_GEOMETRY>
size_t
SparseDataFactory<BOX_GEOMETRY>::getSizeOfMemory(
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   NULL_USE(box);
   return tbox::MemoryUtilities::align(
      sizeof(SparseData<BOX_GEOMETRY>));
}

template<typename BOX_GEOMETRY>
bool
SparseDataFactory<BOX_GEOMETRY>::validCopyTo(
   const boost::shared_ptr<PatchDataFactory>& dst_pdf) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_pdf);
   bool valid_copy = false;

   if (!valid_copy) {

      boost::shared_ptr<SparseDataFactory<BOX_GEOMETRY> > idf(
         dst_pdf,
         boost::detail::dynamic_cast_tag());

      if (idf) {
         valid_copy = true;
      }
   }
   return valid_copy;
}

} // end namespace pdat
} // end namespace SAMRAI

#endif
#endif
