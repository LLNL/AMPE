/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeIndex_C
#define included_pdat_NodeIndex_C

#include "SAMRAI/pdat/NodeIndex.h"

namespace SAMRAI {
namespace pdat {

std::vector<hier::IntVector> NodeIndex::s_offsets[tbox::Dimension::
                                                  MAXIMUM_DIMENSION_VALUE];
bool NodeIndex::s_offsets_are_set[tbox::Dimension::MAXIMUM_DIMENSION_VALUE] = { false };

NodeIndex::NodeIndex(
   const tbox::Dimension& dim):
   hier::Index(dim)
{
   setOffsets();
}

NodeIndex::NodeIndex(
   const hier::Index& rhs,
   const Corner corner):
   hier::Index(rhs.getDim())
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(* this, rhs);

   setOffsets();
   hier::IntVector::operator = (
      rhs + s_offsets[getDim().getValue() - 1][(int)corner]);
}

NodeIndex::NodeIndex(
   const hier::Index& rhs,
   const hier::IntVector& corner):
   hier::Index(rhs.getDim())
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < getDim().getValue(); i++) {
      TBOX_ASSERT(corner(i) == 0 || corner(i) == 1);
   }
#endif
   setOffsets();
   hier::IntVector::operator = (
      rhs + corner);
}

NodeIndex::NodeIndex(
   const NodeIndex& rhs):
   hier::Index(rhs)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);

   setOffsets();
}

NodeIndex::~NodeIndex()
{
}

void
NodeIndex::setOffsets()
{
   const tbox::Dimension& dim(getDim());
   int dim_index = dim.getValue() - 1;
   if (!s_offsets_are_set[dim_index]) {
      s_offsets[dim_index] = std::vector<hier::IntVector>(
            2 << tbox::Dimension::MAXIMUM_DIMENSION_VALUE,
            hier::IntVector(dim));
      for (int i = 0; i < (1 << dim.getValue()); i++) {
         hier::IntVector offset(dim, 0);

         offset(0) = i % 2;
         for (int j = 1; j < dim.getValue(); j++) {
            offset(j) = (i / (1 << j)) % 2;
         }
         s_offsets[dim_index][i] = offset;
      }
      s_offsets_are_set[dim_index] = true;
   }
}

}
}
#endif
