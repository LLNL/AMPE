/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/

#ifndef included_hier_Index_C
#define included_hier_Index_C

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

namespace SAMRAI {
namespace hier {

Index * Index::s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
Index * Index::s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

Index * Index::s_mins[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
Index * Index::s_maxs[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

tbox::StartupShutdownManager::Handler
Index::s_initialize_finalize_handler(
   Index::initializeCallback,
   0,
   0,
   Index::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

Index::Index():
   IntVector()
{
}

Index::Index(
   const tbox::Dimension& dim):
   IntVector(dim)
{
   // an explicit setting Invalid is allowed.
   TBOX_DIM_ASSERT((!dim.isValid()) || (
         dim.getValue() > 0 && dim <= tbox::Dimension::getMaxDimension()));
}

Index::Index(
   const tbox::Dimension& dim,
   const int i):
   IntVector(dim, i)
{
   // an explicit setting Invalid is allowed.
   TBOX_DIM_ASSERT((!dim.isValid()) || (
         dim >= tbox::Dimension(1) && dim <= tbox::Dimension::getMaxDimension()));
}

Index::Index(
   const int i,
   const int j):
   IntVector(tbox::Dimension(2))
{
   TBOX_DIM_ASSERT(
      tbox::Dimension::getMaxDimension() >= tbox::Dimension(2));

   (*this)[0] = i;
   if (tbox::Dimension::MAXIMUM_DIMENSION_VALUE > 1) {
      (*this)[1] = j;
   }
}

Index::Index(
   const int i,
   const int j,
   const int k):
   IntVector(tbox::Dimension(3))
{
   TBOX_DIM_ASSERT(tbox::Dimension::getMaxDimension() >= tbox::Dimension(3));

   (*this)[0] = i;
   if (tbox::Dimension::MAXIMUM_DIMENSION_VALUE > 1) {
      (*this)[1] = j;
   }

   if (tbox::Dimension::MAXIMUM_DIMENSION_VALUE > 2) {
      (*this)[2] = k;
   }

}

Index::Index(
   const tbox::Array<int>& a):
   IntVector(a)
{
}

Index::Index(
   const tbox::Dimension& dim,
   const int array[]):
   IntVector(dim, array)
{
}

Index::Index(
   const Index& rhs):
   IntVector(rhs)
{
}

Index::Index(
   const IntVector& rhs):
   IntVector(rhs)
{
}

Index::~Index()
{
}

void
Index::initializeCallback()
{
   for (unsigned short d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      s_zeros[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)), 0);
      s_ones[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)), 1);

      s_mins[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)),
            tbox::MathUtilities<int>::getMin());
      s_maxs[d] = new Index(tbox::Dimension(static_cast<unsigned short>(d + 1)),
            tbox::MathUtilities<int>::getMax());
   }
}

void
Index::finalizeCallback()
{
   for (int d = 0; d < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++d) {
      delete s_zeros[d];
      delete s_ones[d];

      delete s_mins[d];
      delete s_maxs[d];
   }
}

}
}

#endif
