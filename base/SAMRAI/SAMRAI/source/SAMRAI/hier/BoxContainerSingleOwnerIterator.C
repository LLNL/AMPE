/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Special iterator for BoxContainer.
 *
 ************************************************************************/
#ifndef included_hier_BoxContainerSingleOwnerIterator_C
#define included_hier_BoxContainerSingleOwnerIterator_C

#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"

namespace SAMRAI {
namespace hier {

BoxContainerSingleOwnerIterator::BoxContainerSingleOwnerIterator(
   const BoxContainer& mapped_boxes,
   const int& owner_rank,
   bool begin):
   d_mapped_boxes(&mapped_boxes),
   d_owner_rank(owner_rank),
   d_iter(begin ? d_mapped_boxes->begin() : d_mapped_boxes->end())
{
   if (begin)
   {
      while (d_iter != d_mapped_boxes->end() &&
             d_iter->getOwnerRank() != d_owner_rank) {
         ++d_iter;
      }
   }
}

BoxContainerSingleOwnerIterator::~BoxContainerSingleOwnerIterator()
{
   d_mapped_boxes = NULL;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

BoxContainerSingleOwnerIterator&
BoxContainerSingleOwnerIterator::operator ++ ()
{
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

BoxContainerSingleOwnerIterator
BoxContainerSingleOwnerIterator::operator ++ (
   int)
{
   BoxContainerSingleOwnerIterator saved = *this;
   do {
      ++d_iter;
   } while (d_iter != d_mapped_boxes->end() &&
            d_iter->getOwnerRank() != d_owner_rank);
   return saved;
}

}
}
#endif
