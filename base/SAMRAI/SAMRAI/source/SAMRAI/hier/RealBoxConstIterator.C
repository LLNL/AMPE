/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Iterator over real Boxes in a BoxContainer.
 *
 ************************************************************************/
#ifndef included_hier_RealBoxConstIterator_C
#define included_hier_RealBoxConstIterator_C

#include "SAMRAI/hier/RealBoxConstIterator.h"

namespace SAMRAI {
namespace hier {

RealBoxConstIterator::RealBoxConstIterator(
   const BoxContainer& mapped_boxes,
   bool begin):
   d_mapped_boxes(&mapped_boxes),
   d_ni(begin ? d_mapped_boxes->begin() : d_mapped_boxes->end())
{
   if (begin)
   {
      while (d_ni != d_mapped_boxes->end() && d_ni->isPeriodicImage()) {
         ++d_ni;
      }
   }
}

RealBoxConstIterator::~RealBoxConstIterator()
{
   d_mapped_boxes = NULL;
}

/*
 ****************************************************************************
 * Pre-increment operator.
 ****************************************************************************
 */

RealBoxConstIterator&
RealBoxConstIterator::operator ++ ()
{
   do {
      ++d_ni;
   } while (d_ni != d_mapped_boxes->end() && d_ni->isPeriodicImage());
   return *this;
}

/*
 ****************************************************************************
 * Post-increment operator.
 ****************************************************************************
 */

RealBoxConstIterator
RealBoxConstIterator::operator ++ (
   int)
{
   RealBoxConstIterator saved = *this;
   do {
      ++d_ni;
   } while (d_ni != d_mapped_boxes->end() && d_ni->isPeriodicImage());
   return saved;
}

}
}
#endif
