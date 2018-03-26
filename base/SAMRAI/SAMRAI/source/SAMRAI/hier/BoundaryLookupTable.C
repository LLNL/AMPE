/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Lookup table to aid in BoundaryBox construction
 *
 ************************************************************************/

#ifndef included_hier_BoundaryLookupTable_C
#define included_hier_BoundaryLookupTable_C

#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

BoundaryLookupTable *
BoundaryLookupTable::s_lookup_table_instance[tbox::Dimension::MAXIMUM_DIMENSION_VALUE
] = { (BoundaryLookupTable *)NULL };

tbox::StartupShutdownManager::Handler
BoundaryLookupTable::s_finalize_handler(
   0,
   0,
   0,
   BoundaryLookupTable::finalizeCallback,
   tbox::StartupShutdownManager::priorityBoundaryLookupTable);

/*
 *************************************************************************
 *
 * Lookup table constructor and destructor.
 *
 *************************************************************************
 */

BoundaryLookupTable::BoundaryLookupTable(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   if (d_table[0].isNull()) {
      int factrl[tbox::Dimension::MAXIMUM_DIMENSION_VALUE + 1];
      factrl[0] = 1;
      for (int i = 1; i <= d_dim.getValue(); i++) factrl[i] = i * factrl[i - 1];
      d_ncomb.resizeArray(d_dim.getValue());
      d_max_li.resizeArray(d_dim.getValue());
      for (int codim = 1; codim <= d_dim.getValue(); codim++) {
         int cdm1 = codim - 1;
         d_ncomb[cdm1] = factrl[d_dim.getValue()]
            / (factrl[codim] * factrl[d_dim.getValue() - codim]);

         tbox::Array<int> work;
         work.resizeArray(codim * d_ncomb[cdm1]);

         int recursive_work[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
         int recursive_work_lvl = 0;
         int* recursive_work_ptr;
         buildTable(work.getPointer(), codim, 1, recursive_work,
            recursive_work_lvl, recursive_work_ptr);

         d_table[cdm1].resizeArray(d_ncomb[cdm1]);
         for (int j = 0; j < d_ncomb[cdm1]; j++) {
            d_table[cdm1][j].resizeArray(codim);
            for (int k = 0; k < codim; k++) {
               d_table[cdm1][j][k] = work[j * codim + k] - 1;
            }
         }

         d_max_li[cdm1] = d_ncomb[cdm1] * (1 << codim);
      }
   }

   buildBoundaryDirectionVectors();
}

BoundaryLookupTable::~BoundaryLookupTable()
{
}

/*
 *************************************************************************
 *
 * Recursive function that computes the combinations in the lookup
 * table.
 *
 *************************************************************************
 */

void
BoundaryLookupTable::buildTable(
   int* table,
   int codim,
   int ibeg,
   int(&work)[tbox::Dimension::MAXIMUM_DIMENSION_VALUE],
   int& lvl,
   int *& ptr)
{
   lvl++;
   if (lvl == 1) ptr = table;
   int iend = d_dim.getValue() - codim + lvl;
   for (int i = ibeg; i <= iend; i++) {
      work[lvl - 1] = i;
      if (lvl != codim) {
         buildTable(ptr, codim, i + 1, work, lvl, ptr);
      } else {
         for (int j = 0; j < codim; j++) {
            *(ptr + j) = work[j];
         }
         ptr += codim;
      }
   }
   lvl--;
}

/*
 *************************************************************************
 *
 * Build table of IntVectors indication locations of boundaries relative
 * to a patch.
 *
 *************************************************************************
 */

void
BoundaryLookupTable::buildBoundaryDirectionVectors()
{

   d_bdry_dirs.resizeArray(d_dim.getValue());

   for (int i = 0; i < d_dim.getValue(); i++) {
      d_bdry_dirs[i].resizeArray(d_max_li[i], IntVector::getZero(d_dim));
      int codim = i + 1;

      for (int loc = 0; loc < d_max_li[i]; loc++) {
         const tbox::Array<int>& dirs = getDirections(loc, codim);

         for (int d = 0; d < dirs.size(); d++) {

            if (isUpper(loc, codim, d)) {

               d_bdry_dirs[i][loc](dirs[d]) = 1;

            } else {

               d_bdry_dirs[i][loc](dirs[d]) = -1;

            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Free statics.
 *
 *************************************************************************
 */
void
BoundaryLookupTable::finalizeCallback()
{
   for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++i) {
      if (s_lookup_table_instance[i]) {
         delete s_lookup_table_instance[i];
      }
      s_lookup_table_instance[i] = ((BoundaryLookupTable *)NULL);
   }

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
