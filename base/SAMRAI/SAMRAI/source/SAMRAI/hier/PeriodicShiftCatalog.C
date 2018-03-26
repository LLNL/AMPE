/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Set of edges incident from a mapped_box_level of a distributed box graph.
 *
 ************************************************************************/
#ifndef included_hier_PeriodicShiftCatalog_C
#define included_hier_PeriodicShiftCatalog_C

#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include <iomanip>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

PeriodicShiftCatalog * PeriodicShiftCatalog::s_periodic_shift_catalog_instance[
   tbox::Dimension::MAXIMUM_DIMENSION_VALUE] = { 0 };

tbox::StartupShutdownManager::Handler
PeriodicShiftCatalog::s_finalize_handler(
   0,
   0,
   0,
   PeriodicShiftCatalog::finalizeCallback,
   40);

/*
 ***********************************************************************
 ***********************************************************************
 */

const PeriodicShiftCatalog*
PeriodicShiftCatalog::getCatalog(
   const tbox::Dimension& dim)
{
   if (s_periodic_shift_catalog_instance[dim.getValue() - 1] == NULL) {
      s_periodic_shift_catalog_instance[dim.getValue() - 1] = new PeriodicShiftCatalog(dim);
   }
   return s_periodic_shift_catalog_instance[dim.getValue() - 1];
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
PeriodicShiftCatalog::finalizeCallback()
{
   for (int i = 0; i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; ++i) {
      if (s_periodic_shift_catalog_instance[i] != NULL) {
         delete s_periodic_shift_catalog_instance[i];
         s_periodic_shift_catalog_instance[i] = NULL;
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

PeriodicShiftCatalog::PeriodicShiftCatalog(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_shifts(1, IntVector::getZero(dim)),
   d_opposite_number(1),
   d_zero_shift_number(0)
{
   d_opposite_number[0] = 0;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

PeriodicShiftCatalog::~PeriodicShiftCatalog()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
PeriodicShiftCatalog::setShifts(
   const tbox::Dimension& dim,
   const std::vector<IntVector>& shifts)
{
   getCatalog(dim);  // Causes singleton object creation if it does not yet exist.

   const int dim_index(dim.getValue() - 1);

   std::vector<IntVector> tmp_shifts;
   s_periodic_shift_catalog_instance[dim_index]->d_opposite_number.clear();

   const IntVector& zero_shift(IntVector::getZero(dim));

   // The first position is the zero-shift and its own opposite.
   s_periodic_shift_catalog_instance[dim_index]->d_opposite_number.push_back(
      static_cast<PeriodicId>(static_cast<int>(tmp_shifts.size())));
   tmp_shifts.push_back(zero_shift);

   /*
    * Add the shifts in shifts, avoiding duplicate shifts.
    */
   for (std::vector<IntVector>::const_iterator vi = shifts.begin();
        vi != shifts.end(); ++vi) {

      std::vector<IntVector>::const_iterator vj;
      for (vj = tmp_shifts.begin(); vj != tmp_shifts.end(); ++vj) {
         if (*vi == *vj) {
            break;
         }
      }

      if (vj == tmp_shifts.end()) {
         /*
          * Shift *vi doesn't already exists.  Add it, add it's
          * reverse, and note the positions in d_opposite_number.
          */
         PeriodicId tmpId1(static_cast<int>(tmp_shifts.size()));
         PeriodicId tmpId2(static_cast<int>(tmp_shifts.size() + 1));
         s_periodic_shift_catalog_instance[dim_index]->d_opposite_number.push_back(tmpId2);
         s_periodic_shift_catalog_instance[dim_index]->d_opposite_number.push_back(tmpId1);
         tmp_shifts.push_back(*vi);
         tmp_shifts.push_back(-(*vi));
      }

   }

   s_periodic_shift_catalog_instance[dim_index]->d_shifts = tmp_shifts;
   s_periodic_shift_catalog_instance[dim_index]->d_zero_shift_number = 0;

   if (1) {
      // Write out the shift catalog to log file.
      tbox::plog << "\n\nPeriodicShiftCatalog has "
                 << s_periodic_shift_catalog_instance[dim_index]->d_shifts.size()
                 << " shifts:\n";
      tbox::plog << "Shift   Opposite\n";
      tbox::plog << "Number  Shift     Shift\n";
      for (size_t i = 0;
           i < s_periodic_shift_catalog_instance[dim_index]->d_shifts.size();
           ++i) {
         tbox::plog << std::setw(3) << i << "      "
                    << std::setw(3)
                    << s_periodic_shift_catalog_instance[dim_index]->d_opposite_number[i]
                    << "      "
                    << s_periodic_shift_catalog_instance[dim_index]->d_shifts[i] << "\n";
      }
      tbox::plog << "\n\n";
   }

   TBOX_ASSERT(s_periodic_shift_catalog_instance[dim_index]->d_shifts.size() ==
      s_periodic_shift_catalog_instance[dim_index]->d_opposite_number.size());
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
PeriodicShiftCatalog::initializeShiftsByIndexDirections(
   const IntVector& shift_distance_along_index_directions)
{

   const tbox::Dimension& dim(shift_distance_along_index_directions.getDim());

   // Compute the number of shifts, 3^DIM.
   int num_shift = 3;
   for (int d = 1; d < dim.getValue(); ++d) {
      num_shift *= 3;
   }

   const IntVector& zero_shift(IntVector::getZero(dim));

   std::vector<IntVector> shifts;
   shifts.clear();
   shifts.insert(shifts.begin(), num_shift, zero_shift);

   /*
    * Compute the direction of each shift, a value of -1, 0 or 1
    * independent of the period.  A positive value means shift in
    * the positive direction, etc.
    *
    * After getting the signs, multiply in the distances.
    */
   for (int s = 0; s < num_shift; ++s) {
      int s1 = s;
      for (int d = 0; d < dim.getValue(); ++d) {
         shifts[s](d) = s1 % 3 - 1;
         s1 /= 3;
      }
   }

   for (unsigned int s = 0; s < shifts.size(); ++s) {
      shifts[s] *= shift_distance_along_index_directions;
   }

   setShifts(dim, shifts);
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
