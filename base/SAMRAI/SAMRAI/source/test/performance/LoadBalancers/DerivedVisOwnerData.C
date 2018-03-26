/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   TreeLoadBalancer test.
 *
 ************************************************************************/
#include "DerivedVisOwnerData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

using namespace SAMRAI;

DerivedVisOwnerData::DerivedVisOwnerData()
{
}

DerivedVisOwnerData::~DerivedVisOwnerData()
{
}

bool DerivedVisOwnerData::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_id) const
{
   NULL_USE(patch);
   NULL_USE(depth_id);

   if (variable_name == "Owner") {
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      double owner = mpi.getRank();
      int i, size = region.size();
      for (i = 0; i < size; ++i) buffer[i] = owner;
   } else {
      // Did not register this name.
      TBOX_ERROR(
         "Unregistered variable name '" << variable_name << "' in\n"
                                        <<
         "DerivedVisOwnerData::packDerivedPatchDataIntoDoubleBuffer");
   }

   return true;
}
