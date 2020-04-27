// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"


using namespace SAMRAI;


void copyDepthSideData(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                       const int dst_id, const int dst_depth, const int src_id,
                       const int src_depth)
{
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::SideData<double> > dst(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 p->getPatchData(dst_id)));
         TBOX_ASSERT(dst);

         std::shared_ptr<pdat::SideData<double> > src(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 p->getPatchData(src_id)));
         TBOX_ASSERT(src);

         dst->copyDepth(dst_depth, *src, src_depth);
      }
   }
}

void copyDepthCellData(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                       const int dst_id, const int dst_depth, const int src_id,
                       const int src_depth)
{
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::CellData<double> > dst(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 p->getPatchData(dst_id)));
         TBOX_ASSERT(dst);

         std::shared_ptr<pdat::CellData<double> > src(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 p->getPatchData(src_id)));
         TBOX_ASSERT(src);

         dst->copyDepth(dst_depth, *src, src_depth);
      }
   }
}

int checkForNans(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                 const int data_id)
{
#ifndef NDEBUG
   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
   const double norm_data = mathops.L1Norm(data_id);

   if (norm_data != norm_data) return 1;
#else
   (void)hierarchy;
   (void)data_id;
#endif
   return 0;
}

int checkSideDataForNans(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, const int data_id)
{
#ifndef NDEBUG
   math::HierarchySideDataOpsReal<double> mathops(hierarchy);
   const double norm_data = mathops.L1Norm(data_id);

   if (norm_data != norm_data) return 1;
#else
   (void)hierarchy;
   (void)data_id;
#endif
   return 0;
}
