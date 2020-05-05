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
#include "MinIntCoarsen.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>

using namespace SAMRAI;


MinIntCoarsen::MinIntCoarsen() : hier::CoarsenOperator("BASE_MIN_COARSEN")
{
   d_name_id = "MINIMUM_COARSEN";
}

MinIntCoarsen::~MinIntCoarsen() {}

bool MinIntCoarsen::findCoarsenOperator(
    const std::shared_ptr<hier::Variable>& var,
    const std::string& op_name) const
{
   const std::shared_ptr<pdat::CellVariable<int> > cast_var(
       SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<int>, hier::Variable>(var));
   if (cast_var && (op_name == d_name_id)) {
      return (true);
   } else {
      return (false);
   }
}

const std::string& MinIntCoarsen::getOperatorName() const
{
   return (d_name_id);
}

int MinIntCoarsen::getOperatorPriority() const { return (0); }

hier::IntVector MinIntCoarsen::getStencilWidth(const tbox::Dimension& dim) const
{
   return (hier::IntVector(dim, 0));
}


// Take the minimum positive integer of cell-centered integer patch
// data.  Negative values are ignored, and if all fine values are
// negative, the coarse value will be set to -1.

void MinIntCoarsen::coarsen(hier::Patch& coarse, const hier::Patch& fine,
                            const int dst_component, const int src_component,
                            const hier::Box& coarse_box,
                            const hier::IntVector& ratio) const
{
   std::shared_ptr<pdat::CellData<int> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<int> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
           coarse.getPatchData(dst_component)));
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(fdata);
   assert(cdata);
   assert(cdata->getDepth() == fdata->getDepth());
#endif

   pdat::CellIterator cend(pdat::CellGeometry::end(coarse_box));
   for (pdat::CellIterator ci(pdat::CellGeometry::begin(coarse_box));
        ci != cend; ++ci) {
      pdat::CellIndex coarse_cell = *ci;
      hier::Index fine_lower = coarse_cell * ratio;
      hier::Index fine_upper = fine_lower + ratio - 1;

      hier::Box fine_box(fine_lower, fine_upper, hier::BlockId::zero());

      (*cdata)(coarse_cell) = -1;

      pdat::CellIterator fend(pdat::CellGeometry::end(fine_box));
      for (pdat::CellIterator fi(pdat::CellGeometry::begin(fine_box));
           fi != fend; ++fi) {
         pdat::CellIndex fine_cell = *fi;

         int fine_val = (*fdata)(fine_cell);
         int coarse_val = (*cdata)(coarse_cell);

         if (fine_val >= 0) {
            if (coarse_val < 0 || fine_val < coarse_val) {
               (*cdata)(coarse_cell) = fine_val;
            }
         }
      }
   }
}
