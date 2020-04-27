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
#include "LinearMeltingTemperatureStrategy.h"
#include "ConcFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"

LinearMeltingTemperatureStrategy::LinearMeltingTemperatureStrategy(
    const double Tref, const double c0, const double liquidus_slope,
    const int concentration_id, const int equilibrium_temperature_id)
    : d_Tref(Tref),
      d_c0(c0),
      d_liquidus_slope(liquidus_slope),
      d_concentration_id(concentration_id),
      d_equilibrium_temperature_id(equilibrium_temperature_id)
{
   tbox::plog << "uses LinearMeltingTemperatureStrategy with Tref = " << Tref
              << "..." << std::endl;

   assert(d_concentration_id >= 0);
   assert(d_equilibrium_temperature_id >= 0);
   assert(d_Tref >= 0.);
   assert(d_Tref < 100000.);
}

void LinearMeltingTemperatureStrategy::evaluate(hier::Patch& patch)
{
   assert(d_equilibrium_temperature_id >= 0);
   assert(d_concentration_id >= 0);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_equilibrium_temperature_id)));
   std::shared_ptr<pdat::CellData<double> > concentration(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_concentration_id)));

   assert(temperature);
   assert(concentration);

#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> ops;
   double l2rhs = ops.L2Norm(concentration, patch.getBox());
   assert(l2rhs == l2rhs);
#endif

   FORT_LINEARMELTINGLINE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                          ifirst(2), ilast(2),
#endif
                          temperature->getPointer(),
                          temperature->getGhostCellWidth()[0],
                          concentration->getPointer(),
                          concentration->getGhostCellWidth()[0], d_Tref, d_c0,
                          d_liquidus_slope);

#ifdef DEBUG_CHECK_ASSERTIONS
   l2rhs = ops.L2Norm(temperature, patch.getBox());
   assert(l2rhs == l2rhs);
#endif
}
