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
#include "BiasDoubleWellBeckermannFreeEnergyStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "QuatFort.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cassert>

using namespace SAMRAI;


BiasDoubleWellBeckermannFreeEnergyStrategy::
    BiasDoubleWellBeckermannFreeEnergyStrategy(
        const double alpha, MeltingTemperatureStrategy* meltingTstrat)
    : d_alpha(alpha), d_meltingTstrat(meltingTstrat)
{
   tbox::plog << "BiasDoubleWellBeckermannFreeFreeEnergyStrategy:" << std::endl;
   tbox::plog << "alpha =" << d_alpha << std::endl;
}

//=======================================================================

void BiasDoubleWellBeckermannFreeEnergyStrategy::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int eta_id, const int conc_id, const int f_l_id,
    const int f_a_id, const int f_b_id, const int rhs_id)
{
   assert(phase_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);
   assert(d_meltingTstrat != NULL);

   (void)conc_id;  // unused
   (void)f_l_id;   // unused
   (void)f_a_id;   // unused

   boost::shared_ptr<pdat::CellData<double> > phase(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   boost::shared_ptr<pdat::CellData<double> > temp(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(temp);

   boost::shared_ptr<pdat::CellData<double> > rhs(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));
   assert(rhs);

   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   // evaluate melting temperature field
   //(which may depend on composition in linear model for instance)
   d_meltingTstrat->evaluate(patch);

   boost::shared_ptr<pdat::CellData<double> > eq_temp(
       BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_meltingTstrat->equilibrium_temperature_id())));
   assert(eq_temp);
#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> ops;
   double l2rhs = ops.L2Norm(eq_temp, patch.getBox());
   assert(l2rhs == l2rhs);
#endif


   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   COMPUTERHSBIASWELLBECKERMANN(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                ifirst(2), ilast(2),
#endif
                                phase->getPointer(),
                                phase->getGhostCellWidth()[0],
                                temp->getPointer(),
                                temp->getGhostCellWidth()[0], d_alpha,
                                eq_temp->getPointer(), 0, rhs->getPointer(), 0);
}
