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
#include "PhaseConcentrationsStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

PhaseConcentrationsStrategy::PhaseConcentrationsStrategy(
    const int conc_l_id, const int conc_a_id, const int conc_b_id,
    const bool with_third_phase)
    : d_conc_l_id(conc_l_id),
      d_conc_a_id(conc_a_id),
      d_conc_b_id(conc_b_id),
      d_with_third_phase(with_third_phase)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
}

void PhaseConcentrationsStrategy::computePhaseConcentrations(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int eta_id,
    const int concentration_id)
{
   // tbox::pout<<"CALPHADFreeEnergyStrategy::computePhaseConcentrations()"<<endl;

   assert(temperature_id >= 0);
   assert(phase_id >= 0);
   assert(concentration_id >= 0);
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   if (d_with_third_phase) {
      assert(eta_id >= 0);
      assert(d_conc_b_id >= 0);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
   assert(cellops.max(phase_id) == cellops.max(phase_id));
   double maxphi = cellops.max(phase_id);
   double minphi = cellops.min(phase_id);
   assert(maxphi >= 0.);
   assert(maxphi < 1.1);
   assert(minphi >= -0.1);
   assert(minphi <= 1.);
   double maxc = cellops.max(concentration_id);
   assert(maxc == maxc);
#endif

   const int maxl = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double> > temperature(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));
         assert(temperature);

         std::shared_ptr<pdat::CellData<double> > phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));
         assert(phi);
#ifdef DEBUG_CHECK_ASSERTIONS
         const hier::Box& pbox = patch->getBox();

         SAMRAI::math::PatchCellDataNormOpsReal<double> ops;
         double l2phi = ops.L2Norm(phi, pbox);
         assert(l2phi == l2phi);
         assert(l2phi >= 0.);
         assert(l2phi < 1000.);
#endif

         std::shared_ptr<pdat::CellData<double> > eta;
         if (d_with_third_phase) {
            eta = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(eta_id));
         }

         std::shared_ptr<pdat::CellData<double> > concentration(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(concentration_id)));
         assert(concentration);

         std::shared_ptr<pdat::CellData<double> > c_l(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_l_id)));
         assert(c_l);

         std::shared_ptr<pdat::CellData<double> > c_a(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_a_id)));
         assert(c_a);

         std::shared_ptr<pdat::CellData<double> > c_b;
         if (d_with_third_phase) {
            c_b = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(d_conc_b_id));
            assert(c_b);
         }

         computePhaseConcentrationsOnPatch(temperature, phi, eta, concentration,
                                           c_l, c_a, c_b, patch);
      }
   }
}
