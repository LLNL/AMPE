// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "PhaseConcentrationsStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

PhaseConcentrationsStrategy::PhaseConcentrationsStrategy(const int conc_l_id,
                                                         const int conc_a_id,
                                                         const int conc_b_id)
    : d_conc_l_id(conc_l_id), d_conc_a_id(conc_a_id), d_conc_b_id(conc_b_id)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
}

void PhaseConcentrationsStrategy::computePhaseConcentrations(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int concentration_id)
{
   // tbox::pout<<"CALPHADFreeEnergyStrategy::computePhaseConcentrations()"<<endl;

   assert(temperature_id >= 0);
   assert(phase_id >= 0);
   assert(concentration_id >= 0);
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
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
   int nits = 0;

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
         if (d_conc_b_id > -1) {
            c_b =
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_conc_b_id));
            assert(c_b);
         }

         int ret =
             computePhaseConcentrationsOnPatch(temperature, phi, concentration,
                                               c_l, c_a, c_b, patch);
         if (ret < 0) tbox::SAMRAI_MPI::abort();
         nits += ret;
      }
   }
#if 0
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   mpi.AllReduce(&nits, 1, MPI_SUM);

   tbox::pout << "computePhaseConcentrations: total nits = " << nits
              << std::endl;
#endif
}
