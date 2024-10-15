// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "MultiOrderPEnergyEvaluationStrategy.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"

using namespace SAMRAI;

MultiOrderPEnergyEvaluationStrategy::MultiOrderPEnergyEvaluationStrategy(
    const QuatModelParameters& model_parameters, const int phase_id,
    const int weight_id, const int energy_diag_id)
    : d_model_parameters(model_parameters),
      d_phase_id(phase_id),
      d_weight_id(weight_id),
      d_energy_diag_id(energy_diag_id)
{
}

void MultiOrderPEnergyEvaluationStrategy::evaluatePairEnergy(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   d_total_energy.resize(d_model_parameters.norderpA() *
                         d_model_parameters.norderpA());
   for (auto& e : d_total_energy)
      e = 0.;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         evaluatePairEnergy(patch, d_total_energy);
      }
   }

   sumReduction(&d_total_energy[0], d_total_energy.size());
}

void MultiOrderPEnergyEvaluationStrategy::printPairEnergy(std::ostream& os)
{
   const double threshold = 1.e-8;

   int n = d_model_parameters.norderpA();
   assert(d_total_energy.size() == n * n);

   os << "Interfacial pair energies:" << std::endl;
   for (int i = 0; i < n; i++) {
      os << i << ":";
      for (int j = 0; j < n; j++) {
         const double e = d_total_energy[i * n + j];
         if (std::abs(e) > threshold)
            os << " (" << j << ") " << d_total_energy[i * n + j];
      }
      os << std::endl;
   }
}

void MultiOrderPEnergyEvaluationStrategy::evaluatePairEnergy(
    std::shared_ptr<hier::Patch> patch, std::vector<double>& total_energy)
{
   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_id)));
   std::shared_ptr<pdat::CellData<double> > weight(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_weight_id)));

   assert(phase);
   assert(weight);

   int per_cell = 0;
   double* ptr_energy = nullptr;
   if (d_model_parameters.with_visit_energy_output()) {
      per_cell = 1;
      assert(d_energy_diag_id >= 0);
      std::shared_ptr<pdat::CellData<double> > energy(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_energy_diag_id)));
      ptr_energy = energy->getPointer();
   }

   assert(weight->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   // printf("orient_interp_func_type2 =%s\n",
   // d_model_parameters.orient_interp_func_type2().c_str()[0]);
   PHIPHI_INTERFACIALENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                            ifirst(2), ilast(2),
#endif
                            phase->getPointer(), phase->getGhostCellWidth()[0],
                            d_model_parameters.norderpA(),
                            d_model_parameters.gamma(),
                            d_model_parameters.m_moelans2011(),
                            weight->getPointer(), ptr_energy, &total_energy[0],
                            per_cell);
}
