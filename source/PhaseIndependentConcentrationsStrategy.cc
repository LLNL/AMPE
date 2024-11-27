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
#include "PhaseIndependentConcentrationsStrategy.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"

int PhaseIndependentConcentrationsStrategy::computePhaseConcentrationsOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_concentration,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_b,
    std::shared_ptr<hier::Patch> patch)
{
   (void)cd_temperature;
   (void)cd_phi;
   (void)cd_c_b;

   assert(cd_concentration);
   assert(cd_c_l);
   assert(cd_c_a);

   const hier::Box& pbox = patch->getBox();

   math::PatchCellDataOpsReal<double> cellops;

   cellops.copyData(cd_c_l, cd_concentration, pbox);
   cellops.copyData(cd_c_a, cd_concentration, pbox);

   return 0;
}
