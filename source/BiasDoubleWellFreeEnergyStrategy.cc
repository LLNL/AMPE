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
//
#include "BiasDoubleWellFreeEnergyStrategy.h"
#include "SAMRAI/pdat/CellData.h"

#include <cassert>

using namespace SAMRAI;


BiasDoubleWellFreeEnergyStrategy::BiasDoubleWellFreeEnergyStrategy() {}

//=======================================================================

void BiasDoubleWellFreeEnergyStrategy::addDrivingForceEta(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int eta_id, const int conc_id, const int f_l_id,
    const int f_a_id, const int f_b_id, const int rhs_id)
{
   (void)patch;
   (void)temperature_id;
   (void)phase_id;
   (void)eta_id;
   (void)conc_id;
   (void)f_l_id;
   (void)f_a_id;
   (void)f_b_id;
   (void)rhs_id;
   TBOX_ERROR(
       "BiasDoubleWellFreeEnergyStrategy::addDrivingForceEta() not "
       "implemented!!!\n");
}

//=======================================================================

void BiasDoubleWellFreeEnergyStrategy::computePhaseConcentrations(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int eta_id,
    const int concentration_id)
{
   (void)temperature_id;
   (void)phase_id;
   (void)eta_id;
   (void)concentration_id;
}

//=======================================================================

void BiasDoubleWellFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseL(
    const double temp, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_l;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

//=======================================================================

void BiasDoubleWellFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseA(
    const double temp, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_a;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}

//=======================================================================

void BiasDoubleWellFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseB(
    const double temp, const std::vector<double>& c_b,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;
   (void)c_b;
   (void)use_internal_units;

   d2fdc2.assign(d2fdc2.size(), 0.);
}
