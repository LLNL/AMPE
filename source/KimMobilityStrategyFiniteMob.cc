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
#include "KimMobilityStrategyFiniteMob.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"

template <class FreeEnergyType>
KimMobilityStrategyFiniteMob<FreeEnergyType>::KimMobilityStrategyFiniteMob(
    QuatModel* quat_model, const int conc_l_id, const int conc_s_id,
    const int temp_id, const double interface_mobility, const double epsilon,
    const double phase_well_scale,
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions)
    : KimMobilityStrategy<FreeEnergyType>(quat_model, conc_l_id, conc_s_id, -1,
                                          temp_id, energy_interp_func_type,
                                          conc_interp_func_type, conc_db,
                                          ncompositions)
{
   assert(epsilon > 0.);
   assert(phase_well_scale >= 0.);
   assert(interface_mobility > 0.);

   const double xi = epsilon / sqrt(16. * phase_well_scale);

   d_alpha = 3. * sqrt(2.) * xi / interface_mobility;
}

template class KimMobilityStrategyFiniteMob<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>;
template class KimMobilityStrategyFiniteMob<
    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>;
template class KimMobilityStrategyFiniteMob<
    Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary>;
template class KimMobilityStrategyFiniteMob<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>;
