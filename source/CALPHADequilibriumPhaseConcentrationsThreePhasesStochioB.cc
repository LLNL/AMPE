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
#include "Database2JSON.h"
namespace pt = boost::property_tree;

#include "CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

// Thermo4PFM
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB ::
    CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB(
        const double concStochioB, const int conc_l_scratch_id,
        const int conc_a_scratch_id, const int conc_b_scratch_id,
        const int conc_l_ref_id, const int conc_a_ref_id,
        const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : CALPHADequilibriumPhaseConcentrationsStrategy<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>(
          conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
          conc_l_ref_id, conc_a_ref_id, conc_b_ref_id, energy_interp_func_type,
          Thermo4PFM::ConcInterpolationType::LINEAR, false, calphad_pt,
          newton_db, ncompositions),
      d_concStochioB(concStochioB)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy = std::unique_ptr<
       Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>(
       new Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB(
           concStochioB, calphad_pt, newton_pt, energy_interp_func_type,
           Thermo4PFM::ConcInterpolationType::LINEAR));
}
