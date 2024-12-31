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
#include "KKSdiluteEquilibriumPhaseConcentrationsStrategy.h"

#include "Database2JSON.h"

namespace pt = boost::property_tree;

KKSdiluteEquilibriumPhaseConcentrationsStrategy::
    KKSdiluteEquilibriumPhaseConcentrationsStrategy(
        const int conc_l_id, const int conc_a_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db)
    : EquilibriumPhaseConcentrationsBinary(conc_l_id, conc_a_id,
                                           conc_interp_func_type)
{
   pt::ptree troot;
   copyDatabase(conc_db, troot);
   d_fenergy.reset(new Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary(
       troot, energy_interp_func_type, conc_interp_func_type));
}
