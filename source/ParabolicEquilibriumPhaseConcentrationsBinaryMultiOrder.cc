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
#include "ParabolicEquilibriumPhaseConcentrationsBinaryMultiOrder.h"
#include "ParabolicFreeEnergyFunctionsBinary.h"
#include "FuncFort.h"
#include "ParabolicTools.h"

ParabolicEquilibriumPhaseConcentrationsBinaryMultiOrder::
    ParabolicEquilibriumPhaseConcentrationsBinaryMultiOrder(
        const int conc_l_id, const int conc_a_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrder(conc_l_id, conc_a_id,
                                                     conc_db)
{
   std::shared_ptr<tbox::Database> input_db = conc_db->getDatabase("Parabolic");

   double coeffL[3][2];
   readParabolicData(input_db, "Liquid", coeffL);

   double coeffA[3][2];
   readParabolicData(input_db, "PhaseA", coeffA);

   double Tref = input_db->getDouble("Tref");

   d_fenergy.reset(new Thermo4PFM::ParabolicFreeEnergyFunctionsBinary(
       Tref, coeffL, coeffA, energy_interp_func_type, conc_interp_func_type));
}
