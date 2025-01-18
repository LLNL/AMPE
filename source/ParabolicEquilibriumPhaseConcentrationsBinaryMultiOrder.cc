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
   std::shared_ptr<tbox::Database> liquid_db = input_db->getDatabase("Liquid");
   coeffL[0][0] = liquid_db->getDouble("a0");
   coeffL[0][1] = liquid_db->getDouble("a1");
   coeffL[1][0] = liquid_db->getDouble("b0");
   coeffL[1][1] = liquid_db->getDouble("b1");
   coeffL[2][0] = liquid_db->getDouble("c0");
   coeffL[2][1] = liquid_db->getDouble("c1");

   double coeffA[3][2];
   std::shared_ptr<tbox::Database> phasea_db = input_db->getDatabase("PhaseA");
   coeffA[0][0] = phasea_db->getDouble("a0");
   coeffA[0][1] = phasea_db->getDouble("a1");
   coeffA[1][0] = phasea_db->getDouble("b0");
   coeffA[1][1] = phasea_db->getDouble("b1");
   coeffA[2][0] = phasea_db->getDouble("c0");
   coeffA[2][1] = phasea_db->getDouble("c1");

   double Tref = input_db->getDouble("Tref");

   d_fenergy.reset(new Thermo4PFM::ParabolicFreeEnergyFunctionsBinary(
       Tref, coeffL, coeffA, energy_interp_func_type, conc_interp_func_type));
}
