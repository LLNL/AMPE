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
#include "ParabolicEquilibriumThreePhasesBinaryMultiOrder.h"
#include "FuncFort.h"


ParabolicEquilibriumThreePhasesBinaryMultiOrder::
    ParabolicEquilibriumThreePhasesBinaryMultiOrder(
        const short norderp_A, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const QuatModelParameters& model_parameters,
        std::shared_ptr<tbox::Database> conc_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases(
          norderp_A, conc_l_id, conc_a_id, conc_b_id, model_parameters, conc_db)
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
   coeffL[0][0] = phasea_db->getDouble("a0");
   coeffL[0][1] = phasea_db->getDouble("a1");
   coeffL[1][0] = phasea_db->getDouble("b0");
   coeffL[1][1] = phasea_db->getDouble("b1");
   coeffL[2][0] = phasea_db->getDouble("c0");
   coeffL[2][1] = phasea_db->getDouble("c1");

   double coeffB[3][2];
   std::shared_ptr<tbox::Database> phaseb_db = input_db->getDatabase("PhaseB");
   coeffB[0][0] = phaseb_db->getDouble("a0");
   coeffB[0][1] = phaseb_db->getDouble("a1");
   coeffB[1][0] = phaseb_db->getDouble("b0");
   coeffB[1][1] = phaseb_db->getDouble("b1");
   coeffB[2][0] = phaseb_db->getDouble("c0");
   coeffB[2][1] = phaseb_db->getDouble("c1");

   double Tref = input_db->getDouble("Tref");

   d_fenergy.reset(new Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase(
       Tref, coeffL, coeffA, coeffB, model_parameters.energy_interp_func_type(),
       Thermo4PFM::ConcInterpolationType::LINEAR));
}

ParabolicEquilibriumThreePhasesBinaryMultiOrder::
    ~ParabolicEquilibriumThreePhasesBinaryMultiOrder(){};
