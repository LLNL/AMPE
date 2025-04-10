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
#include "ParabolicTools.h"

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
   readParabolicData(input_db, "Liquid", coeffL);

   double coeffA[3][2];
   readParabolicData(input_db, "PhaseA", coeffA);

   double coeffB[3][2];
   readParabolicData(input_db, "PhaseB", coeffB);

   double Tref = input_db->getDouble("Tref");

   d_fenergy.reset(new Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase(
       Tref, coeffL, coeffA, coeffB, model_parameters.energy_interp_func_type(),
       Thermo4PFM::ConcInterpolationType::LINEAR));
}

ParabolicEquilibriumThreePhasesBinaryMultiOrder::
    ~ParabolicEquilibriumThreePhasesBinaryMultiOrder(){};
