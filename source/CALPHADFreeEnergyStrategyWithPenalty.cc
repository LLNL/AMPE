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
#include "CALPHADFreeEnergyStrategyWithPenalty.h"
#include "CALPHADConcSolverBinaryWithPenalty.h"
#include "CALPHADFreeEnergyFunctionsWithPenaltyBinary.h"

#include "SAMRAI/tbox/InputManager.h"

#include <string>


using namespace SAMRAI;

CALPHADFreeEnergyStrategyWithPenalty::CALPHADFreeEnergyStrategyWithPenalty(
    std::shared_ptr<tbox::Database> calphad_db,
    std::shared_ptr<tbox::Database> newton_db,
    const EnergyInterpolationType energy_interp_func_type,
    const ConcInterpolationType conc_interp_func_type,
    MolarVolumeStrategy* mvstrategy, const int conc_l_id, const int conc_a_id,
    const int conc_b_id, const int ncompositions, const bool with_third_phase)
    : CALPHADFreeEnergyStrategyBinary<CALPHADFreeEnergyFunctionsBinary>(
          calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type,
          mvstrategy, conc_l_id, conc_a_id, conc_b_id, with_third_phase)
{
   tbox::pout << "CALPHADFreeEnergyStrategyWithPenalty()..." << std::endl;
   const short n = with_third_phase ? 3 : 2;
   d_penalty_parameters.resize(n);
   for (short i = 0; i < n; i++)
      d_penalty_parameters[i].resize(6);

   std::string namemixL("PenaltyPhaseL");
   std::shared_ptr<tbox::Database> mixL_db = calphad_db->getDatabase(namemixL);
   mixL_db->getDoubleArray("Left", &d_penalty_parameters[0][0], 3);
   mixL_db->getDoubleArray("Right", &d_penalty_parameters[0][3], 3);

   std::string namemixA("PenaltyPhaseA");
   std::shared_ptr<tbox::Database> mixA_db = calphad_db->getDatabase(namemixA);
   mixA_db->getDoubleArray("Left", &d_penalty_parameters[1][0], 3);
   mixA_db->getDoubleArray("Right", &d_penalty_parameters[1][3], 3);

   if (with_third_phase) {
      std::string namemixB("PenaltyPhaseB");
      std::shared_ptr<tbox::Database> mixB_db =
          calphad_db->getDatabase(namemixB);
      mixB_db->getDoubleArray("Left", &d_penalty_parameters[2][0], 3);
      mixB_db->getDoubleArray("Right", &d_penalty_parameters[2][3], 3);
   }

   setup(calphad_db, newton_db);
}

//=======================================================================

void CALPHADFreeEnergyStrategyWithPenalty::setup(
    std::shared_ptr<tbox::Database> calphad_db,
    std::shared_ptr<tbox::Database> newton_db)
{
   tbox::pout << "CALPHADFreeEnergyStrategyWithPenalty::setupSolver()..."
              << std::endl;

   d_calphad_fenergy.reset(new CALPHADFreeEnergyFunctionsWithPenaltyBinary(
       calphad_db, newton_db, d_energy_interp_func_type,
       d_conc_interp_func_type, d_with_third_phase));
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyStrategyWithPenalty::computeCeqT(const double temperature,
                                                       const PhaseIndex pi0,
                                                       const PhaseIndex pi1,
                                                       double* ceq)
{
   assert(temperature > 0.);

   return d_calphad_fenergy->computeCeqT(temperature, pi0, pi1, ceq);
}

double CALPHADFreeEnergyStrategyWithPenalty::computeValFreeEnergyLiquid(
    const double temperature, const double conc, const bool gp)
{
   const double f1 =
       d_calphad_fenergy->computeFreeEnergy(temperature, &conc,
                                            PhaseIndex::phaseL, gp);
   const double f2 =
       d_calphad_fenergy->computePenalty(PhaseIndex::phaseL, conc);

   return (f1 + f2) * d_mv_strategy->computeInvMolarVolume(temperature, &conc,
                                                           PhaseIndex::phaseL);
}

double CALPHADFreeEnergyStrategyWithPenalty::computeValFreeEnergySolidA(
    const double temperature, const double conc, const bool gp)
{
   const double f1 =
       d_calphad_fenergy->computeFreeEnergy(temperature, &conc,
                                            PhaseIndex::phaseA, gp);
   const double f2 =
       d_calphad_fenergy->computePenalty(PhaseIndex::phaseA, conc);

   return (f1 + f2) * d_mv_strategy->computeInvMolarVolume(temperature, &conc,
                                                           PhaseIndex::phaseA);
}

double CALPHADFreeEnergyStrategyWithPenalty::computeValFreeEnergySolidB(
    const double temperature, const double conc, const bool gp)
{
   const double f1 =
       d_calphad_fenergy->computeFreeEnergy(temperature, &conc,
                                            PhaseIndex::phaseB, gp);
   const double f2 =
       d_calphad_fenergy->computePenalty(PhaseIndex::phaseB, conc);

   return (f1 + f2) * d_mv_strategy->computeInvMolarVolume(temperature, &conc,
                                                           PhaseIndex::phaseB);
}
