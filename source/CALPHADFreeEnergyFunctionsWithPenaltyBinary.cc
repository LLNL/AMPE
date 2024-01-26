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
#include "CALPHADFreeEnergyFunctionsWithPenaltyBinary.h"
#include "CALPHADConcSolverBinaryWithPenalty.h"
#include "CALPHADEqConcSolverBinaryWithPenalty.h"
#include "PhysicalConstants.h"


CALPHADFreeEnergyFunctionsWithPenaltyBinary::
    CALPHADFreeEnergyFunctionsWithPenaltyBinary(
        std::shared_ptr<SAMRAI::tbox::Database> calphad_db,
        std::shared_ptr<SAMRAI::tbox::Database> newton_db,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        const bool with_third_phase)
    : CALPHADFreeEnergyFunctionsBinary(calphad_db, newton_db,
                                       energy_interp_func_type,
                                       conc_interp_func_type, with_third_phase)
{
   const short n = with_third_phase ? 3 : 2;
   d_penalty_parameters.resize(n);
   for (short i = 0; i < n; i++)
      d_penalty_parameters[i].resize(6);

   readParameters(calphad_db);

   setupSolver(newton_db);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsWithPenaltyBinary::setupSolver(
    std::shared_ptr<tbox::Database> newton_db)
{
   tbox::pout << "CALPHADFreeEnergyFunctionsWithPenaltyBinary::setupSolver()..."
              << std::endl;

   if (d_solver != NULL) delete d_solver;
   d_solver =
       new CALPHADConcentrationSolverBinaryWithPenalty(d_with_third_phase,
                                                       d_penalty_parameters);

   readNewtonparameters(newton_db);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsWithPenaltyBinary::readParameters(
    std::shared_ptr<SAMRAI::tbox::Database> calphad_db)
{
   std::string namemixL("PenaltyPhaseL");
   std::shared_ptr<tbox::Database> mixL_db = calphad_db->getDatabase(namemixL);
   mixL_db->getDoubleArray("Left", &d_penalty_parameters[0][0], 3);
   mixL_db->getDoubleArray("Right", &d_penalty_parameters[0][3], 3);

   std::string namemixA("PenaltyPhaseA");
   std::shared_ptr<tbox::Database> mixA_db = calphad_db->getDatabase(namemixA);
   mixA_db->getDoubleArray("Left", &d_penalty_parameters[1][0], 3);
   mixA_db->getDoubleArray("Right", &d_penalty_parameters[1][3], 3);

   if (d_with_third_phase) {
      std::string namemixB("PenaltyPhaseB");
      std::shared_ptr<tbox::Database> mixB_db =
          calphad_db->getDatabase(namemixB);
      mixB_db->getDoubleArray("Left", &d_penalty_parameters[2][0], 3);
      mixB_db->getDoubleArray("Right", &d_penalty_parameters[2][3], 3);
   }
}

//-----------------------------------------------------------------------

double CALPHADFreeEnergyFunctionsWithPenaltyBinary::computeFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    const bool gp)
{
   double fe =
       CALPHADFreeEnergyFunctionsBinary::computeFreeEnergy(temperature, conc,
                                                           pi, false);

   double extra_energy = computePenalty(pi, conc[0]);

   // subtract -mu*c to get grand potential
   if (gp) {
      double deriv;
      computeDerivFreeEnergy(temperature, conc, pi, &deriv);
      fe -= deriv * conc[0];
   }

   return fe + extra_energy;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsWithPenaltyBinary::computeDerivFreeEnergy(
    const double temperature, const double* const conc, const PhaseIndex pi,
    double* deriv)
{
   double fe;
   CALPHADFreeEnergyFunctionsBinary::computeDerivFreeEnergy(temperature, conc,
                                                            pi, &fe);

   double extra_energy = computeDerivPenalty(pi, conc[0]);

   deriv[0] = fe + extra_energy;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsWithPenaltyBinary::
    computeSecondDerivativeFreeEnergy(const double temp,
                                      const std::vector<double>& conc,
                                      const PhaseIndex pi,
                                      std::vector<double>& d2fdc2)
{
   CALPHADFreeEnergyFunctionsBinary::computeSecondDerivativeFreeEnergy(temp,
                                                                       &conc[0],
                                                                       pi,
                                                                       d2fdc2);

   double extra_energy = compute2ndDerivPenalty(pi, conc[0]);

   d2fdc2[0] += extra_energy;
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyFunctionsWithPenaltyBinary::computeCeqT(
    const double temperature, const PhaseIndex pi0, const PhaseIndex pi1,
    double* ceq)
{
   assert(temperature > 0.);
   assert(d_penalty_parameters.size() > 0);

   int N = 2;

   double* fA = new double[N];
   double* fB = new double[N];
   double* L0 = new double[N];
   double* L1 = new double[N];
   double* L2 = new double[N];
   double* L3 = new double[N];

   setupValuesForTwoPhasesSolver(temperature, L0, L1, L2, L3, fA, fB, pi0, pi1);
   std::vector<std::vector<double> > penalty_parameters;
   penalty_parameters.push_back(d_penalty_parameters[static_cast<int>(pi0)]);
   penalty_parameters.push_back(d_penalty_parameters[static_cast<int>(pi1)]);
   double RTinv = 1.0 / (gas_constant_R_JpKpmol * temperature);
   CALPHADEqConcentrationSolverBinaryWithPenalty eq_solver;
   int ret =
       eq_solver.ComputeConcentrationWithPenalty(ceq, RTinv, L0, L1, L2, L3, fA,
                                                 fB, penalty_parameters);

   if (ret >= 0) {
      switch (pi0) {
         case PhaseIndex::phaseL:
            tbox::pout << "CALPHAD with Penalty, c_eq phaseL=" << ceq[0]
                       << std::endl;
            d_ceq_l = ceq[0];
            break;
         case PhaseIndex::phaseA:
            tbox::pout << "CALPHAD with Penalty, c_eq phaseA=" << ceq[0]
                       << std::endl;
            d_ceq_a = ceq[0];
            break;
         case PhaseIndex::phaseB:
            tbox::pout << "CALPHAD with Penalty, c_eq phaseB=" << ceq[0]
                       << std::endl;
            d_ceq_b = ceq[0];
            break;
      }

      switch (pi1) {
         case PhaseIndex::phaseL:
            tbox::pout << "CALPHAD with Penalty, c_eq phaseL=" << ceq[1]
                       << std::endl;
            d_ceq_l = ceq[1];
            break;
         case PhaseIndex::phaseA:
            tbox::pout << "CALPHAD with Penalty, c_eq phaseA=" << ceq[1]
                       << std::endl;
            d_ceq_a = ceq[1];
            break;
         case PhaseIndex::phaseB:
            tbox::pout << "CALPHAD with Penalty, c_eq phaseB=" << ceq[1]
                       << std::endl;
            d_ceq_b = ceq[1];
            break;
      }
   }

   delete[] fA;
   delete[] fB;
   delete[] L0;
   delete[] L1;
   delete[] L2;
   delete[] L3;

   return (ret >= 0);
}
