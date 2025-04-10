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
#include "ParabolicFreeEnergyMultiOrderBinaryThreePhase.h"
#include "ParabolicTools.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

ParabolicFreeEnergyMultiOrderBinaryThreePhase::
    ParabolicFreeEnergyMultiOrderBinaryThreePhase(
        std::shared_ptr<tbox::Database> input_db,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        const short norderp_A, MolarVolumeStrategy* mvstrategy,
        const int conc_l_id, const int conc_a_id, const int conc_b_id)
    : FreeEnergyStrategyBinary(Thermo4PFM::EnergyInterpolationType::LINEAR,
                               conc_interp_func_type, conc_l_id, conc_a_id,
                               conc_b_id, false),
      d_mv_strategy(mvstrategy)
{
   tbox::plog << "ParabolicFreeEnergyMultiOrderBinaryThreePhase..."
              << std::endl;

   assert(norderp_A > 0);
   assert(conc_b_id >= 0);

   double coeffL[3][2];
   readParabolicData(input_db, "Liquid", coeffL);

   double coeffA[3][2];
   readParabolicData(input_db, "PhaseA", coeffA);

   double coeffB[3][2];
   readParabolicData(input_db, "PhaseB", coeffB);

   double Tref = input_db->getDouble("Tref");

   d_parabolic_fenergy.reset(
       new Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase(
           Tref, coeffL, coeffA, coeffB,
           Thermo4PFM::EnergyInterpolationType::LINEAR,
           Thermo4PFM::ConcInterpolationType::LINEAR));

   d_multiorder_driving_force.reset(
       new MultiOrderBinaryThreePhasesDrivingForce(this, norderp_A));

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "ParabolicFreeEnergyMultiOrderBinaryThreePhase:" <<
   // std::endl; tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;
}

ParabolicFreeEnergyMultiOrderBinaryThreePhase::
    ~ParabolicFreeEnergyMultiOrderBinaryThreePhase(){};

bool ParabolicFreeEnergyMultiOrderBinaryThreePhase::computeCeqT(
    const double temperature, const Thermo4PFM::PhaseIndex pi0,
    const Thermo4PFM::PhaseIndex pi1, double* ceq)
{
   TBOX_ERROR(
       "ParabolicFreeEnergyMultiOrderBinaryThreePhase::computeCeqT() not "
       "implemented");
   return false;
}


//=======================================================================

double ParabolicFreeEnergyMultiOrderBinaryThreePhase::computeFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi,
    const bool gp)
{
   assert(d_mv_strategy != nullptr);

   double f = d_parabolic_fenergy->computeFreeEnergy(temperature, c_i, pi, gp);
   f *= d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
   return f;
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinaryThreePhase::computeDerivFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_parabolic_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   deriv *= d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
   return deriv;
}


void ParabolicFreeEnergyMultiOrderBinaryThreePhase::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   (void)time;
   (void)f_b_id;

   d_multiorder_driving_force->addDrivingForce(patch, temperature_id, phase_id,
                                               d_conc_l_id, d_conc_a_id,
                                               d_conc_b_id, f_l_id, f_a_id,
                                               f_b_id, rhs_id);
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinaryThreePhase::computeMuL(
    const double t, const double c0)
{
   double c = c0;
   double mu;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseL,
                                               &mu);
   return mu *
          d_mv_strategy->computeInvMolarVolume(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseL);
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinaryThreePhase::computeMuA(
    const double t, const double c0)
{
   double c = c0;
   double mu;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseA,
                                               &mu);
   return mu *
          d_mv_strategy->computeInvMolarVolume(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseA);
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinaryThreePhase::computeMuB(
    const double t, const double c0)
{
   double c = c0;
   double mu;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseB,
                                               &mu);
   return mu *
          d_mv_strategy->computeInvMolarVolume(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseB);
}

//=======================================================================

void ParabolicFreeEnergyMultiOrderBinaryThreePhase::
    computeSecondDerivativeEnergyPhaseL(const double temp,
                                        const std::vector<double>& c_l,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_parabolic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c_l[0], Thermo4PFM::PhaseIndex::phaseL, &d2fdc2[0]);
   if (use_internal_units)
      for (short i = 0; i < 3; i++)
         d2fdc2[i] *= d_mv_strategy->computeInvMolarVolume(
             temp, &c_l[0], Thermo4PFM::PhaseIndex::phaseL);
}

//=======================================================================

void ParabolicFreeEnergyMultiOrderBinaryThreePhase::
    computeSecondDerivativeEnergyPhaseA(const double temp,
                                        const std::vector<double>& c_a,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_parabolic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c_a[0], Thermo4PFM::PhaseIndex::phaseA, &d2fdc2[0]);
   if (use_internal_units)
      for (short i = 0; i < 3; i++)
         d2fdc2[i] *= d_mv_strategy->computeInvMolarVolume(
             temp, &c_a[0], Thermo4PFM::PhaseIndex::phaseA);
}

//=======================================================================

void ParabolicFreeEnergyMultiOrderBinaryThreePhase::
    computeSecondDerivativeEnergyPhaseB(const double temp,
                                        const std::vector<double>& c_b,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_parabolic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c_b[0], Thermo4PFM::PhaseIndex::phaseB, &d2fdc2[0]);
   if (use_internal_units)
      for (short i = 0; i < 3; i++)
         d2fdc2[i] *= d_mv_strategy->computeInvMolarVolume(
             temp, &c_b[0], Thermo4PFM::PhaseIndex::phaseB);
}
