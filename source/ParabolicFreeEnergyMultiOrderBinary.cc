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
#include "ParabolicFreeEnergyMultiOrderBinary.h"
#include "ParabolicTools.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

ParabolicFreeEnergyMultiOrderBinary::ParabolicFreeEnergyMultiOrderBinary(
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    MolarVolumeStrategy* mvstrategy, const int conc_l_id, const int conc_a_id,
    std::shared_ptr<tbox::Database> conc_db)
    : FreeEnergyStrategyBinary(energy_interp_func_type, conc_interp_func_type,
                               conc_l_id, conc_a_id, -1, false),
      d_mv_strategy(mvstrategy)
{
   std::shared_ptr<tbox::Database> input_db = conc_db->getDatabase("Parabolic");

   double coeffL[3][2];
   readParabolicData(input_db, "Liquid", coeffL);

   double coeffA[3][2];
   std::shared_ptr<tbox::Database> phasea_db = input_db->getDatabase("PhaseA");
   coeffA[0][0] = phasea_db->getDouble("a0");
   coeffA[0][1] = phasea_db->getDouble("a1");
   coeffA[1][0] = phasea_db->getDouble("b0");
   coeffA[1][1] = phasea_db->getDouble("b1");
   coeffA[2][0] = phasea_db->getDouble("c0");
   coeffA[2][1] = phasea_db->getDouble("c1");

   double Tref = input_db->getDouble("Tref");

   d_parabolic_fenergy.reset(new Thermo4PFM::ParabolicFreeEnergyFunctionsBinary(
       Tref, coeffL, coeffA, energy_interp_func_type, conc_interp_func_type));

   d_multiorder_driving_force.reset(new MultiOrderBinaryDrivingForce(this));
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinary::computeFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi,
    const bool gp)
{
   double f = d_parabolic_fenergy->computeFreeEnergy(temperature, c_i, pi);
   return f * d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinary::computeDerivFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_parabolic_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   return deriv * d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
}

//=======================================================================

void ParabolicFreeEnergyMultiOrderBinary::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   (void)time;
   (void)f_b_id;

   d_multiorder_driving_force->addDrivingForce(patch, temperature_id, phase_id,
                                               d_conc_l_id, d_conc_a_id, f_l_id,
                                               f_a_id, rhs_id);
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinary::computeMuL(const double t,
                                                       const double c_l)
{
   double deriv;
   double conc = c_l;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseL,
                                               &deriv);

   return deriv *
          d_mv_strategy->computeInvMolarVolume(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseL);
}

//=======================================================================

double ParabolicFreeEnergyMultiOrderBinary::computeMuA(const double t,
                                                       const double c_a)
{
   double deriv;
   double conc = c_a;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseA,
                                               &deriv);

   return deriv *
          d_mv_strategy->computeInvMolarVolume(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseA);
}

//=======================================================================

void ParabolicFreeEnergyMultiOrderBinary::computeSecondDerivativeEnergyPhaseL(
    const double temperature, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temperature;

   double c = c_l[0];
   d_parabolic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c, Thermo4PFM::PhaseIndex::phaseL, &d2fdc2[0]);
   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temperature, &c_l[0],
                                               Thermo4PFM::PhaseIndex::phaseL);
}
//=======================================================================

void ParabolicFreeEnergyMultiOrderBinary::computeSecondDerivativeEnergyPhaseA(
    const double temperature, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temperature;
   double c = c_a[0];
   d_parabolic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c, Thermo4PFM::PhaseIndex::phaseA, &d2fdc2[0]);
   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temperature, &c_a[0],
                                               Thermo4PFM::PhaseIndex::phaseA);
}
