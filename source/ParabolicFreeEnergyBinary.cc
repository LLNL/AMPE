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
#include "ConcFort.h"
#include "ParabolicFreeEnergyBinary.h"
#include "MolarVolumeStrategy.h"
#include "ParabolicFreeEnergyFunctionsBinary.h"
#include "FuncFort.h"


#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;


#include <cassert>

//=======================================================================

ParabolicFreeEnergyBinary::ParabolicFreeEnergyBinary(
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    MolarVolumeStrategy* mvstrategy, const int conc_l_id, const int conc_a_id,
    std::shared_ptr<tbox::Database> conc_db)
    : FreeEnergyStrategyBinary(energy_interp_func_type, conc_interp_func_type,
                               conc_l_id, conc_a_id, -1, false),
      d_mv_strategy(mvstrategy)
{
   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "ParabolicFreeEnergyBinary:" << std::endl;
   // tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;
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

   d_parabolic_fenergy.reset(new Thermo4PFM::ParabolicFreeEnergyFunctionsBinary(
       Tref, coeffL, coeffA, energy_interp_func_type, conc_interp_func_type));
}

//=======================================================================
bool ParabolicFreeEnergyBinary::computeCeqT(const double temperature,
                                            const Thermo4PFM::PhaseIndex pi0,
                                            const Thermo4PFM::PhaseIndex pi1,
                                            double* ceq)
{
   TBOX_ERROR("computeCeqT not implented for that class");
   return false;
}

//=======================================================================

double ParabolicFreeEnergyBinary ::computeFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi,
    const bool gp)
{
   double f = d_parabolic_fenergy->computeFreeEnergy(temperature, c_i, pi, gp);
   f *= d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
   return f;
}

//=======================================================================

double ParabolicFreeEnergyBinary::computeDerivFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_parabolic_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   deriv *= d_mv_strategy->computeInvMolarVolume(temperature, c_i, pi);
   return deriv;
}

//=======================================================================

void ParabolicFreeEnergyBinary::computeDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type_saved(
       d_energy_interp_func_type);
   // use linear interpolation function to get driving force
   // without polynomial of phi factor
   d_energy_interp_func_type = Thermo4PFM::EnergyInterpolationType::LINEAR,

   FreeEnergyStrategy::computeDrivingForce(time, patch, temperature_id,
                                           phase_id, conc_id, f_l_id, f_a_id,
                                           f_b_id, rhs_id);

   d_energy_interp_func_type = d_energy_interp_func_type_saved;
};

//=======================================================================

double ParabolicFreeEnergyBinary::computeMuA(const double t, const double c)
{
   double mu;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseA,
                                               &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseA);

   return mu;
}

//=======================================================================

double ParabolicFreeEnergyBinary::computeMuL(const double t, const double c)
{
   double mu;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseL,
                                               &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseL);

   return mu;
}

double ParabolicFreeEnergyBinary::computeMuB(const double t, const double c)
{
   double mu;
   d_parabolic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseB,
                                               &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseB);

   return mu;
}


//=======================================================================

void ParabolicFreeEnergyBinary::computeSecondDerivativeEnergyPhaseL(
    const double temp, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   d_parabolic_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c_l[0], Thermo4PFM::PhaseIndex::phaseL, d2fdc2.data());

   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temp, &c_l[0],
                                               Thermo4PFM::PhaseIndex::phaseL);
}

//=======================================================================

void ParabolicFreeEnergyBinary::computeSecondDerivativeEnergyPhaseA(
    const double temp, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   d_parabolic_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c_a[0], Thermo4PFM::PhaseIndex::phaseA, d2fdc2.data());

   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temp, &c_a[0],
                                               Thermo4PFM::PhaseIndex::phaseA);
}
