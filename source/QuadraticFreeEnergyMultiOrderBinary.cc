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
#include "QuadraticFreeEnergyMultiOrderBinary.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

QuadraticFreeEnergyMultiOrderBinary::QuadraticFreeEnergyMultiOrderBinary(
    std::shared_ptr<tbox::Database> input_db,
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    const double vml, const double vma, const int conc_l_id,
    const int conc_a_id)
    : FreeEnergyStrategyBinary(energy_interp_func_type, conc_interp_func_type,
                               conc_l_id, conc_a_id, -1, false),
      d_conc_l_id(conc_l_id),
      d_conc_a_id(conc_a_id)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);

   // From J/mol -> pJ/um**3
   d_energy_conv_factor_L = 1.e-6 / vml;
   d_energy_conv_factor_A = 1.e-6 / vma;

   double Tref = input_db->getDouble("T_ref");

   double A_liquid = input_db->getDouble("A_liquid");
   double Ceq_liquid = input_db->getDouble("Ceq_liquid");
   double m_liquid = input_db->getDouble("m_liquid");

   double A_solid_A = input_db->getDouble("A_solid");
   double Ceq_solid_A = input_db->getDouble("Ceq_solid");
   double m_solid = input_db->getDouble("m_solid");

   d_quadratic_fenergy.reset(new Thermo4PFM::QuadraticFreeEnergyFunctionsBinary(
       Tref, A_liquid, Ceq_liquid, m_liquid, A_solid_A, Ceq_solid_A, m_solid,
       energy_interp_func_type, Thermo4PFM::ConcInterpolationType::LINEAR));

   d_multiorder_driving_force.reset(new MultiOrderBinaryDrivingForce(this));

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "QuadraticFreeEnergyMultiOrderBinary:" << std::endl;
   // tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;
}

//=======================================================================

double QuadraticFreeEnergyMultiOrderBinary::computeFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi,
    const bool gp)
{
   double f = d_quadratic_fenergy->computeFreeEnergy(temperature, c_i, pi);
   return f * energyFactor(pi);
}

//=======================================================================

double QuadraticFreeEnergyMultiOrderBinary::computeDerivFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_quadratic_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   return deriv * energyFactor(pi);
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderBinary::addDrivingForce(
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

double QuadraticFreeEnergyMultiOrderBinary::computeMuL(const double t,
                                                       const double c_l)
{
   double deriv;
   double conc = c_l;
   d_quadratic_fenergy->computeDerivFreeEnergy(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseL,
                                               &deriv);

   return deriv * d_energy_conv_factor_L;
}

//=======================================================================

double QuadraticFreeEnergyMultiOrderBinary::computeMuA(const double t,
                                                       const double c_a)
{
   double deriv;
   double conc = c_a;
   d_quadratic_fenergy->computeDerivFreeEnergy(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseA,
                                               &deriv);

   return deriv * d_energy_conv_factor_A;
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderBinary::computeSecondDerivativeEnergyPhaseL(
    const double temperature, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temperature;

   double c = c_l[0];
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c, Thermo4PFM::PhaseIndex::phaseL, &d2fdc2[0]);
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_L;
}
//=======================================================================

void QuadraticFreeEnergyMultiOrderBinary::computeSecondDerivativeEnergyPhaseA(
    const double temperature, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temperature;
   double c = c_a[0];
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c, Thermo4PFM::PhaseIndex::phaseA, &d2fdc2[0]);
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_A;
}
