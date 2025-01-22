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
#include "QuadraticFreeEnergyMultiOrderBinaryThreePhase.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

QuadraticFreeEnergyMultiOrderBinaryThreePhase::
    QuadraticFreeEnergyMultiOrderBinaryThreePhase(
        std::shared_ptr<tbox::Database> input_db,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        const short norderp_A, const double vml, const double vma,
        const double vmb, const int conc_l_id, const int conc_a_id,
        const int conc_b_id)
    : FreeEnergyStrategyBinary(energy_interp_func_type, conc_interp_func_type,
                               conc_l_id, conc_a_id, conc_b_id, false),
      d_norderp_A(norderp_A)
{
   tbox::plog << "QuadraticFreeEnergyMultiOrderBinaryThreePhase..."
              << std::endl;

   d_energy_conv_factor_L = 1.e-6 / vml;
   d_energy_conv_factor_A = 1.e-6 / vma;
   d_energy_conv_factor_B = 1.e-6 / vmb;

   double A_liquid = input_db->getDouble("A_liquid");
   double Ceq_liquid = input_db->getDouble("Ceq_liquid");

   double A_solid_A = input_db->getDouble("A_solid_A");
   double Ceq_solid_A = input_db->getDouble("Ceq_solid_A");

   double A_solid_B = input_db->getDouble("A_solid_B");
   double Ceq_solid_B = input_db->getDouble("Ceq_solid_B");

   d_quadratic_fenergy.reset(
       new Thermo4PFM::QuadraticFreeEnergyFunctionsBinaryThreePhase(
           A_liquid, Ceq_liquid, A_solid_A, Ceq_solid_A, A_solid_B, Ceq_solid_B,
           energy_interp_func_type, Thermo4PFM::ConcInterpolationType::LINEAR));

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "QuadraticFreeEnergyMultiOrderBinaryThreePhase:" <<
   // std::endl; tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;

   d_driving_force.reset(
       new MultiOrderBinaryThreePhasesDrivingForce(this, norderp_A));
}

QuadraticFreeEnergyMultiOrderBinaryThreePhase::
    ~QuadraticFreeEnergyMultiOrderBinaryThreePhase(){};

//=======================================================================

double QuadraticFreeEnergyMultiOrderBinaryThreePhase::computeFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi,
    const bool gp)
{

   double f = d_quadratic_fenergy->computeFreeEnergy(temperature, c_i, pi);
   return f * energyFactor(pi);
}

double QuadraticFreeEnergyMultiOrderBinaryThreePhase::computeDerivFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_quadratic_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   return deriv * energyFactor(pi);
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderBinaryThreePhase::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   (void)time;

   d_driving_force->addDrivingForce(patch, temperature_id, phase_id,
                                    d_conc_l_id, d_conc_a_id, d_conc_b_id,
                                    f_l_id, f_a_id, f_b_id, rhs_id);
}

//=======================================================================

double QuadraticFreeEnergyMultiOrderBinaryThreePhase::computeMuL(
    const double t, const double c0)
{
   double mu;
   double c = c0;
   d_quadratic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseL,
                                               &mu);
   return mu * d_energy_conv_factor_L;
}

//=======================================================================

double QuadraticFreeEnergyMultiOrderBinaryThreePhase::computeMuA(
    const double t, const double c0)
{
   double mu;
   double c = c0;
   d_quadratic_fenergy->computeDerivFreeEnergy(t, &c,
                                               Thermo4PFM::PhaseIndex::phaseA,
                                               &mu);
   return mu * d_energy_conv_factor_A;
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderBinaryThreePhase::
    computeSecondDerivativeEnergyPhaseL(const double temp,
                                        const std::vector<double>& c_l,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c_l[0], Thermo4PFM::PhaseIndex::phaseL, &d2fdc2[0]);
   if (use_internal_units)
      for (short i = 0; i < 3; i++)
         d2fdc2[i] *= d_energy_conv_factor_L;
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderBinaryThreePhase::
    computeSecondDerivativeEnergyPhaseA(const double temp,
                                        const std::vector<double>& c_a,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c_a[0], Thermo4PFM::PhaseIndex::phaseA, &d2fdc2[0]);
   if (use_internal_units)
      for (short i = 0; i < 3; i++)
         d2fdc2[i] *= d_energy_conv_factor_A;
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderBinaryThreePhase::
    computeSecondDerivativeEnergyPhaseB(const double temp,
                                        const std::vector<double>& c_b,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c_b[0], Thermo4PFM::PhaseIndex::phaseB, &d2fdc2[0]);
   if (use_internal_units)
      for (short i = 0; i < 3; i++)
         d2fdc2[i] *= d_energy_conv_factor_B;
}
