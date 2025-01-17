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
#include "SAMRAI/tbox/InputManager.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Index.h"

#include "FuncFort.h"
#include "ConcFort.h"
#include "QuadraticFreeEnergyBinary.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
using namespace SAMRAI;

#include <cassert>

#include <vector>


//=======================================================================

QuadraticFreeEnergyBinary::QuadraticFreeEnergyBinary(
    std::shared_ptr<tbox::Database> input_db,
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    const double vml, const double vma, const int conc_l_id,
    const int conc_a_id)
    : FreeEnergyStrategyBinary(energy_interp_func_type, conc_interp_func_type,
                               conc_l_id, conc_a_id, -1, false)
{
   d_energy_interp_func_type = energy_interp_func_type;

   d_vm_L = vml;
   d_vm_A = vma;

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // j/mol -> pj/mumcube = 1.e-6 / d_vm;

   d_energy_conv_factor_L = 1.e-6 / d_vm_L;
   d_energy_conv_factor_A = 1.e-6 / d_vm_A;

   tbox::plog << "QuadraticFreeEnergyBinary:" << std::endl;
   tbox::plog << "Molar volume L =" << d_vm_L << std::endl;
   tbox::plog << "Molar volume A =" << d_vm_A << std::endl;

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

   // print database just read
   tbox::plog << "Quadratic database..." << std::endl;
   input_db->printClassData(tbox::plog);
}

//=======================================================================

double QuadraticFreeEnergyBinary::computeFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi,
    const bool gp)
{

   double f = d_quadratic_fenergy->computeFreeEnergy(temperature, c_i, pi);
   return f * energyFactor(pi);
}

//=======================================================================

double QuadraticFreeEnergyBinary::computeDerivFreeEnergy(
    const double temperature, double* c_i, const Thermo4PFM::PhaseIndex pi)
{
   double deriv;
   d_quadratic_fenergy->computeDerivFreeEnergy(temperature, c_i, pi, &deriv);
   return deriv * energyFactor(pi);
}

//=======================================================================

void QuadraticFreeEnergyBinary::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int conc_id, const int f_l_id, const int f_a_id,
    const int f_b_id, const int rhs_id)
{
   (void)time;
   (void)f_b_id;

   assert(conc_id >= 0);
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   std::shared_ptr<pdat::CellData<double> > t(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(t);

   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_l_id)));
   assert(fl);

   std::shared_ptr<pdat::CellData<double> > fa(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_a_id)));
   assert(fa);

   std::shared_ptr<pdat::CellData<double> > c_l(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_l_id)));
   assert(c_l);

   std::shared_ptr<pdat::CellData<double> > c_a(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_a_id)));
   assert(c_a);

   std::shared_ptr<pdat::CellData<double> > rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));

   assert(rhs);
   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   const hier::Box& pbox = patch.getBox();

   addDrivingForceOnPatch(rhs, t, phase, fl, fa, c_l, c_a, pbox);
}

//=======================================================================

void QuadraticFreeEnergyBinary::addDrivingForceOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_rhs,
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_f_l,
    std::shared_ptr<pdat::CellData<double> > cd_f_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox)
{
   double* ptr_rhs = cd_rhs->getPointer();
   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_phi = cd_phi->getPointer();
   double* ptr_f_l = cd_f_l->getPointer();
   double* ptr_f_a = cd_f_a->getPointer();
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();

   const hier::Box& rhs_gbox = cd_rhs->getGhostBox();
   int imin_rhs = rhs_gbox.lower(0);
   int jmin_rhs = rhs_gbox.lower(1);
   int jp_rhs = rhs_gbox.numberCells(0);
   int kmin_rhs = 0;
   int kp_rhs = 0;
#if (NDIM == 3)
   kmin_rhs = rhs_gbox.lower(2);
   kp_rhs = jp_rhs * rhs_gbox.numberCells(1);
#endif

   const hier::Box& temp_gbox = cd_temperature->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   // Assuming phi and concentration all have same box
   const hier::Box& pf_gbox = cd_phi->getGhostBox();
   int imin_pf = pf_gbox.lower(0);
   int jmin_pf = pf_gbox.lower(1);
   int jp_pf = pf_gbox.numberCells(0);
   int kmin_pf = 0;
   int kp_pf = 0;
#if (NDIM == 3)
   kmin_pf = pf_gbox.lower(2);
   kp_pf = jp_pf * pf_gbox.numberCells(1);
#endif

   // Assuming f_l, f_a, and f_b all have same box
   const hier::Box& f_i_gbox = cd_f_l->getGhostBox();
   int imin_f_i = f_i_gbox.lower(0);
   int jmin_f_i = f_i_gbox.lower(1);
   int jp_f_i = f_i_gbox.numberCells(0);
   int kmin_f_i = 0;
   int kp_f_i = 0;
#if (NDIM == 3)
   kmin_f_i = f_i_gbox.lower(2);
   kp_f_i = jp_f_i * f_i_gbox.numberCells(1);
#endif

   // Assuming c_l, c_a, and c_b all have same box
   const hier::Box& c_i_gbox = cd_c_l->getGhostBox();
   int imin_c_i = c_i_gbox.lower(0);
   int jmin_c_i = c_i_gbox.lower(1);
   int jp_c_i = c_i_gbox.numberCells(0);
   int kmin_c_i = 0;
   int kp_c_i = 0;
#if (NDIM == 3)
   kmin_c_i = c_i_gbox.lower(2);
   kp_c_i = jp_c_i * c_i_gbox.numberCells(1);
#endif

   int imin = pbox.lower(0);
   int imax = pbox.upper(0);
   int jmin = pbox.lower(1);
   int jmax = pbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = pbox.lower(2);
   kmax = pbox.upper(2);
#endif

   const char interp = Thermo4PFM::energyInterpChar(d_energy_interp_func_type);

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_rhs = (ii - imin_rhs) + (jj - jmin_rhs) * jp_rhs +
                                (kk - kmin_rhs) * kp_rhs;

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;

            const int idx_f_i = (ii - imin_f_i) + (jj - jmin_f_i) * jp_f_i +
                                (kk - kmin_f_i) * kp_f_i;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            double t = ptr_temp[idx_temp];
            double phi = ptr_phi[idx_pf];
            double f_l = ptr_f_l[idx_f_i];
            double f_a = ptr_f_a[idx_f_i];
            double c_l = ptr_c_l[idx_c_i];
            double c_a = ptr_c_a[idx_c_i];

            double mu = computeMuL(t, c_l);

            double hphi_prime = DERIV_INTERP_FUNC(phi, &interp);

            ptr_rhs[idx_rhs] += hphi_prime * ((f_l - f_a) - mu * (c_l - c_a));
         }
      }
   }
}

//=======================================================================

double QuadraticFreeEnergyBinary::computeMuL(const double t, const double c_l)
{
   double deriv;
   double conc = c_l;
   d_quadratic_fenergy->computeDerivFreeEnergy(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseL,
                                               &deriv);

   return deriv * d_energy_conv_factor_L;
}

double QuadraticFreeEnergyBinary::computeMuA(const double t, const double c_a)
{
   double deriv;
   double conc = c_a;
   d_quadratic_fenergy->computeDerivFreeEnergy(t, &conc,
                                               Thermo4PFM::PhaseIndex::phaseA,
                                               &deriv);

   return deriv * d_energy_conv_factor_A;
}

//=======================================================================

void QuadraticFreeEnergyBinary::computeSecondDerivativeEnergyPhaseL(
    const double temp, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   assert(c_l.size() == 1);

   double c = c_l[0];
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c, Thermo4PFM::PhaseIndex::phaseL, &d2fdc2[0]);
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_L;
}

//=======================================================================

void QuadraticFreeEnergyBinary::computeSecondDerivativeEnergyPhaseA(
    const double temp, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   assert(c_a.size() == 1);

   double c = c_a[0];
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c, Thermo4PFM::PhaseIndex::phaseA, &d2fdc2[0]);
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_A;
}
