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
#include "QuatParams.h"
#include "QuadraticFreeEnergyStrategy.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
using namespace SAMRAI;

#include <cassert>

#include <vector>


//=======================================================================

QuadraticFreeEnergyStrategy::QuadraticFreeEnergyStrategy(
    std::shared_ptr<tbox::Database> input_db,
    const EnergyInterpolationType energy_interp_func_type, const double vml,
    const double vma, const double vmb, const int conc_l_id,
    const int conc_a_id, const int conc_b_id, const bool with_third_phase)
{
   assert(conc_l_id >= 0);
   assert(conc_a_id >= 0);

   d_with_third_phase = with_third_phase;

   d_energy_interp_func_type = energy_interp_func_type;

   d_vm_L = vml;
   d_vm_A = vma;
   d_vm_B = vmb;

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // j/mol -> pj/mumcube = 1.e-6 / d_vm;

   d_energy_conv_factor_L = 1.e-6 / d_vm_L;
   d_energy_conv_factor_A = 1.e-6 / d_vm_A;
   if (d_with_third_phase) d_energy_conv_factor_B = 1.e-6 / d_vm_B;

   tbox::plog << "QuadraticFreeEnergyStrategy:" << std::endl;
   tbox::plog << "Molar volume L =" << d_vm_L << std::endl;
   tbox::plog << "Molar volume A =" << d_vm_A << std::endl;
   if (d_with_third_phase)
      tbox::plog << "Molar volume B =" << d_vm_B << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;

   d_A_liquid = input_db->getDouble("A_liquid");
   d_Ceq_liquid = input_db->getDouble("Ceq_liquid");

   if (!d_with_third_phase) {
      d_A_solid_A = input_db->getDouble("A_solid");
      d_Ceq_solid_A = input_db->getDouble("Ceq_solid");
      d_A_solid_B = 1.0;
      d_Ceq_solid_B = 1.0;
   } else {
      d_A_solid_A = input_db->getDouble("A_solid_A");
      d_A_solid_B = input_db->getDouble("A_solid_B");
      d_Ceq_solid_A = input_db->getDouble("Ceq_solid_A");
      d_Ceq_solid_B = input_db->getDouble("Ceq_solid_B");
   }

   // print database just read
   tbox::plog << "Quadratic database..." << std::endl;
   input_db->printClassData(tbox::plog);

   d_conc_l_id = conc_l_id;
   d_conc_a_id = conc_a_id;
   d_conc_b_id = conc_b_id;
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeFreeEnergyLiquid(
    hier::Patch& patch, const int temperature_id, const int fl_id,
    const bool gp)
{
   assert(fl_id >= 0);
   assert(temperature_id >= 0.);

   assert(d_conc_l_id >= 0);

   computeFreeEnergy(patch, temperature_id, d_A_liquid, d_Ceq_liquid, fl_id,
                     d_conc_l_id, d_energy_conv_factor_L);
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeFreeEnergySolidA(
    hier::Patch& patch, const int temperature_id, const int fs_id,
    const bool gp)
{
   assert(fs_id >= 0);
   assert(temperature_id >= 0.);

   assert(d_conc_a_id >= 0);

   computeFreeEnergy(patch, temperature_id, d_A_solid_A, d_Ceq_solid_A, fs_id,
                     d_conc_a_id, d_energy_conv_factor_A);
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeFreeEnergySolidB(
    hier::Patch& patch, const int temperature_id, const int fs_id,
    const bool gp)
{
   assert(fs_id >= 0);
   assert(temperature_id >= 0.);

   assert(d_conc_b_id >= 0);

   computeFreeEnergy(patch, temperature_id, d_A_solid_B, d_Ceq_solid_B, fs_id,
                     d_conc_b_id, d_energy_conv_factor_B);
}

//=======================================================================

// \/\/ No temperature dependence yet

double QuadraticFreeEnergyStrategy::computeFreeEnergy(
    const double temperature, const double conc, const double A,
    const double Ceq, const double energy_factor, const bool gp) const
{
   (void)temperature;

   double d = conc - Ceq;

   double fe = A * d * d;
   fe *= energy_factor;

   // subtract -mu*c to get grand potential
   if (gp)
      fe -= computeDerivFreeEnergy(temperature, conc, A, Ceq, energy_factor) *
            conc;

   return fe;
}

//=======================================================================

double QuadraticFreeEnergyStrategy::computeDerivFreeEnergy(
    const double temperature, const double conc, const double A,
    const double Ceq, const double energy_factor) const
{
   (void)temperature;

   double mu = 2. * A * (conc - Ceq);

   return mu * energy_factor;
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeFreeEnergy(
    hier::Patch& patch, const int temperature_id, const double A,
    const double Ceq, const int f_id, const int conc_i_id,
    const double energy_factor)
{
   assert(temperature_id >= 0);
   assert(f_id >= 0);
   assert(conc_i_id >= 0);

   const hier::Box& pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));

   std::shared_ptr<pdat::CellData<double> > f(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_id)));

   std::shared_ptr<pdat::CellData<double> > c_i(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_i_id)));

   computeFreeEnergy(pbox, temperature, A, Ceq, f, c_i, energy_factor);
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeDerivFreeEnergy(
    hier::Patch& patch, const int temperature_id, const double A,
    const double Ceq, const int df_id, const int conc_i_id,
    const double energy_factor)
{
   assert(temperature_id >= 0);
   assert(df_id >= 0);
   assert(conc_i_id >= 0);

   const hier::Box& pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));

   std::shared_ptr<pdat::CellData<double> > df(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(df_id)));

   std::shared_ptr<pdat::CellData<double> > c_i(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_i_id)));

   computeDerivFreeEnergy(pbox, temperature, A, Ceq, df, c_i, energy_factor);
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeFreeEnergy(
    const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
    const double A, const double Ceq,
    std::shared_ptr<pdat::CellData<double> > cd_free_energy,
    std::shared_ptr<pdat::CellData<double> > cd_conc_i,
    const double energy_factor)
{
   double* ptr_temp = cd_temp->getPointer();
   double* ptr_f = cd_free_energy->getPointer();
   double* ptr_c_i = cd_conc_i->getPointer();

   const hier::Box& temp_gbox = cd_temp->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   const hier::Box& f_gbox = cd_free_energy->getGhostBox();
   int imin_f = f_gbox.lower(0);
   int jmin_f = f_gbox.lower(1);
   int jp_f = f_gbox.numberCells(0);
   int kmin_f = 0;
   int kp_f = 0;
#if (NDIM == 3)
   kmin_f = f_gbox.lower(2);
   kp_f = jp_f * f_gbox.numberCells(1);
#endif

   const hier::Box& c_i_gbox = cd_conc_i->getGhostBox();
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

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_f =
                (ii - imin_f) + (jj - jmin_f) * jp_f + (kk - kmin_f) * kp_f;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            double t = ptr_temp[idx_temp];

            double c_i = ptr_c_i[idx_c_i];

            ptr_f[idx_f] = computeFreeEnergy(t, c_i, A, Ceq, energy_factor);
         }
      }
   }
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeDerivFreeEnergy(
    const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
    const double A, const double Ceq,
    std::shared_ptr<pdat::CellData<double> > cd_free_energy,
    std::shared_ptr<pdat::CellData<double> > cd_conc_i,
    const double energy_factor)
{
   double* ptr_temp = cd_temp->getPointer();
   double* ptr_f = cd_free_energy->getPointer();
   double* ptr_c_i = cd_conc_i->getPointer();

   const hier::Box& temp_gbox = cd_temp->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   const hier::Box& f_gbox = cd_free_energy->getGhostBox();
   int imin_f = f_gbox.lower(0);
   int jmin_f = f_gbox.lower(1);
   int jp_f = f_gbox.numberCells(0);
   int kmin_f = 0;
   int kp_f = 0;
#if (NDIM == 3)
   kmin_f = f_gbox.lower(2);
   kp_f = jp_f * f_gbox.numberCells(1);
#endif

   const hier::Box& c_i_gbox = cd_conc_i->getGhostBox();
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

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_f =
                (ii - imin_f) + (jj - jmin_f) * jp_f + (kk - kmin_f) * kp_f;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            double t = ptr_temp[idx_temp];

            double c_i = ptr_c_i[idx_c_i];

            ptr_f[idx_f] =
                computeDerivFreeEnergy(t, c_i, A, Ceq, energy_factor);
         }
      }
   }
}

//=======================================================================

void QuadraticFreeEnergyStrategy::addDrivingForce(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int eta_id, const int conc_id, const int f_l_id,
    const int f_a_id, const int f_b_id, const int rhs_id)
{
   (void)time;

   assert(conc_id >= 0);
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);
   if (d_with_third_phase) {
      assert(eta_id >= 0);
      assert(f_b_id >= 0);
   }

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

   std::shared_ptr<pdat::CellData<double> > eta;
   std::shared_ptr<pdat::CellData<double> > fb;
   std::shared_ptr<pdat::CellData<double> > c_b;
   if (d_with_third_phase) {
      eta = std::dynamic_pointer_cast<pdat::CellData<double>, hier::PatchData>(
          patch.getPatchData(eta_id));
      assert(eta);
      fb = std::dynamic_pointer_cast<pdat::CellData<double>, hier::PatchData>(
          patch.getPatchData(f_b_id));
      assert(fb);

      c_b = std::dynamic_pointer_cast<pdat::CellData<double>, hier::PatchData>(
          patch.getPatchData(d_conc_b_id));
      assert(c_b);
   }

   const hier::Box& pbox = patch.getBox();

   addDrivingForceOnPatch(rhs, t, phase, eta, fl, fa, fb, c_l, c_a, c_b, pbox);
}

//=======================================================================

void QuadraticFreeEnergyStrategy::addDrivingForceOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_rhs,
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_eta,
    std::shared_ptr<pdat::CellData<double> > cd_f_l,
    std::shared_ptr<pdat::CellData<double> > cd_f_a,
    std::shared_ptr<pdat::CellData<double> > cd_f_b,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox)
{
   double* ptr_rhs = cd_rhs->getPointer();
   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_phi = cd_phi->getPointer();
   double* ptr_eta = NULL;
   double* ptr_f_l = cd_f_l->getPointer();
   double* ptr_f_a = cd_f_a->getPointer();
   double* ptr_f_b = NULL;
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();
   double* ptr_c_b = NULL;

   if (d_with_third_phase) {
      ptr_eta = cd_eta->getPointer();
      ptr_f_b = cd_f_b->getPointer();
      ptr_c_b = cd_c_b->getPointer();
   }

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

   // Assuming phi, eta, and concentration all have same box
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

   const char interp = energyInterpChar(d_energy_interp_func_type);

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
            double eta = 0.0;
            double f_l = ptr_f_l[idx_f_i];
            double f_a = ptr_f_a[idx_f_i];
            double f_b = 0.0;
            double c_l = ptr_c_l[idx_c_i];
            double c_a = ptr_c_a[idx_c_i];
            double c_b = 0.0;

            double mu = computeMu(t, c_l);

            double hphi_prime = DERIV_INTERP_FUNC(phi, &interp);

            double heta = 0.0;

            if (d_with_third_phase) {
               eta = ptr_eta[idx_pf];
               f_b = ptr_f_b[idx_f_i];
               c_b = ptr_c_b[idx_c_i];

               heta = INTERP_FUNC(eta, &interp);
            }

            ptr_rhs[idx_rhs] +=
                hphi_prime * ((f_l - (1.0 - heta) * f_a - heta * f_b) -
                              mu * (c_l - (1.0 - heta) * c_a - heta * c_b));
         }
      }
   }
}

//=======================================================================

double QuadraticFreeEnergyStrategy::computeMu(const double t, const double c_l)
{
   const double A = d_A_liquid;
   const double Ceq = d_Ceq_liquid;

   double mu = computeDerivFreeEnergy(t, c_l, A, Ceq, d_energy_conv_factor_L);

   return mu;
}

//=======================================================================

void QuadraticFreeEnergyStrategy::addDrivingForceEta(
    const double time, hier::Patch& patch, const int temperature_id,
    const int phase_id, const int eta_id, const int conc_id, const int f_l_id,
    const int f_a_id, const int f_b_id, const int rhs_id)
{
   (void)time;

   assert(conc_id >= 0);
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);
   assert(eta_id >= 0);
   assert(f_b_id >= 0);

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);

   std::shared_ptr<pdat::CellData<double> > eta(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(eta_id)));
   assert(eta);

   std::shared_ptr<pdat::CellData<double> > conc(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_id)));
   assert(conc);

   std::shared_ptr<pdat::CellData<double> > t(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(t);

   std::shared_ptr<pdat::CellData<double> > f_l(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_l_id)));
   assert(f_l);

   std::shared_ptr<pdat::CellData<double> > f_a(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_a_id)));
   assert(f_a);

   std::shared_ptr<pdat::CellData<double> > f_b(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_b_id)));
   assert(f_b);

   std::shared_ptr<pdat::CellData<double> > c_l(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_l_id)));
   assert(c_l);

   std::shared_ptr<pdat::CellData<double> > c_a(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_a_id)));
   assert(c_a);

   std::shared_ptr<pdat::CellData<double> > c_b(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_b_id)));
   assert(c_b);

   std::shared_ptr<pdat::CellData<double> > rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(rhs_id)));
   assert(rhs);

   const hier::Box& pbox = patch.getBox();

   addDrivingForceEtaOnPatch(rhs, t, phase, eta, f_l, f_a, f_b, c_l, c_a, c_b,
                             pbox);
}

//=======================================================================

void QuadraticFreeEnergyStrategy::addDrivingForceEtaOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_rhs,
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_eta,
    std::shared_ptr<pdat::CellData<double> > cd_f_l,
    std::shared_ptr<pdat::CellData<double> > cd_f_a,
    std::shared_ptr<pdat::CellData<double> > cd_f_b,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox)
{
   double* ptr_rhs = cd_rhs->getPointer();
   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_phi = cd_phi->getPointer();
   double* ptr_eta = cd_eta->getPointer();
   double* ptr_f_a = cd_f_a->getPointer();
   double* ptr_f_b = cd_f_b->getPointer();
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();
   double* ptr_c_b = cd_c_b->getPointer();

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

   // Assuming phi, eta, and concentration all have same box
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

   const char interpf = energyInterpChar(d_energy_interp_func_type);

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
            double eta = ptr_eta[idx_pf];
            double f_a = ptr_f_a[idx_f_i];
            double f_b = ptr_f_b[idx_f_i];
            double c_l = ptr_c_l[idx_c_i];
            double c_a = ptr_c_a[idx_c_i];
            double c_b = ptr_c_b[idx_c_i];

            double mu = computeMu(t, c_l);

            double hphi = INTERP_FUNC(phi, &interpf);

            double heta_prime = DERIV_INTERP_FUNC(eta, &interpf);

            ptr_rhs[idx_rhs] +=
                hphi * heta_prime * ((f_a - f_b) - mu * (c_a - c_b));
         }
      }
   }
}

//=======================================================================
#if 0
double QuadraticFreeEnergyStrategy::computeLocalInvD2fDc2(
   const double c,
   const double hphi,
   const double heta,
   const double temp )
{
   // delta c for finite difference calculation of d2f/dc2
   const double deltac = 0.001;

   double c_l =
      computeLiquidConcentration(
            hphi,heta,c-0.5*deltac,
            d_A_liquid, d_A_solid_A, d_A_solid_B,
            d_Ceq_liquid, d_Ceq_solid_A, d_Ceq_solid_B);
   const double mu_minus = computeMu( temp, c_l );
   
   c_l =
      computeLiquidConcentration(
            hphi,heta,c+0.5*deltac,
            d_A_liquid, d_A_solid_A, d_A_solid_B,
            d_Ceq_liquid, d_Ceq_solid_A, d_Ceq_solid_B);
   const double mu_plus = computeMu( temp, c_l );
   
   // (d2f/dc2)^-1
   return deltac/(mu_plus-mu_minus);

}
#endif

//=======================================================================

double QuadraticFreeEnergyStrategy::computeLiquidConcentration(
    const double hphi, const double heta, const double c) const
{
   return (c -
           hphi * (1.0 - heta) *
               (d_Ceq_solid_A - (d_A_liquid / d_A_solid_A) * d_Ceq_liquid) -
           hphi * heta *
               (d_Ceq_solid_B - (d_A_liquid / d_A_solid_B) * d_Ceq_liquid)) /
          ((1.0 - hphi) + hphi * (1.0 - heta) * (d_A_liquid / d_A_solid_A) +
           hphi * heta * (d_A_liquid / d_A_solid_B));
}

//=======================================================================

double QuadraticFreeEnergyStrategy::computeSolidAConcentration(
    const double hphi, const double heta, const double c) const
{
   return (c -
           (1.0 - hphi) *
               (d_Ceq_liquid - (d_A_solid_A / d_A_liquid) * d_Ceq_solid_A) -
           hphi * heta *
               (d_Ceq_solid_B - (d_A_solid_A / d_A_solid_B) * d_Ceq_solid_A)) /
          ((1.0 - hphi) * (d_A_solid_A / d_A_liquid) + hphi * (1.0 - heta) +
           hphi * heta * (d_A_solid_A / d_A_solid_B));
}

//=======================================================================

double QuadraticFreeEnergyStrategy::computeSolidBConcentration(
    const double hphi, const double heta, const double c) const
{
   return (c -
           (1.0 - hphi) *
               (d_Ceq_liquid - (d_A_solid_B / d_A_liquid) * d_Ceq_solid_B) -
           hphi * (1.0 - heta) *
               (d_Ceq_solid_A - (d_A_solid_B / d_A_solid_A) * d_Ceq_solid_B)) /
          ((1.0 - hphi) * (d_A_solid_B / d_A_liquid) +
           hphi * (1.0 - heta) * (d_A_solid_B / d_A_solid_A) + hphi * heta);
}

//=======================================================================

double QuadraticFreeEnergyStrategy::computeLiquidConcentration(
    const double hphi, const double heta, const double c, const double Al,
    const double Aa, const double Ab, const double Ceql, const double CeqA,
    const double CeqB) const
{
   return (c - hphi * (1.0 - heta) * (CeqA - (Al / Aa) * Ceql) -
           hphi * heta * (CeqB - (Al / Ab) * Ceql)) /
          ((1.0 - hphi) + hphi * (1.0 - heta) * (Al / Aa) +
           hphi * heta * (Al / Ab));
}

//=======================================================================

double QuadraticFreeEnergyStrategy::computeSolidAConcentration(
    const double hphi, const double heta, const double c, const double Al,
    const double Aa, const double Ab, const double Ceql, const double CeqA,
    const double CeqB) const
{
   return (c - (1.0 - hphi) * (Ceql - (Aa / Al) * CeqA) -
           hphi * heta * (CeqB - (Aa / Ab) * CeqA)) /
          ((1.0 - hphi) * (Aa / Al) + hphi * (1.0 - heta) +
           hphi * heta * (Aa / Ab));
}

//=======================================================================

double QuadraticFreeEnergyStrategy::computeSolidBConcentration(
    const double hphi, const double heta, const double c, const double Al,
    const double Aa, const double Ab, const double Ceql, const double CeqA,
    const double CeqB) const
{
   return (c - (1.0 - hphi) * (Ceql - (Ab / Al) * CeqB) -
           hphi * (1.0 - heta) * (CeqA - (Ab / Aa) * CeqB)) /
          ((1.0 - hphi) * (Ab / Al) + hphi * (1.0 - heta) * (Ab / Aa) +
           hphi * heta);
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseL(
    const double temp, const std::vector<double>& c_l,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;

   assert(c_l.size() == 1);
   d2fdc2[0] = 2. * d_A_liquid;
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_L;
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseA(
    const double temp, const std::vector<double>& c_a,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;

   assert(c_a.size() == 1);
   d2fdc2[0] = 2. * d_A_solid_A;
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_A;
}

//=======================================================================

void QuadraticFreeEnergyStrategy::computeSecondDerivativeEnergyPhaseB(
    const double temp, const std::vector<double>& c_b,
    std::vector<double>& d2fdc2, const bool use_internal_units)
{
   (void)temp;

   assert(c_b.size() == 1);
   d2fdc2[0] = 2. * d_A_solid_B;
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_B;
}
