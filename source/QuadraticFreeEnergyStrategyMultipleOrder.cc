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
#ifdef HAVE_THERMO4PFM

#include "QuadraticFreeEnergyStrategyMultipleOrder.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

QuadraticFreeEnergyStrategyMultipleOrder::
    QuadraticFreeEnergyStrategyMultipleOrder(
        const int norderp, std::shared_ptr<tbox::Database> input_db,
        const EnergyInterpolationType energy_interp_func_type, const double vml,
        const double vma, const int conc_l_id, const int conc_a_id)
    : d_norderp(norderp),
      d_vm_L(vml),
      d_vm_A(vma),
      d_energy_interp_func_type(energy_interp_func_type),
      d_conc_l_id(conc_l_id),
      d_conc_a_id(conc_a_id)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);

   d_energy_conv_factor_L = 1.e-6 / d_vm_L;
   d_energy_conv_factor_A = 1.e-6 / d_vm_A;

   d_A_liquid = input_db->getDouble("A_liquid");
   d_Ceq_liquid = input_db->getDouble("Ceq_liquid");

   d_A_solid_A = input_db->getDouble("A_solid");
   d_Ceq_solid_A = input_db->getDouble("Ceq_solid");

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "QuadraticFreeEnergyStrategyMultipleOrder:" << std::endl;
   // tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;
}

//=======================================================================

double evaluateDerivFreeEnergy(const double conc, const double A,
                               const double Ceq, const double energy_factor)
{
   double mu = 2. * A * (conc - Ceq);

   return mu * energy_factor;
}

//=======================================================================

double evaluateFreeEnergy(const double conc, const double A, const double Ceq,
                          const double energy_factor, const bool gp = false)
{
   double d = conc - Ceq;

   double fe = A * d * d;
   fe *= energy_factor;

   // subtract -mu*c to get grand potential
   if (gp) fe -= evaluateDerivFreeEnergy(conc, A, Ceq, energy_factor) * conc;

   return fe;
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder ::computeDerivFreeEnergyLiquid(
    hier::Patch& patch, const int dfl_id)
{
   assert(dfl_id >= 0);

   computeDerivFreeEnergy(patch, d_A_liquid, d_Ceq_liquid, dfl_id, d_conc_l_id,
                          d_energy_conv_factor_L);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder ::computeDerivFreeEnergySolidA(
    hier::Patch& patch, const int dfa_id)
{
   assert(dfa_id >= 0);

   computeDerivFreeEnergy(patch, d_A_solid_A, d_Ceq_solid_A, dfa_id,
                          d_conc_a_id, d_energy_conv_factor_A);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder ::computeFreeEnergyLiquid(
    hier::Patch& patch, const int fl_id, const bool gp)
{
   assert(fl_id >= 0);

   assert(d_conc_l_id >= 0);

   computeFreeEnergy(patch, d_A_liquid, d_Ceq_liquid, fl_id, d_conc_l_id,
                     d_energy_conv_factor_L);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder ::computeFreeEnergySolidA(
    hier::Patch& patch, const int fa_id, const bool gp)
{
   assert(fa_id >= 0);

   computeFreeEnergy(patch, d_A_solid_A, d_Ceq_solid_A, fa_id, d_conc_a_id,
                     d_energy_conv_factor_A);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder ::computeFreeEnergy(
    hier::Patch& patch, const double A, const double Ceq, const int f_id,
    const int conc_i_id, const double energy_factor)
{
   assert(f_id >= 0);
   assert(conc_i_id >= 0);

   const hier::Box& pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<double> > f(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_id)));

   std::shared_ptr<pdat::CellData<double> > c_i(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_i_id)));

   computeFreeEnergy(pbox, A, Ceq, f, c_i, energy_factor);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder ::computeDerivFreeEnergy(
    hier::Patch& patch, const double A, const double Ceq, const int df_id,
    const int conc_i_id, const double energy_factor)
{
   assert(df_id >= 0);
   assert(conc_i_id >= 0);

   const hier::Box& pbox = patch.getBox();

   std::shared_ptr<pdat::CellData<double> > df(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(df_id)));

   std::shared_ptr<pdat::CellData<double> > c_i(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(conc_i_id)));

   computeDerivFreeEnergy(pbox, A, Ceq, df, c_i, energy_factor);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder::computeFreeEnergy(
    const hier::Box& pbox, const double A, const double Ceq,
    std::shared_ptr<pdat::CellData<double> > cd_free_energy,
    std::shared_ptr<pdat::CellData<double> > cd_conc_i,
    const double energy_factor)
{
   double* ptr_f = cd_free_energy->getPointer();
   double* ptr_c_i = cd_conc_i->getPointer();

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

            const int idx_f =
                (ii - imin_f) + (jj - jmin_f) * jp_f + (kk - kmin_f) * kp_f;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            double c_i = ptr_c_i[idx_c_i];

            ptr_f[idx_f] = evaluateFreeEnergy(c_i, A, Ceq, energy_factor);
         }
      }
   }
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder::computeDerivFreeEnergy(
    const hier::Box& pbox, const double A, const double Ceq,
    std::shared_ptr<pdat::CellData<double> > cd_free_energy,
    std::shared_ptr<pdat::CellData<double> > cd_conc_i,
    const double energy_factor)
{
   double* ptr_f = cd_free_energy->getPointer();
   double* ptr_c_i = cd_conc_i->getPointer();

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

            const int idx_f =
                (ii - imin_f) + (jj - jmin_f) * jp_f + (kk - kmin_f) * kp_f;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            double c_i = ptr_c_i[idx_c_i];

            ptr_f[idx_f] = evaluateFreeEnergy(c_i, A, Ceq, energy_factor);
         }
      }
   }
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder::addDrivingForce(
    const double time, hier::Patch& patch, const int phase_id,
    const int conc_id, const int f_l_id, const int f_a_id, const int rhs_id)
{
   (void)time;

   assert(conc_id >= 0);
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(rhs_id >= 0);
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);
   assert(phase->getDepth() == 3);

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

   const hier::Box& pbox(patch.getBox());

   addDrivingForceOnPatch(rhs, phase, fl, fa, c_l, c_a, pbox);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder::addDrivingForceOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_rhs,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_f_l,
    std::shared_ptr<pdat::CellData<double> > cd_f_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox)
{
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

   // Assuming f_l, f_a, all have same ghost box
   assert(cd_f_l->getGhostCellWidth()[0] == cd_f_a->getGhostCellWidth()[0]);
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

   // Assuming c_l, c_a, all have same ghost box
   assert(cd_c_l->getGhostCellWidth()[0] == cd_c_a->getGhostCellWidth()[0]);
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

   std::vector<double> rhs(d_norderp);

   std::vector<double*> ptr_rhs(d_norderp);
   for (short i = 0; i < d_norderp; i++)
      ptr_rhs[i] = cd_rhs->getPointer(i);

   std::vector<double*> ptr_phi(d_norderp);
   for (short i = 0; i < d_norderp; i++)
      ptr_phi[i] = cd_phi->getPointer(i);

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_rhs = (ii - imin_rhs) + (jj - jmin_rhs) * jp_rhs +
                                (kk - kmin_rhs) * kp_rhs;

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;

            const int idx_f_i = (ii - imin_f_i) + (jj - jmin_f_i) * jp_f_i +
                                (kk - kmin_f_i) * kp_f_i;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;

            double fl = ptr_f_l[idx_f_i];
            double fa = ptr_f_a[idx_f_i];

            double cl = ptr_c_l[idx_c_i];
            double ca = ptr_c_a[idx_c_i];

            double mu = computeMu(cl);

            //
            // see Moelans, Acta Mat 59 (2011)
            //

            // driving forces
            double dfs = (fa - mu * ca);
            double dfl = (fl - mu * cl);

            // interpolation polynomials
            double hphis = 0.;
            for (short i = 1; i < d_norderp; i++)
               hphis += ptr_phi[i][idx_pf] * ptr_phi[i][idx_pf];

            double hphil = ptr_phi[0][idx_pf] * ptr_phi[0][idx_pf];

            const double sum2 = hphil + hphis;
            const double sum2inv = 1. / sum2;

            hphis *= sum2inv;
            hphil *= sum2inv;

            // rhs from solid part
            rhs[0] = -2. * ptr_phi[0][idx_pf] * hphis * dfs * sum2inv;
            for (short i = 1; i < d_norderp; i++)
               rhs[i] = 2. * ptr_phi[i][idx_pf] * hphil * dfs * sum2inv;

            // rhs from liquid part
            rhs[0] = -2. * ptr_phi[0][idx_pf] * hphil * dfl * sum2inv;
            for (short i = 1; i < d_norderp; i++)
               rhs[i] = 2. * ptr_phi[i][idx_pf] * hphis * dfl * sum2inv;

            // add to rhs a component to satisfies the constraint on phi
            double corr = 0.;
            for (short i = 1; i < d_norderp; i++)
               corr += rhs[i];
            corr /= (double)d_norderp;

            for (short i = 1; i < d_norderp; i++)
               ptr_rhs[i][idx_rhs] -= (rhs[i] - corr);
         }
      }
   }
}

//=======================================================================

double QuadraticFreeEnergyStrategyMultipleOrder::computeMu(const double c_l)
{
   const double A = d_A_liquid;
   const double Ceq = d_Ceq_liquid;

   return evaluateDerivFreeEnergy(c_l, A, Ceq, d_energy_conv_factor_L);
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder::
    defaultComputeSecondDerivativeEnergyPhaseL(const std::vector<double>& c_l,
                                               std::vector<double>& d2fdc2,
                                               const bool use_internal_units)
{
   d2fdc2[0] = 2. * d_A_liquid;

   if (use_internal_units) d2fdc2[0] *= d2fdc2[0] *= d_energy_conv_factor_L;
}

//=======================================================================

void QuadraticFreeEnergyStrategyMultipleOrder::
    defaultComputeSecondDerivativeEnergyPhaseA(const std::vector<double>& c_a,
                                               std::vector<double>& d2fdc2,
                                               const bool use_internal_units)
{
   assert(c_a.size() == 1);
   d2fdc2[0] = 2. * d_A_solid_A;
   if (use_internal_units) d2fdc2[0] *= d_energy_conv_factor_A;
}

#endif
