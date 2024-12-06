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
#include "QuadraticFreeEnergyMultiOrderTernaryThreePhase.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

QuadraticFreeEnergyMultiOrderTernaryThreePhase::
    QuadraticFreeEnergyMultiOrderTernaryThreePhase(
        std::shared_ptr<tbox::Database> input_db,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const short norderp_A, const double vml, const double vma,
        const double vmb, const int conc_l_id, const int conc_a_id,
        const int conc_b_id)
    : FreeEnergyStrategyThreePhase(input_db, vml, vma, vmb, conc_l_id,
                                   conc_a_id, conc_b_id),
      d_norderp_A(norderp_A)
{
   tbox::plog << "QuadraticFreeEnergyMultiOrderTernaryThreePhase..."
              << std::endl;

   double A_liquid[2];
   input_db->getDoubleArray("A_liquid", &A_liquid[0], 2);
   double Ceq_liquid[2];
   input_db->getDoubleArray("Ceq_liquid", &Ceq_liquid[0], 2);

   double A_solid_A[2];
   input_db->getDoubleArray("A_solid_A", &A_solid_A[0], 2);
   double Ceq_solid_A[2];
   input_db->getDoubleArray("Ceq_solid_A", &Ceq_solid_A[0], 2);

   double A_solid_B[2];
   input_db->getDoubleArray("A_solid_B", &A_solid_B[0], 2);
   double Ceq_solid_B[2];
   input_db->getDoubleArray("Ceq_solid_B", &Ceq_solid_B[0], 2);

   d_quadratic_fenergy.reset(
       new Thermo4PFM::QuadraticFreeEnergyFunctionsTernaryThreePhase(
           A_liquid, Ceq_liquid, A_solid_A, Ceq_solid_A, A_solid_B, Ceq_solid_B,
           energy_interp_func_type, Thermo4PFM::ConcInterpolationType::LINEAR));

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "QuadraticFreeEnergyMultiOrderTernaryThreePhase:" <<
   // std::endl; tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;
}

//=======================================================================

QuadraticFreeEnergyMultiOrderTernaryThreePhase::
    ~QuadraticFreeEnergyMultiOrderTernaryThreePhase()
{
}


void QuadraticFreeEnergyMultiOrderTernaryThreePhase::computeFreeEnergy(
    const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
    std::shared_ptr<pdat::CellData<double> > cd_free_energy,
    std::shared_ptr<pdat::CellData<double> > cd_conc_i,
    Thermo4PFM::PhaseIndex pi, const double energy_factor)
{
   assert(cd_conc_i->getDepth() == 2);

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

   const size_t coffset = c_i_gbox.size();

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

            // get local two components of ternary alloy
            double c_i[2] = {ptr_c_i[idx_c_i], ptr_c_i[idx_c_i + coffset]};

            ptr_f[idx_f] = d_quadratic_fenergy->computeFreeEnergy(t, c_i, pi);
            ptr_f[idx_f] *= energy_factor;
         }
      }
   }
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderTernaryThreePhase::addDrivingForceOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_rhs,
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_f_l,
    std::shared_ptr<pdat::CellData<double> > cd_f_a,
    std::shared_ptr<pdat::CellData<double> > cd_f_b,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox)
{
   assert(cd_temperature);
   assert(cd_f_b);
   assert(cd_c_b);
   assert(cd_f_l->getGhostCellWidth()[0] == cd_f_a->getGhostCellWidth()[0]);
   assert(cd_f_l->getGhostCellWidth()[0] == cd_f_b->getGhostCellWidth()[0]);
   assert(cd_c_l->getGhostCellWidth()[0] == cd_c_a->getGhostCellWidth()[0]);
   assert(cd_c_l->getGhostCellWidth()[0] == cd_c_b->getGhostCellWidth()[0]);
   assert(cd_phi->getDepth() > 1);
   assert(cd_phi->getDepth() == cd_rhs->getDepth());

   const int norderp = cd_phi->getDepth();

   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_f_l = cd_f_l->getPointer();
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

   std::vector<double> rhs(norderp);

   std::vector<double*> ptr_rhs(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_rhs[i] = cd_rhs->getPointer(i);

   std::vector<double*> ptr_phi(norderp);
   for (short i = 0; i < norderp; i++)
      ptr_phi[i] = cd_phi->getPointer(i);

   // number of cells for each field
   const size_t coffset = c_i_gbox.size();
   // tbox::plog<<"d_norderp_A = "<<d_norderp_A<<std::endl;

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
            double fl = ptr_f_l[idx_f_i];
            double fa = ptr_f_a[idx_f_i];
            double fb = ptr_f_b[idx_f_i];
            double cl[2] = {ptr_c_l[idx_c_i], ptr_c_l[idx_c_i + coffset]};
            double ca[2] = {ptr_c_a[idx_c_i], ptr_c_a[idx_c_i + coffset]};
            double cb[2] = {ptr_c_b[idx_c_i], ptr_c_b[idx_c_i + coffset]};

            double muA[2];
            computeMuA(t, ca[0], ca[1], muA);
            double muB[2];
            computeMuB(t, cb[0], cb[1], muB);
            double muL[2];
            computeMuL(t, cl[0], cl[1], muL);

            //
            // see Moelans, Acta Mat 59 (2011)
            //

            // driving forces
            double dfa = (fa - muA[0] * ca[0] - muA[1] * ca[1]);
            double dfb = (fb - muB[0] * cb[0] - muB[1] * cb[1]);
            double dfl = (fl - muL[0] * cl[0] - muL[1] * cl[1]);
            assert(!std::isnan(dfa));

            // interpolation polynomials
            double hphiA = 0.;
            for (short i = 0; i < d_norderp_A; i++)
               hphiA += ptr_phi[i][idx_pf] * ptr_phi[i][idx_pf];
            assert(!std::isnan(hphiA));

            double hphiB = 0.;
            for (short i = d_norderp_A; i < norderp - 1; i++)
               hphiB += ptr_phi[i][idx_pf] * ptr_phi[i][idx_pf];
            assert(!std::isnan(hphiB));

            double hphil =
                ptr_phi[norderp - 1][idx_pf] * ptr_phi[norderp - 1][idx_pf];

            const double sum2 = hphil + hphiA + hphiB;
            assert(sum2 > 0.);
            const double sum2inv = 1. / sum2;

            hphil *= sum2inv;
            hphiA *= sum2inv;
            hphiB *= sum2inv;

            assert(!std::isnan(hphiA));
            assert(!std::isnan(hphiB));

            // solid phase A order parameters
            for (short i = 0; i < d_norderp_A; i++)
               rhs[i] = 2. * ptr_phi[i][idx_pf] *
                        ((1. - hphiA) * dfa - hphil * dfl - hphiB * dfb) *
                        sum2inv;
            // solid phase B order parameters
            for (short i = d_norderp_A; i < norderp - 1; i++)
               rhs[i] = 2. * ptr_phi[i][idx_pf] *
                        ((1. - hphiB) * dfb - hphil * dfl - hphiA * dfa) *
                        sum2inv;
            // liquid phase order parameter
            rhs[norderp - 1] =
                2. * ptr_phi[norderp - 1][idx_pf] *
                ((1. - hphil) * dfl - hphiA * dfa - hphiB * dfb) * sum2inv;
            for (short i = 0; i < norderp; i++)
               assert(!std::isnan(rhs[i]));

            for (short i = 0; i < norderp; i++)
               ptr_rhs[i][idx_rhs] -= (rhs[i]);
         }
      }
   }
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderTernaryThreePhase::computeMuL(const double t,
                                                                const double c0,
                                                                const double c1,
                                                                double* mu)
{
   double c[2] = {c0, c1};
   d_quadratic_fenergy->computeDerivFreeEnergy(t, c,
                                               Thermo4PFM::PhaseIndex::phaseL,
                                               mu);
   mu[0] *= d_energy_conv_factor_L;
   mu[1] *= d_energy_conv_factor_L;
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderTernaryThreePhase::computeMuA(const double t,
                                                                const double c0,
                                                                const double c1,
                                                                double* mu)
{
   double c[2] = {c0, c1};
   d_quadratic_fenergy->computeDerivFreeEnergy(t, c,
                                               Thermo4PFM::PhaseIndex::phaseA,
                                               mu);
   mu[0] *= d_energy_conv_factor_A;
   mu[1] *= d_energy_conv_factor_A;
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderTernaryThreePhase::computeMuB(const double t,
                                                                const double c0,
                                                                const double c1,
                                                                double* mu)
{
   double c[2] = {c0, c1};
   d_quadratic_fenergy->computeDerivFreeEnergy(t, c,
                                               Thermo4PFM::PhaseIndex::phaseB,
                                               mu);
   mu[0] *= d_energy_conv_factor_B;
   mu[1] *= d_energy_conv_factor_B;
}

//=======================================================================

void QuadraticFreeEnergyMultiOrderTernaryThreePhase::
    computeSecondDerivativeEnergyPhaseL(const std::vector<double>& c_l,
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

void QuadraticFreeEnergyMultiOrderTernaryThreePhase::
    computeSecondDerivativeEnergyPhaseA(const std::vector<double>& c_a,
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

void QuadraticFreeEnergyMultiOrderTernaryThreePhase::
    computeSecondDerivativeEnergyPhaseB(const std::vector<double>& c_b,
                                        std::vector<double>& d2fdc2,
                                        const bool use_internal_units)
{
   d_quadratic_fenergy->computeSecondDerivativeFreeEnergy(
       0., &c_b[0], Thermo4PFM::PhaseIndex::phaseB, &d2fdc2[0]);
   if (use_internal_units)
      for (short i = 0; i < 3; i++)
         d2fdc2[i] *= d_energy_conv_factor_B;
}
