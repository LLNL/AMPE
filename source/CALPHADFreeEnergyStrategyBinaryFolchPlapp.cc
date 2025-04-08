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
#include "CALPHADFreeEnergyStrategyBinaryFolchPlapp.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "MolarVolumeStrategy.h"
#include "TiltingFolchPlapp2005.h"

#include <boost/property_tree/json_parser.hpp>
#include "Database2JSON.h"
namespace pt = boost::property_tree;

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

template <class FreeEnergyFunctionType>
CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    CALPHADFreeEnergyStrategyBinaryFolchPlapp(
        pt::ptree calphad_db, std::shared_ptr<tbox::Database> newton_db,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        MolarVolumeStrategy* mvstrategy, const int conc_l_id,
        const int conc_a_id, const int conc_b_id)
    : ConcFreeEnergyStrategy(conc_l_id, conc_a_id, conc_b_id),
      d_energy_interp_func_type(energy_interp_func_type),
      d_conc_interp_func_type(conc_interp_func_type),
      d_mv_strategy(mvstrategy)
{
   assert(conc_l_id >= 0);
   assert(conc_a_id >= 0);
   assert(conc_b_id >= 0);

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   // d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1
   // tbox::plog << "CALPHADFreeEnergyStrategyBinaryFolchPlapp:" << std::endl;
   // tbox::plog << "Molar volume L =" << vml << std::endl;
   // tbox::plog << "Molar volume A =" << vma << std::endl;
   // tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << std::endl;

   setup(calphad_db, newton_db);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::setup(
    pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db)
{
   tbox::pout << "CALPHADFreeEnergyStrategyBinaryFolchPlapp<"
                 "FreeEnergyFunctionType>::setup()"
              << std::endl;
   // newton_db->printClassData(std::cout);
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   pt::json_parser::write_json(std::cout, newton_pt);

   d_calphad_fenergy.reset(new FreeEnergyFunctionType(calphad_pt, newton_pt,
                                                      d_energy_interp_func_type,
                                                      d_conc_interp_func_type));
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<
    FreeEnergyFunctionType>::computeFreeEnergy(hier::Patch& patch,
                                               const int temperature_id,
                                               const int f_id,
                                               const int conc_i_id,
                                               const Thermo4PFM::PhaseIndex pi,
                                               const bool gp)
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

   computeFreeEnergy(pbox, temperature, f, c_i, pi, gp);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    computeDerivFreeEnergy(hier::Patch& patch, const int temperature_id,
                           const int df_id, const int conc_i_id,
                           const Thermo4PFM::PhaseIndex pi)
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

   computeDerivFreeEnergy(pbox, temperature, df, c_i, pi);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    computeFreeEnergy(const hier::Box& pbox,
                      std::shared_ptr<pdat::CellData<double> > cd_temp,
                      std::shared_ptr<pdat::CellData<double> > cd_free_energy,
                      std::shared_ptr<pdat::CellData<double> > cd_conc_i,
                      const Thermo4PFM::PhaseIndex pi, const bool gp)
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
                d_calphad_fenergy->computeFreeEnergy(t, &c_i, pi, gp);
            ptr_f[idx_f] *= d_mv_strategy->computeInvMolarVolume(t, &c_i, pi);
         }
      }
   }
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    computeDerivFreeEnergy(
        const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
        std::shared_ptr<pdat::CellData<double> > cd_free_energy,
        std::shared_ptr<pdat::CellData<double> > cd_conc_i,
        const Thermo4PFM::PhaseIndex pi)
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

            d_calphad_fenergy->computeDerivFreeEnergy(t, &c_i, pi,
                                                      &ptr_f[idx_f]);
            ptr_f[idx_f] *= d_mv_strategy->computeInvMolarVolume(t, &c_i, pi);
         }
      }
   }
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<
    FreeEnergyFunctionType>::addDrivingForce(const double time,
                                             hier::Patch& patch,
                                             const int temperature_id,
                                             const int phase_id,
                                             const int conc_id,
                                             const int f_l_id, const int f_a_id,
                                             const int f_b_id, const int rhs_id)
{
   (void)time;

   assert(conc_id >= 0);
   assert(phase_id >= 0);
   assert(f_l_id >= 0);
   assert(f_a_id >= 0);
   assert(f_b_id >= 0);
   assert(rhs_id >= 0);
   assert(temperature_id >= 0);
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   assert(d_conc_b_id >= 0);

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_id)));
   assert(phase);
   assert(phase->getDepth() == 3);

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

   std::shared_ptr<pdat::CellData<double> > fb(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(f_b_id)));
   assert(fb);

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
   assert(rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   const hier::Box& pbox(patch.getBox());

   addDrivingForceOnPatch(rhs, t, phase, fl, fa, fb, c_l, c_a, c_b, pbox);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    addDrivingForceOnPatch(
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
   double* ptr_rhs0 = cd_rhs->getPointer(0);
   double* ptr_rhs1 = cd_rhs->getPointer(1);
   double* ptr_rhs2 = cd_rhs->getPointer(2);

   double* ptr_temp = cd_temperature->getPointer();
   double* ptr_phi0 = cd_phi->getPointer(0);
   double* ptr_phi1 = cd_phi->getPointer(1);
   double* ptr_phi2 = cd_phi->getPointer(2);
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

   // Assuming f_l, f_a, and f_b all have same ghost box
   assert(cd_f_l->getGhostCellWidth()[0] == cd_f_a->getGhostCellWidth()[0]);
   assert(cd_f_l->getGhostCellWidth()[0] == cd_f_b->getGhostCellWidth()[0]);
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

   // Assuming c_l, c_a, and c_b all have same ghost box
   assert(cd_c_l->getGhostCellWidth()[0] == cd_c_a->getGhostCellWidth()[0]);
   assert(cd_c_l->getGhostCellWidth()[0] == cd_c_b->getGhostCellWidth()[0]);
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

            double phi0 = ptr_phi0[idx_pf];
            double phi1 = ptr_phi1[idx_pf];
            double phi2 = ptr_phi2[idx_pf];

            double fl = ptr_f_l[idx_f_i];
            double fa = ptr_f_a[idx_f_i];
            double fb = ptr_f_b[idx_f_i];

            double cl = ptr_c_l[idx_c_i];
            double ca = ptr_c_a[idx_c_i];
            double cb = ptr_c_b[idx_c_i];

            double t = ptr_temp[idx_temp];
            double mu = computeMuA(t, ca);
            // double muL = computeMuL(t, cl);
            // double muB = computeMuB(t, cb);
            // if (fabs(mu - muL) > 1.e-2 * fabs(mu) ||
            //    fabs(muB - muL) > 1.e-2 * fabs(muB)) {
            //   double c = phi0 * cl + phi1 * ca + phi2 * cb;
            //   std::cout << "WARNING: c=" << c << ", ca=" << ca << ", cl=" <<
            //   cl
            //             << ", cb=" << cb << ", phi0=" << phi0
            //             << ", phi1=" << phi1 << ", phi2=" << phi2
            //             << ", muA=" << mu << ", muL=" << muL << ", muB=" <<
            //             muB
            //             << std::endl;
            //}
            double rhs0 = (TiltingFolchPlapp2005::dg0dp0(phi0, phi1, phi2) *
                               (fl - mu * cl) +
                           TiltingFolchPlapp2005::dg0dp1(phi1, phi0, phi2) *
                               (fa - mu * ca) +
                           TiltingFolchPlapp2005::dg0dp1(phi2, phi0, phi1) *
                               (fb - mu * cb));

            double rhs1 = (TiltingFolchPlapp2005::dg0dp0(phi1, phi0, phi2) *
                               (fa - mu * ca) +
                           TiltingFolchPlapp2005::dg0dp1(phi0, phi1, phi2) *
                               (fl - mu * cl) +
                           TiltingFolchPlapp2005::dg0dp1(phi2, phi1, phi0) *
                               (fb - mu * cb));

            double rhs2 = (TiltingFolchPlapp2005::dg0dp0(phi2, phi1, phi0) *
                               (fb - mu * cb) +
                           TiltingFolchPlapp2005::dg0dp1(phi1, phi2, phi0) *
                               (fa - mu * ca) +
                           TiltingFolchPlapp2005::dg0dp1(phi0, phi2, phi1) *
                               (fl - mu * cl));
            // add to rhs a component that satisfies the constraint on phi
            double corr = (rhs0 + rhs1 + rhs2) / 3.;
            ptr_rhs0[idx_rhs] -= (rhs0 - corr);
            ptr_rhs1[idx_rhs] -= (rhs1 - corr);
            ptr_rhs2[idx_rhs] -= (rhs2 - corr);
         }
      }
   }
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinaryFolchPlapp<
    FreeEnergyFunctionType>::computeMuA(const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseA,
                                             &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseA);

   return mu;
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinaryFolchPlapp<
    FreeEnergyFunctionType>::computeMuL(const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseL,
                                             &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseL);

   return mu;
}

//=======================================================================

template <class FreeEnergyFunctionType>
double CALPHADFreeEnergyStrategyBinaryFolchPlapp<
    FreeEnergyFunctionType>::computeMuB(const double t, const double c)
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t, &c,
                                             Thermo4PFM::PhaseIndex::phaseB,
                                             &mu);
   mu *= d_mv_strategy->computeInvMolarVolume(t, &c,
                                              Thermo4PFM::PhaseIndex::phaseB);

   return mu;
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    defaultComputeSecondDerivativeEnergyPhaseL(const double temp,
                                               const std::vector<double>& c_l,
                                               std::vector<double>& d2fdc2,
                                               const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c_l[0], Thermo4PFM::PhaseIndex::phaseL, d2fdc2.data());

   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temp, &c_l[0],
                                               Thermo4PFM::PhaseIndex::phaseL);
}

//=======================================================================

template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    defaultComputeSecondDerivativeEnergyPhaseA(const double temp,
                                               const std::vector<double>& c_a,
                                               std::vector<double>& d2fdc2,
                                               const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c_a[0], Thermo4PFM::PhaseIndex::phaseA, d2fdc2.data());

   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temp, &c_a[0],
                                               Thermo4PFM::PhaseIndex::phaseA);
}

//=======================================================================
template <class FreeEnergyFunctionType>
void CALPHADFreeEnergyStrategyBinaryFolchPlapp<FreeEnergyFunctionType>::
    defaultComputeSecondDerivativeEnergyPhaseB(const double temp,
                                               const std::vector<double>& c_b,
                                               std::vector<double>& d2fdc2,
                                               const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(
       temp, &c_b[0], Thermo4PFM::PhaseIndex::phaseB, d2fdc2.data());

   if (use_internal_units)
      d2fdc2[0] *=
          d_mv_strategy->computeInvMolarVolume(temp, &c_b[0],
                                               Thermo4PFM::PhaseIndex::phaseB);
}

template class CALPHADFreeEnergyStrategyBinaryFolchPlapp<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>;
template class CALPHADFreeEnergyStrategyBinaryFolchPlapp<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>;
