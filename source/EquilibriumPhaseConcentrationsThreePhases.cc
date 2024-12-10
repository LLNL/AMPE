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
#include "Database2JSON.h"
namespace pt = boost::property_tree;

#include "EquilibriumPhaseConcentrationsThreePhases.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/tbox/IEEE.h"

#include <omp.h>

template <class FreeEnergyType>
EquilibriumPhaseConcentrationsThreePhases<FreeEnergyType>::
    EquilibriumPhaseConcentrationsThreePhases(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
}

template <>
EquilibriumPhaseConcentrationsThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>::
    EquilibriumPhaseConcentrationsThreePhases(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy =
       std::unique_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>(
           new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(
               calphad_pt, newton_pt, energy_interp_func_type,
               Thermo4PFM::ConcInterpolationType::LINEAR));
}

template <>
EquilibriumPhaseConcentrationsThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>::
    EquilibriumPhaseConcentrationsThreePhases(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy =
       std::unique_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>(
           new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl(
               calphad_pt, newton_pt, energy_interp_func_type,
               Thermo4PFM::ConcInterpolationType::LINEAR));
}

template <>
EquilibriumPhaseConcentrationsThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>::
    EquilibriumPhaseConcentrationsThreePhases(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy =
       std::unique_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>(
           new Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase(
               calphad_pt, newton_pt, energy_interp_func_type,
               Thermo4PFM::ConcInterpolationType::LINEAR));
}

template <>
EquilibriumPhaseConcentrationsThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>::
    EquilibriumPhaseConcentrationsThreePhases(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy =
       std::unique_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>(
           new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl(
               calphad_pt, newton_pt, energy_interp_func_type,
               Thermo4PFM::ConcInterpolationType::LINEAR));
}

template <>
EquilibriumPhaseConcentrationsThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>::
    EquilibriumPhaseConcentrationsThreePhases(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy =
       std::unique_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>(
           new Thermo4PFM::CALPHADFreeEnergyFunctionsTernary(
               calphad_pt, newton_pt, energy_interp_func_type,
               Thermo4PFM::ConcInterpolationType::LINEAR));
}


template <class FreeEnergyType>
int EquilibriumPhaseConcentrationsThreePhases<
    FreeEnergyType>::computeAuxilliaryConcentrations(const double temp,
                                                     double* c, double* hphi,
                                                     double* x)
{
   double x_init[4] = {x[0], x[1], x[2], x[3]};

   int status = d_calphad_fenergy->computePhaseConcentrations(temp, c, hphi, x);

#ifndef GPU_OFFLOAD
   if (status < 0 || std::isnan(x[0])) {
      std::cerr << "EquilibriumPhaseConcentrationsThreePhases" << std::endl;
      std::cerr << "computePhaseConcentrations failed for T = " << temp
                << ", hphi = ";
      for (short i = 0; i < 3; i++)
         std::cerr << hphi[i] << ", ";
      std::cerr << "c =" << c[0] << std::endl;
      std::cerr << "cinit = " << x_init[0] << "," << x_init[1] << ","
                << x_init[2] << std::endl;
      std::cerr << "x = " << x[0] << "," << x[1] << "," << x[2] << std::endl;
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      MPI_Abort(mpi.getCommunicator(), -1);
   }
#endif
   return status;
}


template <class FreeEnergyType>
int EquilibriumPhaseConcentrationsThreePhases<FreeEnergyType>::
    computePhaseConcentrationsOnPatch(
        std::shared_ptr<pdat::CellData<double> > cd_te,
        std::shared_ptr<pdat::CellData<double> > cd_pf,
        std::shared_ptr<pdat::CellData<double> > cd_conc,
        std::shared_ptr<pdat::CellData<double> > cd_cl,
        std::shared_ptr<pdat::CellData<double> > cd_ca,
        std::shared_ptr<pdat::CellData<double> > cd_cb,
        std::shared_ptr<hier::Patch> patch)
{
   assert(cd_te);
   assert(cd_pf);
   assert(cd_conc);
   assert(cd_cl);
   assert(cd_ca);
   assert(cd_conc->getDepth() == cd_cl->getDepth());
   assert(cd_conc->getDepth() == cd_ca->getDepth());
   assert(cd_cl->getGhostCellWidth()[0] <= cd_te->getGhostCellWidth()[0]);
   assert(cd_cl->getGhostCellWidth()[0] <= cd_pf->getGhostCellWidth()[0]);
#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> cops;
   double l2n = cops.L2Norm(cd_conc, patch->getBox());
   assert(l2n == l2n);
#endif

   const int nphases = cd_pf->getDepth();
   if (nphases == 3) assert(cd_cb);

   std::shared_ptr<pdat::CellData<double> > cd_cl_ref(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_conc_l_ref_id)));
   assert(cd_cl_ref);
   assert(cd_conc->getDepth() == cd_cl_ref->getDepth());

   std::shared_ptr<pdat::CellData<double> > cd_ca_ref(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_conc_a_ref_id)));
   assert(cd_ca_ref);
   assert(cd_conc->getDepth() == cd_ca_ref->getDepth());

   std::shared_ptr<pdat::CellData<double> > cd_cb_ref;
   if (nphases == 3) {
      cd_cb_ref =
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_conc_b_ref_id));
      assert(cd_cb_ref);
   }

   const double* const ptr_temp = cd_te->getPointer();
   double* ptr_phi[3];  // up to 3 phases
   for (int i = 0; i < nphases; i++)
      ptr_phi[i] = cd_pf->getPointer(i);

   const hier::Box& temp_gbox = cd_te->getGhostBox();
   int imin_te = temp_gbox.lower(0);
   int jmin_te = temp_gbox.lower(1);
   int inc_j_te = temp_gbox.numberCells(0);
   int kmin_te = 0;
   int inc_k_te = 0;
#if (NDIM == 3)
   kmin_te = temp_gbox.lower(2);
   inc_k_te = inc_j_te * temp_gbox.numberCells(1);
#endif

   // Assuming phi and concentration all have same box
   assert(cd_pf->getGhostCellWidth()[0] == cd_conc->getGhostCellWidth()[0]);

   const hier::Box& pf_gbox = cd_pf->getGhostBox();
   int imin_pf = pf_gbox.lower(0);
   int jmin_pf = pf_gbox.lower(1);
   int inc_j_pf = pf_gbox.numberCells(0);
   int kmin_pf = 0;
   int inc_k_pf = 0;
#if (NDIM == 3)
   kmin_pf = pf_gbox.lower(2);
   inc_k_pf = inc_j_pf * pf_gbox.numberCells(1);
#endif

   // Assuming c_l, c_a, and c_b all have same box
   assert(cd_cl->getGhostCellWidth()[0] == cd_ca->getGhostCellWidth()[0]);
   assert(cd_cl->getGhostCellWidth()[0] == cd_cl_ref->getGhostCellWidth()[0]);
   const hier::Box& ci_gbox = cd_cl->getGhostBox();
   int imin_ci = ci_gbox.lower(0);
   int jmin_ci = ci_gbox.lower(1);
   int inc_j_ci = ci_gbox.numberCells(0);
   int kmin_ci = 0;
   int inc_k_ci = 0;
#if (NDIM == 3)
   kmin_ci = ci_gbox.lower(2);
   inc_k_ci = inc_j_ci * ci_gbox.numberCells(1);
#endif

   // loop indexes are based on internal compositions
   // ghost boxes since we need to initialize their ghost values
   int imin[3] = {ci_gbox.lower(0), ci_gbox.lower(1), 0};
   int imax[3] = {ci_gbox.upper(0), ci_gbox.upper(1), 0};
#if (NDIM == 3)
   imin[2] = ci_gbox.lower(2);
   imax[2] = ci_gbox.upper(2);
#endif

   // number of compositions fields (number of species -1)
   const int nc = cd_conc->getDepth();
   int N = 2 * cd_conc->getDepth();
   if (nphases == 3) {
      N += cd_conc->getDepth();
   }

   double* cl = cd_cl->getPointer(0);
   double* ca = cd_ca->getPointer(0);
   double* cb = nullptr;
   double* conc = cd_conc->getPointer(0);
   double* cl_ref = cd_cl_ref->getPointer(0);
   double* ca_ref = cd_ca_ref->getPointer(0);
   double* cb_ref = nullptr;
   if (nphases == 3) {
      cb = cd_cb->getPointer(0);
      cb_ref = cd_cb_ref->getPointer(0);
   }

   // number of cells for each field
   const int ncp = pf_gbox.size();
   const int ncc = ci_gbox.size();

   // Limits computation to one ghost cell if the Patch touches a non-periodic
   // physical boundary.
   if (cd_cl->getGhostCellWidth()[0] > 1) {
      // get CartesianPatchGeometry associated with patch
      std::shared_ptr<geom::CartesianPatchGeometry> pg(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                 hier::PatchGeometry>(
              patch->getPatchGeometry()));

      // The side is specified by an axis direction (direction normal to the
      // side being checked) and a flag specified the upper (side=1) or lower
      // side (side=0).
      if (pg->getTouchesRegularBoundary(0, 0)) {
         imin[0]++;
      }
      if (pg->getTouchesRegularBoundary(0, 1)) {
         imax[0]--;
      }
      if (pg->getTouchesRegularBoundary(1, 0)) {
         imin[1]++;
      }
      if (pg->getTouchesRegularBoundary(1, 1)) {
         imax[1]--;
      }
      if (pg->getTouchesRegularBoundary(2, 0)) {
         imin[2]++;
      }
      if (pg->getTouchesRegularBoundary(2, 1)) {
         imax[2]--;
      }
   }

   int nits = 0;

#ifdef GPU_OFFLOAD
   const int nct = temp_gbox.size();
// clang-format off
#pragma omp target map(to: ptr_temp[:nct]) \
                   map(to: ptr_phi[:ncp]) \
                   map(to: conc[:nc*ncp]) \
                   map(to : cl_ref[:nc*ncc]) \
                   map(to : ca_ref[:nc*ncc]) \
                   map(from : cl[:nc*ncc]) \
                   map(from : ca[:nc*ncc])
   // clang-format on
   {
#pragma omp teams distribute
#endif
      for (int kk = imin[2]; kk <= imax[2]; kk++) {
#ifdef GPU_OFFLOAD
#pragma omp parallel for collapse(2) schedule(static, 1)
#endif
         for (int jj = imin[1]; jj <= imax[1]; jj++) {
            for (int ii = imin[0]; ii <= imax[0]; ii++) {

               int idx_ci = (ii - imin_ci) + (jj - jmin_ci) * inc_j_ci +
                            (kk - kmin_ci) * inc_k_ci;
               int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * inc_j_pf +
                            (kk - kmin_pf) * inc_k_pf;
               int idx_te = (ii - imin_te) + (jj - jmin_te) * inc_j_te +
                            (kk - kmin_te) * inc_k_te;

               const double temp = ptr_temp[idx_te];
               // up to 3 phases
               double phi[3];
               for (short i = 0; i < nphases; i++)
                  phi[i] = ptr_phi[i][idx_pf];
               double hphi[3];
               for (short i = 0; i < nphases; i++)
                  hphi[i] =
                      Thermo4PFM::interp_func(d_conc_interp_func_type, phi[i]);
               double c[2];  // up to 2 components
               for (int ic = 0; ic < nc; ic++) {
                  c[ic] = conc[ic * ncp + idx_pf];
               }
               double x[4];
               for (int ic = 0; ic < nc; ic++) {
                  x[ic] = cl_ref[ic * ncc + idx_ci];
                  x[ic + nc] = ca_ref[ic * ncc + idx_ci];
                  if (nphases == 3) {
                     x[ic + 2 * nc] = cb_ref[ic * ncc + idx_ci];
                  }
               }
               assert(!std::isnan(x[0]));
               assert(!std::isnan(x[nc]));
               // std::cerr << ", c=" << c[0] << std::endl;
               // std::cerr << "c_ref=" << cl_ref[idx_ci] << "," <<
               // ca_ref[idx_ci]
               //          << "," << cb_ref[idx_ci];
               // std::cerr << ", hphi =";
               // for (short i = 0; i < nphases; i++)
               //   std::cerr << hphi[i] << ", ";

               // compute cL, cS
               int status = computeAuxilliaryConcentrations(temp, c, hphi, x);

               nits += status;

               // std::cout << "phi=" << phi[0] << "," << phi[1] << "," <<
               // phi[2]
               //          << "c=" << c[0] << ", x=" << x[0] << "," << x[1] <<
               //          ","
               //          << x[2] << std::endl;

               // set cell values with cL and cS just computed
               for (int ic = 0; ic < nc; ic++) {
                  cl[ic * ncc + idx_ci] = x[ic];
                  ca[ic * ncc + idx_ci] = x[ic + nc];
                  if (nphases == 3) {
                     cb[ic * ncc + idx_ci] = x[ic + 2 * nc];
                  }

               }  // ic

            }  // ii

         }  // jj

      }  // kk
#ifdef GPU_OFFLOAD
   }
#endif
   return nits;
}

template class EquilibriumPhaseConcentrationsThreePhases<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>;
