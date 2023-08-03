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
#include "Database2JSON.h"
namespace pt = boost::property_tree;
#endif

#include "CALPHADequilibriumPhaseConcentrationsStrategy.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#ifdef HAVE_THERMO4PFM
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#else
#include "FuncFort.h"
#endif

#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/tbox/IEEE.h"

#include <omp.h>

template <>
CALPHADequilibriumPhaseConcentrationsStrategy<
    CALPHADFreeEnergyFunctionsBinary>::
    CALPHADequilibriumPhaseConcentrationsStrategy(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        const bool with_third_phase,
#ifdef HAVE_THERMO4PFM
        pt::ptree calphad_pt,
#else
        std::shared_ptr<tbox::Database> calphad_db,
#endif
        std::shared_ptr<tbox::Database> newton_db, const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id, with_third_phase),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
#ifdef HAVE_THERMO4PFM
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
#endif
   d_calphad_fenergy = std::unique_ptr<CALPHADFreeEnergyFunctionsBinary>(
       new CALPHADFreeEnergyFunctionsBinary(
#ifdef HAVE_THERMO4PFM
           calphad_pt, newton_pt, energy_interp_func_type,
           ConcInterpolationType::LINEAR
#else
           calphad_db, newton_db, energy_interp_func_type,
           ConcInterpolationType::LINEAR, with_third_phase
#endif
           ));
}

#ifdef HAVE_THERMO4PFM
template <>
CALPHADequilibriumPhaseConcentrationsStrategy<
    CALPHADFreeEnergyFunctionsBinary2Ph1Sl>::
    CALPHADequilibriumPhaseConcentrationsStrategy(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        const bool with_third_phase, pt::ptree calphad_pt,
        std::shared_ptr<tbox::Database> newton_db, const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id, with_third_phase),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy = std::unique_ptr<CALPHADFreeEnergyFunctionsBinary2Ph1Sl>(
       new CALPHADFreeEnergyFunctionsBinary2Ph1Sl(
           calphad_pt, newton_pt, energy_interp_func_type,
           ConcInterpolationType::LINEAR));
}

template <>
CALPHADequilibriumPhaseConcentrationsStrategy<
    CALPHADFreeEnergyFunctionsBinaryThreePhase>::
    CALPHADequilibriumPhaseConcentrationsStrategy(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        const bool with_third_phase, pt::ptree calphad_pt,
        std::shared_ptr<tbox::Database> newton_db, const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id, with_third_phase),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy =
       std::unique_ptr<CALPHADFreeEnergyFunctionsBinaryThreePhase>(
           new CALPHADFreeEnergyFunctionsBinaryThreePhase(
               calphad_pt, newton_pt, energy_interp_func_type,
               ConcInterpolationType::LINEAR));
}

template <>
CALPHADequilibriumPhaseConcentrationsStrategy<
    CALPHADFreeEnergyFunctionsBinary3Ph2Sl>::
    CALPHADequilibriumPhaseConcentrationsStrategy(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        const bool with_third_phase, pt::ptree calphad_pt,
        std::shared_ptr<tbox::Database> newton_db, const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id, with_third_phase),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy = std::unique_ptr<CALPHADFreeEnergyFunctionsBinary3Ph2Sl>(
       new CALPHADFreeEnergyFunctionsBinary3Ph2Sl(
           calphad_pt, newton_pt, energy_interp_func_type,
           ConcInterpolationType::LINEAR));
}

#endif

template <>
CALPHADequilibriumPhaseConcentrationsStrategy<
    CALPHADFreeEnergyFunctionsTernary>::
    CALPHADequilibriumPhaseConcentrationsStrategy(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        const bool with_third_phase,
#ifdef HAVE_THERMO4PFM
        pt::ptree calphad_pt,
#else
        std::shared_ptr<tbox::Database> calphad_db,
#endif
        std::shared_ptr<tbox::Database> newton_db, const unsigned ncompositions)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id, with_third_phase),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id),
      d_conc_interp_func_type(conc_interp_func_type)
{
#ifdef HAVE_THERMO4PFM
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
#endif
   d_calphad_fenergy = std::unique_ptr<CALPHADFreeEnergyFunctionsTernary>(
       new CALPHADFreeEnergyFunctionsTernary(
#ifdef HAVE_THERMO4PFM
           calphad_pt, newton_pt,
#else
           calphad_db, newton_db,
#endif
           energy_interp_func_type, ConcInterpolationType::LINEAR));
}

template <class FreeEnergyType>
int CALPHADequilibriumPhaseConcentrationsStrategy<FreeEnergyType>::
    computePhaseConcentrationsOnPatch(
        std::shared_ptr<pdat::CellData<double> > cd_te,
        std::shared_ptr<pdat::CellData<double> > cd_pf,
        std::shared_ptr<pdat::CellData<double> > cd_eta,
        std::shared_ptr<pdat::CellData<double> > cd_conc,
        std::shared_ptr<pdat::CellData<double> > cd_cl,
        std::shared_ptr<pdat::CellData<double> > cd_ca,
        std::shared_ptr<pdat::CellData<double> > cd_cb,
        std::shared_ptr<hier::Patch> patch)
{
   assert(cd_te);
   assert(cd_pf);
   assert(!cd_eta);  // not supported
   assert(cd_conc);
   assert(cd_cl);
   assert(cd_ca);
   assert(d_calphad_fenergy != nullptr);
   assert(cd_conc->getDepth() == cd_cl->getDepth());
   assert(cd_conc->getDepth() == cd_ca->getDepth());
   assert(cd_cl->getGhostCellWidth()[0] <= cd_te->getGhostCellWidth()[0]);
   assert(cd_cl->getGhostCellWidth()[0] <= cd_pf->getGhostCellWidth()[0]);
#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> cops;
   double l2n = cops.L2Norm(cd_conc, patch->getBox());
   assert(l2n == l2n);
#endif
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

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
   if (d_with_third_phase || nphases == 3) {
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
   if (d_with_third_phase || nphases == 3) {
      N += cd_conc->getDepth();
   }

   double* cl = cd_cl->getPointer(0);
   double* ca = cd_ca->getPointer(0);
   double* cb = nullptr;
   double* conc = cd_conc->getPointer(0);
   double* cl_ref = cd_cl_ref->getPointer(0);
   double* ca_ref = cd_ca_ref->getPointer(0);
   double* cb_ref = nullptr;
   if (d_with_third_phase || nphases == 3) {
      cb = cd_cb->getPointer(0);
      cb_ref = cd_cb_ref->getPointer(0);
   }

   // number of cells for each field
   const int ncp = pf_gbox.size();
   const int ncc = ci_gbox.size();
   const int nct = temp_gbox.size();

   int nits = 0;

#ifdef GPU_OFFLOAD
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
#ifdef HAVE_THERMO4PFM
               for (short i = 0; i < nphases; i++)
                  hphi[i] =
                      Thermo4PFM::interp_func(d_conc_interp_func_type, phi[i]);
#else
            const char conc_inter = concInterpChar(d_conc_interp_func_type);
            for (short i = 0; i < nphases; i++)
               hphi[i] = INTERP_FUNC(phi[i], &conc_inter);
#endif
               double c[2];  // up to 2 components
               for (int ic = 0; ic < nc; ic++) {
                  c[ic] = conc[ic * ncp + idx_pf];
               }
               double x[4];
               for (int ic = 0; ic < nc; ic++) {
                  x[ic] = cl_ref[ic * ncc + idx_ci];
                  x[ic + nc] = ca_ref[ic * ncc + idx_ci];
                  if (d_with_third_phase || nphases == 3) {
                     x[ic + 2 * nc] = cb_ref[ic * ncc + idx_ci];
                  }
               }
               assert(!std::isnan(x[0]));
               assert(!std::isnan(x[nc]));

               // compute cL, cS
               int status =
                   d_calphad_fenergy->computePhaseConcentrations(temp, c,
#ifdef HAVE_THERMO4PFM
                                                                 hphi,
#else
                                                              hphi[0], 0.,
#endif
                                                                 x);
               if (status < 0) {
                  std::cerr
                      << "computePhaseConcentrations failed for T=" << temp
                      << ", hphi=";
                  for (short i = 0; i < nphases; i++)
                     std::cerr << hphi[i] << ", ";
                  std::cerr << ", c=" << c[0] << std::endl;
                  std::cerr << "c_ref=" << cl_ref[0] << "," << ca_ref[0] << ","
                            << cb_ref[0] << std::endl;
                  std::cerr << "x=" << x[0] << "," << x[1] << "," << x[2]
                            << std::endl;
                  MPI_Abort(mpi.getCommunicator(), -1);
               }
#ifndef GPU_OFFLOAD
               assert(!std::isnan(x[0]));
               /*
                              if(std::isnan(x[0]))
                              {
                                 std::cerr
                                     << "computePhaseConcentrations failed for
                  T=" << temp
                                     << ", hphi=";
                                 for (short i = 0; i < nphases; i++)
                                    std::cerr << hphi[i] << ", ";
                                 std::cerr << "c=" << c[0] << std::endl;
                                 std::cerr << ", c_ref=" << cl_ref[0] << "," <<
                  ca_ref[0] << ","
                                           << cb_ref[0] << std::endl;
                                 std::cerr << ", x=" << x[0] << "," << x[1] <<
                  "," << x[2]
                                           << ", idx_pf="<<idx_pf<<",
                  imin[0]="<<imin[0]
                                           << std::endl;
                                 abort();
                              }
               */
               nits += status;
#endif
               // std::cout << "phi=" << phi[0] << "," << phi[1] << "," <<
               // phi[2]
               //          << "c=" << c[0] << ", x=" << x[0] << "," << x[1] <<
               //          ","
               //          << x[2] << std::endl;

               // set cell values with cL and cS just computed
               for (int ic = 0; ic < nc; ic++) {
                  cl[ic * ncc + idx_ci] = x[ic];
                  ca[ic * ncc + idx_ci] = x[ic + nc];
                  if (d_with_third_phase || nphases == 3) {
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
