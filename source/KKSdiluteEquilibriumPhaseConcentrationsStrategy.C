// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "KKSdiluteEquilibriumPhaseConcentrationsStrategy.h"

#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/tbox/IEEE.h"

using namespace std;

KKSdiluteEquilibriumPhaseConcentrationsStrategy::
    KKSdiluteEquilibriumPhaseConcentrationsStrategy(
        const int conc_l_scratch_id, const int conc_a_scratch_id,
        const int conc_b_scratch_id, const int conc_l_ref_id,
        const int conc_a_ref_id, const int conc_b_ref_id,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        boost::shared_ptr<tbox::Database> conc_db)
    : PhaseConcentrationsStrategy(conc_l_scratch_id, conc_a_scratch_id,
                                  conc_b_scratch_id, false),
      d_conc_l_ref_id(conc_l_ref_id),
      d_conc_a_ref_id(conc_a_ref_id),
      d_conc_b_ref_id(conc_b_ref_id)
{
   d_fenergy =
       new KKSFreeEnergyFunctionDiluteBinary(conc_db, energy_interp_func_type,
                                             conc_interp_func_type);
}

void KKSdiluteEquilibriumPhaseConcentrationsStrategy::
    computePhaseConcentrationsOnPatch(
        boost::shared_ptr<pdat::CellData<double> > cd_te,
        boost::shared_ptr<pdat::CellData<double> > cd_pf,
        boost::shared_ptr<pdat::CellData<double> > cd_eta,
        boost::shared_ptr<pdat::CellData<double> > cd_conc,
        boost::shared_ptr<pdat::CellData<double> > cd_cl,
        boost::shared_ptr<pdat::CellData<double> > cd_ca,
        boost::shared_ptr<pdat::CellData<double> > cd_cb,
        boost::shared_ptr<hier::Patch> patch)
{
   (void)cd_eta;

   assert(cd_te);
   assert(cd_pf);
   assert(cd_conc);
   assert(cd_cl);
   assert(cd_ca);
   assert(d_fenergy != NULL);
   assert(cd_conc->getDepth() == cd_cl->getDepth());
   assert(cd_conc->getDepth() == cd_ca->getDepth());
   assert(cd_cl->getGhostCellWidth()[0] <= cd_te->getGhostCellWidth()[0]);
   assert(cd_cl->getGhostCellWidth()[0] <= cd_pf->getGhostCellWidth()[0]);
#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> cops;
   double l2n = cops.L2Norm(cd_conc, patch->getBox());
   assert(l2n == l2n);
#endif

   //   assert( d_conc_l_ref_id>=0 );
   //   boost::shared_ptr< pdat::CellData<double> > cd_cl_ref (
   //      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
   //         patch->getPatchData( d_conc_l_ref_id) ) );
   //   assert( cd_cl_ref );
   //   assert( cd_conc->getDepth()==cd_cl_ref->getDepth() );

   //   assert( d_conc_a_ref_id>=0 );
   //   boost::shared_ptr< pdat::CellData<double> > cd_ca_ref (
   //      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
   //         patch->getPatchData( d_conc_a_ref_id) ) );
   //   assert( cd_ca_ref );
   //   assert( cd_conc->getDepth()==cd_ca_ref->getDepth() );

   const double* const ptr_temp = cd_te->getPointer();
   const double* const ptr_phi = cd_pf->getPointer();

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

   // Assuming phi, eta, and concentration all have same box
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
   // assert( cd_cl->getGhostCellWidth()[0]==cd_cl_ref->getGhostCellWidth()[0]
   // );
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


   double x[2];
   for (short i = 0; i < 2; i++)
      x[i] = tbox::IEEE::getSignalingNaN();

   int idx_pf = (imin[0] - imin_pf) + (imin[1] - jmin_pf) * inc_j_pf +
                (imin[2] - kmin_pf) * inc_k_pf;
   int idx_te = (imin[0] - imin_te) + (imin[1] - jmin_te) * inc_j_te +
                (imin[2] - kmin_te) * inc_k_te;
   int idx_ci = (imin[0] - imin_ci) + (imin[1] - jmin_ci) * inc_j_ci +
                (imin[2] - kmin_ci) * inc_k_ci;

   for (int kk = imin[2]; kk <= imax[2]; kk++) {
      for (int jj = imin[1]; jj <= imax[1]; jj++) {
         for (int ii = imin[0]; ii <= imax[0]; ii++) {

            const double temp = ptr_temp[idx_te];
            const double phi = ptr_phi[idx_pf];

            double c = cd_conc->getPointer(0)[idx_pf];
            TBOX_ASSERT(c == c);

            //            x[0] = cd_cl_ref->getPointer(0)[idx_ci];
            //            x[1] = cd_ca_ref->getPointer(0)[idx_ci];
            x[0] = cd_cl->getPointer(0)[idx_ci];
            x[1] = cd_ca->getPointer(0)[idx_ci];

            // compute cL, cS
            d_fenergy->computePhaseConcentrations(temp, &c, phi, 0., x);
            //            std::cout<<"cl, cs = "<<x[0]<<","<<x[1]<<std::endl;

            // set cell values with cL and cS just computed
            cd_cl->getPointer(0)[idx_ci] = x[0];
            cd_ca->getPointer(0)[idx_ci] = x[1];

            idx_pf++;
            idx_te++;
            idx_ci++;

         }  // ii

         idx_pf += 2 * (cd_pf->getGhostCellWidth()[0] -
                        cd_cl->getGhostCellWidth()[0]);
         idx_te += 2 * (cd_te->getGhostCellWidth()[0] -
                        cd_cl->getGhostCellWidth()[0]);
      }  // jj

      idx_pf += 2 * inc_j_pf *
                (cd_pf->getGhostCellWidth()[1] - cd_cl->getGhostCellWidth()[1]);
      idx_te += 2 * inc_j_te *
                (cd_te->getGhostCellWidth()[1] - cd_cl->getGhostCellWidth()[1]);
   }  // kk
}
