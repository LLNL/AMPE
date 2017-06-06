// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
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
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
#include "CALPHADequilibriumPhaseConcentrationsStrategy.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"

using namespace std;

CALPHADequilibriumPhaseConcentrationsStrategy::CALPHADequilibriumPhaseConcentrationsStrategy(
      const int conc_l_id,
      const int conc_a_id,
      const int conc_b_id,
      const int conc_l_ref_id,
      const int conc_a_ref_id,
      const int conc_b_ref_id,
      const int ncompositions,
      const string& phase_interp_func_type,
      const string& eta_interp_func_type,
      const string& avg_func_type,
      const bool with_third_phase,
      const double  phase_well_scale,
      const double eta_well_scale,
      const std::string& phase_well_func_type,
      const std::string& eta_well_func_type,
      boost::shared_ptr<tbox::Database> calphad_db,
      boost::shared_ptr<tbox::Database> newton_db):
   PhaseConcentrationsStrategy(
      conc_l_id,
      conc_a_id,
      conc_b_id,
      with_third_phase),
   d_conc_l_ref_id(conc_l_ref_id),
   d_conc_a_ref_id(conc_a_ref_id),
   d_conc_b_ref_id(conc_b_ref_id),
   d_ncompositions(ncompositions)
{
   if( d_ncompositions>1 )
   d_calphad_fenergy = new
      CALPHADFreeEnergyFunctionsTernary(calphad_db,newton_db,
                                 phase_interp_func_type,
                                 avg_func_type,
                                 phase_well_scale,
                                 phase_well_func_type);  
   else
   d_calphad_fenergy = new
      CALPHADFreeEnergyFunctionsBinary(calphad_db,newton_db,
                                 phase_interp_func_type,
                                 eta_interp_func_type,avg_func_type,
                                 with_third_phase,
                                 phase_well_scale,eta_well_scale,
                                 phase_well_func_type,eta_well_func_type);
}

void CALPHADequilibriumPhaseConcentrationsStrategy::computePhaseConcentrationsOnPatch(
   boost::shared_ptr< pdat::CellData<double> > cd_temperature,
   boost::shared_ptr< pdat::CellData<double> > cd_phi,
   boost::shared_ptr< pdat::CellData<double> > cd_eta,
   boost::shared_ptr< pdat::CellData<double> > cd_concentration,
   boost::shared_ptr< pdat::CellData<double> > cd_c_l,
   boost::shared_ptr< pdat::CellData<double> > cd_c_a,
   boost::shared_ptr< pdat::CellData<double> > cd_c_b,
   boost::shared_ptr<hier::Patch > patch )
{
   assert( cd_temperature );
   assert( cd_phi );
   assert( cd_concentration );
   assert( cd_c_l );
   assert( cd_c_a );
   assert( d_calphad_fenergy!=NULL );
   assert( d_ncompositions==cd_concentration->getDepth() );

   const hier::Box& pbox = patch->getBox();

   boost::shared_ptr< pdat::CellData<double> > cd_c_l_ref (
      patch->getPatchData( d_conc_l_ref_id ), boost::detail::dynamic_cast_tag());
   assert( cd_c_l_ref );
   
   boost::shared_ptr< pdat::CellData<double> > cd_c_a_ref (
      patch->getPatchData( d_conc_a_ref_id ), boost::detail::dynamic_cast_tag());
   assert( cd_c_a_ref );
   
   boost::shared_ptr< pdat::CellData<double> > cd_c_b_ref;
   if ( d_with_third_phase ) {
      cd_c_b_ref = boost::dynamic_pointer_cast<pdat::CellData<double>,
                                        hier::PatchData>( patch->getPatchData( d_conc_b_ref_id )); 
      assert( cd_c_b_ref );
   }

   
   const double* const ptr_temp = cd_temperature->getPointer();
   const double* const ptr_phi = cd_phi->getPointer();

   double* ptr_eta = NULL;
   if ( d_with_third_phase ) {
      ptr_eta = cd_eta->getPointer();
   }

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
   assert( cd_phi->getGhostCellWidth()[0]==cd_concentration->getGhostCellWidth()[0] );
   
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

   // Assuming c_l, c_a, and c_b all have same box
   assert( cd_c_l->getGhostCellWidth()[0]==cd_c_a->getGhostCellWidth()[0] );
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
         
   int N = 2*d_ncompositions;
   if ( d_with_third_phase ) {
      N+=d_ncompositions;
   }

   double* x = new double[N];
   double* fA = new double[N];
   double* fB = new double[N];
   double* L0 = new double[N];
   double* L1 = new double[N];
   double* L2 = new double[N];

   for ( int ii = 0; ii < N; ii++ ) {
      x[ii] = 0.0;
      fA[ii] = 0.0;
      fB[ii] = 0.0;
      L0[ii] = 0.0;
      L1[ii] = 0.0;
      L2[ii] = 0.0;
   }

   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;

            const int idx_c_i = (ii - imin_c_i) +
               (jj - jmin_c_i) * jp_c_i + (kk - kmin_c_i) * kp_c_i;

            const double t = ptr_temp[idx_temp];
            
            const double phi = ptr_phi[idx_pf];
            double eta = 0.0;

            if ( d_with_third_phase ) {
               eta = ptr_eta[idx_pf];
            }
            
            vector<double> c(d_ncompositions);
            
            // loop over atomic species
            for(int ic=0;ic<cd_concentration->getDepth();ic++)
            {
               const double* const ptr_conc = cd_concentration->getPointer(ic);
               double* ptr_c_l_ref = cd_c_l_ref->getPointer(ic);
               double* ptr_c_a_ref = cd_c_a_ref->getPointer(ic);
               double* ptr_c_b_ref = NULL;
               if ( d_with_third_phase ) {
                  ptr_c_b_ref = cd_c_b_ref->getPointer(ic);
               }

               c[ic] = ptr_conc[idx_pf];
               x[ic] = ptr_c_l_ref[idx_c_i];
               x[ic+d_ncompositions] = ptr_c_a_ref[idx_c_i];

               if ( d_with_third_phase ) {
                  x[ic+d_ncompositions*2] = ptr_c_b_ref[idx_c_i];

               }
            }

            d_calphad_fenergy->computePhaseConcentrations(
                  t,&c[0],phi,eta,x);

            for(int ic=0;ic<cd_concentration->getDepth();ic++)
            {
               double* ptr_c_l = cd_c_l->getPointer(ic);
               double* ptr_c_a = cd_c_a->getPointer(ic);
               double* ptr_c_b = NULL;
               if ( d_with_third_phase ) {
                  ptr_c_b = cd_c_b->getPointer(ic);
               }
               ptr_c_l[idx_c_i] = x[ic];
               ptr_c_a[idx_c_i] = x[ic+d_ncompositions];
               if ( d_with_third_phase ) {
                  ptr_c_b[idx_c_i] = x[ic+d_ncompositions*2];
               }
            } // ic
         }
      }
   }

   delete[] x;
   delete[] fA;
   delete[] fB;
   delete[] L0;
   delete[] L1;
   delete[] L2;
}
