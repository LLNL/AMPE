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
#include "HBSMequilibriumPhaseConcentrationsStrategy.h"
#include "HBSMFreeEnergyStrategy.h"
#include "FuncFort.h"

using namespace std;

HBSMequilibriumPhaseConcentrationsStrategy::HBSMequilibriumPhaseConcentrationsStrategy(
   const int conc_l_id,
   const int conc_a_id,
   const int conc_b_id,
   const QuatModelParameters& model_parameters,
   boost::shared_ptr<tbox::Database> conc_db):
      d_phase_interp_func_type(model_parameters.energy_interp_func_type()),
      d_eta_interp_func_type(model_parameters.eta_interp_func_type()),
      PhaseConcentrationsStrategy(
         conc_l_id,
         conc_a_id,
         conc_b_id,
         model_parameters.with_third_phase())
{
   d_hbsm_fenergy = new HBSMFreeEnergyStrategy(
               conc_db->getDatabase( "HBSM" ),
               model_parameters.energy_interp_func_type(),
               model_parameters.eta_interp_func_type(),
               model_parameters.molar_volume_liquid(),
               model_parameters.molar_volume_solid_A(),
               model_parameters.molar_volume_solid_B(),
               model_parameters.D_liquid(),
               model_parameters.D_solid_A(),
               model_parameters.D_solid_B(),
               model_parameters.Q0_liquid(),
               model_parameters.Q0_solid_A(),
               model_parameters.Q0_solid_B(),
               conc_l_id,
               conc_a_id,
               conc_b_id,
               model_parameters.with_third_phase() );
}

void HBSMequilibriumPhaseConcentrationsStrategy::computePhaseConcentrationsOnPatch(
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
   
   const hier::Box& pbox = patch->getBox();

   double* ptr_phi = cd_phi->getPointer();
   double* ptr_conc = cd_concentration->getPointer();
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();
   double* ptr_c_b = NULL;
   
   if ( d_with_third_phase ) {
      ptr_c_b = cd_c_b->getPointer();
   }
   double* ptr_eta = NULL;
   if ( d_with_third_phase ) {
      ptr_eta = cd_eta->getPointer();
   }

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
         
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;

            const int idx_c_i = (ii - imin_c_i) +
               (jj - jmin_c_i) * jp_c_i + (kk - kmin_c_i) * kp_c_i;

            double phi = ptr_phi[idx_pf];
            double c = ptr_conc[idx_pf];
            double eta = 0.0;
            if ( d_with_third_phase ) {
               eta = ptr_eta[idx_pf];
            }

            double hphi =
               FORT_INTERP_FUNC(
                  phi,
                  d_phase_interp_func_type.c_str() );

            double heta = 0.0;
            if ( d_with_third_phase ) {
               heta =
                  FORT_INTERP_FUNC(
                     eta,
                     d_eta_interp_func_type.c_str() );
            }

            double c_l = d_hbsm_fenergy->computeLiquidConcentration(
               hphi, heta, c );

            double c_a = d_hbsm_fenergy->computeSolidAConcentration(
               hphi, heta, c );

            ptr_c_l[idx_c_i] = c_l;
            ptr_c_a[idx_c_i] = c_a;
            if ( d_with_third_phase ) {
               double c_b = d_hbsm_fenergy->computeSolidBConcentration(
                  hphi, heta, c );
               ptr_c_b[idx_c_i] = c_b;
            }

         }
      }
   }
}
