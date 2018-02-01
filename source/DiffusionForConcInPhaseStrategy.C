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
#include "DiffusionForConcInPhaseStrategy.h"
#include "CompositionStrategyMobilities.h"
#include "FreeEnergyStrategy.h"

#include <vector>
using namespace std;

void small_mat_mult(const int n, const double* const a, const double* const b, double* c)
{
   for(short k=0;k<n;k++)
   for(short i=0;i<n;i++){
      c[n*k+i]=0.;
      for(short j=0;j<n;j++)
         c[n*k+i]+=a[n*k+j]*b[n*j+k];
   }
}

DiffusionForConcInPhaseStrategy::DiffusionForConcInPhaseStrategy(
   const unsigned short ncompositions,
   const int conc_l_scratch_id,
   const int conc_a_scratch_id,
   const int conc_b_scratch_id,
   const int diffusion_l_id,
   const int diffusion_a_id,
   const int diffusion_b_id,
   const int diffusion_coeff_l_id,
   const int diffusion_coeff_a_id,
   const int diffusion_coeff_b_id,
   const std::string& avg_func_type,
   const std::string& diff_interp_func_type,
   CompositionStrategyMobilities* mobilities_strategy,
   FreeEnergyStrategy* free_energy_strategy):
      d_mobilities_strategy(mobilities_strategy),
      d_free_energy_strategy(free_energy_strategy)
{
   d_ncompositions=ncompositions;
   
   d_conc_l_scratch_id=conc_l_scratch_id;
   d_conc_a_scratch_id=conc_a_scratch_id;
   d_conc_b_scratch_id=conc_b_scratch_id;
   
   d_diffusion_l_id=diffusion_l_id;
   d_diffusion_a_id=diffusion_a_id;
   d_diffusion_b_id=diffusion_b_id;
   
   d_diffusion_coeff_l_id=diffusion_coeff_l_id;
   d_diffusion_coeff_a_id=diffusion_coeff_a_id;
   d_diffusion_coeff_b_id=diffusion_coeff_b_id;

   d_avg_func_type=avg_func_type;
   d_diff_interp_func_type=diff_interp_func_type;

   d_with_third_phase = ( conc_b_scratch_id>=0 );

   d_same_composition_for_third_phase=false;
   
   assert( d_free_energy_strategy!=NULL );

   d_d2f.resize(       d_ncompositions*d_ncompositions);
   d_mobmat.resize(    d_ncompositions*d_ncompositions);
   d_local_dmat.resize(d_ncompositions*d_ncompositions);
}

void DiffusionForConcInPhaseStrategy::computeLocalDiffusionMatrixL(
   const double temperature,
   const vector<double>& c)
{
   d_free_energy_strategy->computeSecondDerivativeEnergyPhaseL(temperature, c, d_d2f, false);
   d_mobilities_strategy->computeDiffusionMobilityPhaseL(c, temperature, d_mobmat);

   small_mat_mult(d_ncompositions,&d_mobmat[0],&d_d2f[0],&d_local_dmat[0]);

   for(short i=0;i<d_ncompositions*d_ncompositions;i++)assert( d_local_dmat[i]==d_local_dmat[i] );
}

void DiffusionForConcInPhaseStrategy::computeLocalDiffusionMatrixA(
   const double temperature,
   const vector<double>& c)
{
   d_free_energy_strategy->computeSecondDerivativeEnergyPhaseA(temperature, c, d_d2f, false);
   d_mobilities_strategy->computeDiffusionMobilityPhaseA(c, temperature, d_mobmat);

   small_mat_mult(d_ncompositions,&d_mobmat[0],&d_d2f[0],&d_local_dmat[0]);

   for(short i=0;i<d_ncompositions*d_ncompositions;i++)assert( d_local_dmat[i]==d_local_dmat[i] );
}

void DiffusionForConcInPhaseStrategy::computeLocalDiffusionMatrixB(
   const double temperature,
   const vector<double>& c)
{
   d_free_energy_strategy->computeSecondDerivativeEnergyPhaseB(temperature, c, d_d2f, false);
   d_mobilities_strategy->computeDiffusionMobilityPhaseB(c, temperature, d_mobmat);

   small_mat_mult(d_ncompositions,&d_mobmat[0],&d_d2f[0],&d_local_dmat[0]);
}

void DiffusionForConcInPhaseStrategy::setDiffusionForConcInPhase(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int phase_id,
   const int eta_id)
{
   //tbox::pout<<"DiffusionForConcInPhaseStrategy::setDiffusionForConcInPhase"<<endl;
   assert( phase_id >= 0 );
   assert( d_diffusion_l_id >= 0 );
   assert( d_diffusion_a_id >= 0 );
   if ( d_with_third_phase ) {
      assert( eta_id >= 0 );
      assert( d_diffusion_b_id >= 0 );
   }
   assert( d_diffusion_coeff_l_id >= 0 );
   assert( d_diffusion_coeff_a_id >= 0 );
   if ( d_with_third_phase ) {
      assert( d_diffusion_coeff_b_id >= 0 );
   }

   const int maxl = hierarchy->getNumberOfLevels();

   for ( int amr_level = 0; amr_level < maxl; amr_level++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( amr_level );

      for ( hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); p++ ) {
         boost::shared_ptr<hier::Patch > patch = *p;
      
         boost::shared_ptr< pdat::CellData<double> > phi (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( phase_id) ) );
         
         boost::shared_ptr< pdat::SideData<double> > diffusion_l (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( d_diffusion_l_id) ) );
         boost::shared_ptr< pdat::SideData<double> > diffusion_a (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( d_diffusion_a_id) ) );
         boost::shared_ptr< pdat::SideData<double> > diffusion_b;
         
         boost::shared_ptr< pdat::SideData<double> > diffusion_coeff_l (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( d_diffusion_coeff_l_id) ) );
         boost::shared_ptr< pdat::SideData<double> > diffusion_coeff_a (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( d_diffusion_coeff_a_id) ) );
         boost::shared_ptr< pdat::SideData<double> > diffusion_coeff_b;

         boost::shared_ptr< pdat::CellData<double> > eta;
         if ( d_with_third_phase ) {
            eta = boost::dynamic_pointer_cast<pdat::CellData<double>,hier::PatchData> ( patch->getPatchData( eta_id ));
            diffusion_b = boost::dynamic_pointer_cast<pdat::SideData<double>,hier::PatchData>( patch->getPatchData( d_diffusion_b_id ));
            diffusion_coeff_b = boost::dynamic_pointer_cast<pdat::SideData<double>,hier::PatchData>( patch->getPatchData( d_diffusion_coeff_b_id ));
         }

         setDiffusionForConcInPhaseOnPatch(
            phi, eta,
            diffusion_coeff_l, diffusion_coeff_a, diffusion_coeff_b,
            diffusion_l, diffusion_a, diffusion_b, 
            patch->getBox() );
      }
   }
}

//-----------------------------------------------------------------------

void DiffusionForConcInPhaseStrategy::setDiffusionCoeffForConcInPhase(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int eta_scratch_id)
{
   //tbox::pout<<"DiffusionForConcInPhaseStrategy::setDiffusionCoeffForConcInPhase"<<endl;
   assert( temperature_id >= 0 );
   assert( d_diffusion_coeff_l_id >= 0 );
   assert( d_diffusion_coeff_a_id >= 0 );
   if ( d_with_third_phase ) {
      assert( d_diffusion_coeff_b_id >= 0 );
   }

   const int maxl = hierarchy->getNumberOfLevels();

   for ( int amr_level = 0; amr_level < maxl; amr_level++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( amr_level );

      for ( hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); p++ ) {
         boost::shared_ptr<hier::Patch > patch = *p;

         boost::shared_ptr< pdat::CellData<double> > temp (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( temperature_id) ) );

         boost::shared_ptr< pdat::CellData<double> > cl (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_conc_l_scratch_id) ) );
         assert( cl );
         boost::shared_ptr< pdat::CellData<double> > ca (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_conc_a_scratch_id) ) );
         assert( ca );
         boost::shared_ptr< pdat::CellData<double> > cb;

         boost::shared_ptr< pdat::SideData<double> > diffusion_coeff_l (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( d_diffusion_coeff_l_id) ) );
         boost::shared_ptr< pdat::SideData<double> > diffusion_coeff_a (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( d_diffusion_coeff_a_id) ) );
         boost::shared_ptr< pdat::SideData<double> > diffusion_coeff_b;

         boost::shared_ptr< pdat::CellData<double> > eta;
         if ( d_with_third_phase ) {
            cb  = boost::dynamic_pointer_cast<pdat::CellData<double>,hier::PatchData>( patch->getPatchData( d_conc_b_scratch_id ));
            assert( cb );
            diffusion_coeff_b = boost::dynamic_pointer_cast<pdat::SideData<double>,hier::PatchData>( patch->getPatchData( d_diffusion_coeff_b_id ));
            eta = boost::dynamic_pointer_cast<pdat::CellData<double>,hier::PatchData>( patch->getPatchData( eta_scratch_id ));
         }

         setDiffusionCoeffForConcInPhaseOnPatch(
            cl, ca, cb, temp, eta,
            diffusion_coeff_l, diffusion_coeff_a, diffusion_coeff_b,
            patch->getBox() );
      }
   }
}

//-----------------------------------------------------------------------

void DiffusionForConcInPhaseStrategy::setDiffusionCoeffForConcInPhaseOnPatch(
   boost::shared_ptr< pdat::CellData<double> > cd_c_l,
   boost::shared_ptr< pdat::CellData<double> > cd_c_a,
   boost::shared_ptr< pdat::CellData<double> > cd_c_b,
   boost::shared_ptr< pdat::CellData<double> > cd_temp,
   boost::shared_ptr< pdat::CellData<double> > cd_eta,
   boost::shared_ptr< pdat::SideData<double> > sd_d_coeff_l,
   boost::shared_ptr< pdat::SideData<double> > sd_d_coeff_a,
   boost::shared_ptr< pdat::SideData<double> > sd_d_coeff_b,
   const hier::Box& pbox )
{
   assert( cd_c_l );
   assert( cd_c_a );
   
   vector<double*> ptr_c_l(d_ncompositions);
   vector<double*> ptr_c_a(d_ncompositions);
   vector<double*> ptr_c_b(d_ncompositions);
   for(unsigned short ic=0;ic<d_ncompositions;ic++){
      ptr_c_l[ic] = cd_c_l->getPointer(ic);
      ptr_c_a[ic] = cd_c_a->getPointer(ic);
   }
   
   double* ptr_eta = NULL;
   if ( d_with_third_phase ) {
      for(unsigned short ic=0;ic<d_ncompositions;ic++)
         ptr_c_b[ic] = cd_c_b->getPointer(ic);
      ptr_eta = cd_eta->getPointer();
   }

   double* ptr_temp = cd_temp->getPointer();

   const int nc2=d_ncompositions*d_ncompositions;
   vector<double*> ptr_dx_coeff_l(nc2);
   vector<double*> ptr_dy_coeff_l(nc2);
   vector<double*> ptr_dz_coeff_l(nc2, NULL);
   for(unsigned short ic=0;ic<d_ncompositions;ic++)
   for(unsigned short jc=0;jc<d_ncompositions;jc++){
      const unsigned ijc=ic+jc*d_ncompositions;
      ptr_dx_coeff_l[ijc] = sd_d_coeff_l->getPointer( 0, ijc );
      ptr_dy_coeff_l[ijc] = sd_d_coeff_l->getPointer( 1, ijc );
      if ( NDIM > 2 ) {
         ptr_dz_coeff_l[ijc] = sd_d_coeff_l->getPointer( 2, ijc );
      }
   }
   
   vector<double*> ptr_dx_coeff_a(nc2);
   vector<double*> ptr_dy_coeff_a(nc2);
   vector<double*> ptr_dz_coeff_a(nc2,NULL);
   for(unsigned short ic=0;ic<d_ncompositions;ic++)
   for(unsigned short jc=0;jc<d_ncompositions;jc++){
      const unsigned ijc=ic+jc*d_ncompositions;
      ptr_dx_coeff_a[ijc] = sd_d_coeff_a->getPointer( 0, ijc );
      ptr_dy_coeff_a[ijc] = sd_d_coeff_a->getPointer( 1, ijc );
      if ( NDIM > 2 ) {
         ptr_dz_coeff_a[ijc] = sd_d_coeff_a->getPointer( 2, ijc );
      }
   }
   
   vector<double*> ptr_dx_coeff_b;
   vector<double*> ptr_dy_coeff_b;
   vector<double*> ptr_dz_coeff_b;
   if ( d_with_third_phase ) {
      ptr_dx_coeff_b.resize(nc2, NULL);
      ptr_dy_coeff_b.resize(nc2, NULL);
      ptr_dz_coeff_b.resize(nc2, NULL);

      for(unsigned short ic=0;ic<d_ncompositions;ic++)
      for(unsigned short jc=0;jc<d_ncompositions;jc++){
         const unsigned ijc=ic+jc*d_ncompositions;
         ptr_dx_coeff_b[ijc] = sd_d_coeff_b->getPointer( 0, ijc );
         ptr_dy_coeff_b[ijc] = sd_d_coeff_b->getPointer( 1, ijc );
         if ( NDIM > 2 ) {
            ptr_dz_coeff_b[ijc] = sd_d_coeff_b->getPointer( 2, ijc );
         }
      }
   }
   
   // Assuming all sd_d_* have same box
   const hier::Box& dcoeff_gbox = sd_d_coeff_l->getGhostBox();
   int imin_dcoeff = dcoeff_gbox.lower(0);
   int jmin_dcoeff = dcoeff_gbox.lower(1);
   int jp_dcoeff = dcoeff_gbox.numberCells(0);
   int kmin_dcoeff = 0;
   int kp_dcoeff = 0;
#if (NDIM == 3)
   kmin_dcoeff = dcoeff_gbox.lower(2);
   kp_dcoeff = jp_dcoeff * dcoeff_gbox.numberCells(1);
#endif

   // assumes all c have same number of ghosts
   const hier::Box& c_gbox = cd_c_l->getGhostBox();
   int imin_c = c_gbox.lower(0);
   int jmin_c = c_gbox.lower(1);
   int jp_c   = c_gbox.numberCells(0);
   int kmin_c = 0;
   int kp_c   = 0;
#if (NDIM == 3)
   kmin_c = c_gbox.lower(2);
   kp_c   = jp_c * c_gbox.numberCells(1);
#endif

   const hier::Box& temp_gbox = cd_temp->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp   = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp   = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp   = jp_temp * temp_gbox.numberCells(1);
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
         
   vector<double> c_l(d_ncompositions);
   vector<double> c_a(d_ncompositions);
   vector<double> c_b(d_ncompositions);
   vector<double> c_s(d_ncompositions);
   
   // X-side
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax+1; ii++ ) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
               (jj - jmin_dcoeff) * (jp_dcoeff + 1) +
               (kk - kmin_dcoeff) * (kp_dcoeff + dcoeff_gbox.numberCells(1) );

            const int idx_c = (ii - imin_c) +
               (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
            const int idxm1_c = idx_c - 1;

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - 1;
            
            const int idx_pf  =idx_c;
            const int idxm1_pf=idxm1_c;

            double temp = 0.5 * ( ptr_temp[idx_temp] + ptr_temp[idxm1_temp] );

            for(unsigned short ic=0;ic<d_ncompositions;ic++){
               assert( ptr_c_l[ic][idx_c]>=0. );
               assert( ptr_c_l[ic][idxm1_c]>=0. );
               c_l[ic] = 0.5 * ( ptr_c_l[ic][idx_c] + ptr_c_l[ic][idxm1_c] );
               c_a[ic] = 0.5 * ( ptr_c_a[ic][idx_c] + ptr_c_a[ic][idxm1_c] );
            }
           
            computeLocalDiffusionMatrixL(temp, c_l); 
            for(int ic=0;ic<nc2;ic++)
               ptr_dx_coeff_l[ic][idx_dcoeff] = d_local_dmat[ic];

            double eta=0.;
            if ( d_with_third_phase ) {
               eta = 0.5*( ptr_eta[idx_pf] + ptr_eta[idxm1_pf] );
               for(unsigned short ic=0;ic<d_ncompositions;ic++){
                  c_b[ic] = 0.5 * ( ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c] );
               }
            }
            
            if( d_with_third_phase && d_same_composition_for_third_phase )
            {
               double heta =
                  FORT_INTERP_FUNC( eta, d_diff_interp_func_type.c_str() );
               for(unsigned short ic=0;ic<d_ncompositions;ic++){
                  c_s[ic] = (1.-heta)*c_a[ic]+heta*c_b[ic];
               }
               computeLocalDiffusionMatrixA(temp, c_s);
               for(int ic=0;ic<nc2;ic++){
                  ptr_dx_coeff_a[ic][idx_dcoeff] = d_local_dmat[ic];
                  ptr_dx_coeff_b[ic][idx_dcoeff] = d_local_dmat[ic];
               }
            }else{
               computeLocalDiffusionMatrixA(temp, c_a);
               for(int ic=0;ic<nc2;ic++)
                  ptr_dx_coeff_a[ic][idx_dcoeff] = d_local_dmat[ic];
            
               if ( d_with_third_phase ) {
                  computeLocalDiffusionMatrixB(temp, c_b);
                  for(int ic=0;ic<nc2;ic++)
                     ptr_dx_coeff_b[ic][idx_dcoeff] = d_local_dmat[ic];
               }
            }
         }
      }
   }

   // Y-side
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax+1; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
               (jj - jmin_dcoeff) * jp_dcoeff +
               (kk - kmin_dcoeff) * (kp_dcoeff + jp_dcoeff);

            const int idx_c = (ii - imin_c) +
               (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
            const int idxm1_c = idx_c - jp_c;

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - jp_temp;

            for(unsigned short ic=0;ic<d_ncompositions;ic++){
               c_l[ic] = 0.5 * ( ptr_c_l[ic][idx_c] + ptr_c_l[ic][idxm1_c] );
               c_a[ic] = 0.5 * ( ptr_c_a[ic][idx_c] + ptr_c_a[ic][idxm1_c] );
               c_b[ic] = 0.0;
            }

            const int idx_pf  =idx_c;
            const int idxm1_pf=idxm1_c;

            double temp = 0.5 * ( ptr_temp[idx_temp] + ptr_temp[idxm1_temp] );

            if ( d_with_third_phase ) {
               for(unsigned short ic=0;ic<d_ncompositions;ic++)
                  c_b[ic] = 0.5 * ( ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c] );
            }

            computeLocalDiffusionMatrixL(temp,c_l);
            for(int ic=0;ic<nc2;ic++)
               ptr_dy_coeff_l[ic][idx_dcoeff] = d_local_dmat[ic];
            
            double eta=0.;
            if ( d_with_third_phase ) {
               eta = 0.5*( ptr_eta[idx_pf] + ptr_eta[idxm1_pf] );
               for(unsigned short ic=0;ic<d_ncompositions;ic++){
                  c_b[ic] = 0.5 * ( ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c] );
               }
            }
            
            if( d_with_third_phase && d_same_composition_for_third_phase )
            {
               double heta =
                  FORT_INTERP_FUNC( eta, d_diff_interp_func_type.c_str() );
               for(unsigned short ic=0;ic<d_ncompositions;ic++){
                  c_s[ic] = (1.-heta)*c_a[ic]+heta*c_b[ic];
               }
               
               computeLocalDiffusionMatrixA(temp, c_s);
               for(int ic=0;ic<nc2;ic++){
                  ptr_dy_coeff_a[ic][idx_dcoeff] = d_local_dmat[ic];
                  ptr_dy_coeff_b[ic][idx_dcoeff] = d_local_dmat[ic];
               }
            }else{
               computeLocalDiffusionMatrixA(temp, c_a);
               for(int ic=0;ic<nc2;ic++)
                  ptr_dy_coeff_a[ic][idx_dcoeff] = d_local_dmat[ic];
            
               if ( d_with_third_phase ) {
                  computeLocalDiffusionMatrixB(temp, c_b);
                  for(int ic=0;ic<nc2;ic++)
                     ptr_dy_coeff_b[ic][idx_dcoeff] = d_local_dmat[ic];
               }
            }
         }
      }
   }

   if ( NDIM > 2 ) {
      // Z-side
      for ( int kk = kmin; kk <= kmax+1; kk++ ) {
         for ( int jj = jmin; jj <= jmax; jj++ ) {
            for ( int ii = imin; ii <= imax; ii++ ) {

               const int idx_dcoeff = (ii - imin_dcoeff) +
                  (jj - jmin_dcoeff) * jp_dcoeff +
                  (kk - kmin_dcoeff) * kp_dcoeff;

               const int idx_c = (ii - imin_c) +
                  (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
               const int idxm1_c = idx_c - kp_c;

               const int idx_temp = (ii - imin_temp) +
                  (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;
               const int idxm1_temp = idx_temp - kp_temp;

               const int idx_pf  =idx_c;
               const int idxm1_pf=idxm1_c;

               for(unsigned short ic=0;ic<d_ncompositions;ic++){
                  c_l[ic] = 0.5 * ( ptr_c_l[ic][idx_c] + ptr_c_l[ic][idxm1_c] );
                  c_a[ic] = 0.5 * ( ptr_c_a[ic][idx_c] + ptr_c_a[ic][idxm1_c] );
                  c_b[ic] = 0.0;
               }

               double temp = 0.5 * ( ptr_temp[idx_temp] + ptr_temp[idxm1_temp] );

               if ( d_with_third_phase ) {
                  for(unsigned short ic=0;ic<d_ncompositions;ic++)
                     c_b[ic] = 0.5 * ( ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c] );
               }

               computeLocalDiffusionMatrixL(temp, c_l);
               for(int ic=0;ic<nc2;ic++)
                  ptr_dz_coeff_l[ic][idx_dcoeff] = d_local_dmat[ic];
               
               double eta=0.;
               if ( d_with_third_phase ) {
                  eta = 0.5*( ptr_eta[idx_pf] + ptr_eta[idxm1_pf] );
                  for(unsigned short ic=0;ic<d_ncompositions;ic++){
                     c_b[ic] = 0.5 * ( ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c] );
                  }
               }
            
               if( d_with_third_phase && d_same_composition_for_third_phase )
               {
                  double heta =
                     FORT_INTERP_FUNC( eta, d_diff_interp_func_type.c_str() );
                  for(unsigned short ic=0;ic<d_ncompositions;ic++){
                     c_s[ic] = (1.-heta)*c_a[ic]+heta*c_b[ic];
                  }
               
                  computeLocalDiffusionMatrixA(temp, c_s);
                  for(int ic=0;ic<nc2;ic++){
                     ptr_dz_coeff_a[ic][idx_dcoeff] = d_local_dmat[ic];
                     ptr_dz_coeff_b[ic][idx_dcoeff] = d_local_dmat[ic];
                  }
               }else{
               
                  computeLocalDiffusionMatrixA(temp, c_a);
                  for(int ic=0;ic<nc2;ic++)
                     ptr_dz_coeff_a[ic][idx_dcoeff] = d_local_dmat[ic];
               
                  if ( d_with_third_phase ) {
                     computeLocalDiffusionMatrixB(temp, c_b);
                     for(int ic=0;ic<nc2;ic++)
                        ptr_dz_coeff_b[ic][idx_dcoeff] = d_local_dmat[ic];
                  }
               }
            }
         }
      }
   }  // if ( NDIM > 2 )
}

//-----------------------------------------------------------------------

void DiffusionForConcInPhaseStrategy::setDiffusionForConcInPhaseOnPatch(
   boost::shared_ptr< pdat::CellData<double> > cd_phi,
   boost::shared_ptr< pdat::CellData<double> > cd_eta,
   boost::shared_ptr< pdat::SideData<double> > sd_d_coeff_l, 
   boost::shared_ptr< pdat::SideData<double> > sd_d_coeff_a,
   boost::shared_ptr< pdat::SideData<double> > sd_d_coeff_b, 
   boost::shared_ptr< pdat::SideData<double> > sd_d_l, // output
   boost::shared_ptr< pdat::SideData<double> > sd_d_a, // output
   boost::shared_ptr< pdat::SideData<double> > sd_d_b, // output
   const hier::Box& pbox )
{
   //tbox::pout<<"DiffusionForConcInPhaseStrategy::setDiffusionForConcInPhaseOnPatch"<<endl;
   assert( cd_phi );
   assert( sd_d_l );
   assert( sd_d_coeff_l );

   const unsigned nc2=d_ncompositions*d_ncompositions;
   vector<double*> ptr_dx_l;ptr_dx_l.resize(nc2);
   vector<double*> ptr_dy_l;ptr_dy_l.resize(nc2);
   vector<double*> ptr_dz_l;ptr_dz_l.resize(nc2, NULL);
   for(unsigned short ic=0;ic<d_ncompositions;ic++)
   for(unsigned short jc=0;jc<d_ncompositions;jc++){
      const unsigned ijc=ic+jc*d_ncompositions;
      ptr_dx_l[ijc] = sd_d_l->getPointer( 0, ijc );
      ptr_dy_l[ijc] = sd_d_l->getPointer( 1, ijc );
      if ( NDIM > 2 ) {
         ptr_dz_l[ijc] = sd_d_l->getPointer( 2, ijc );
      }
   }
   
   vector<double*> ptr_dx_a(nc2);
   vector<double*> ptr_dy_a(nc2);
   vector<double*> ptr_dz_a(nc2, NULL);
   for(unsigned short ic=0;ic<d_ncompositions;ic++)
   for(unsigned short jc=0;jc<d_ncompositions;jc++){
      const unsigned ijc=ic+jc*d_ncompositions;
      ptr_dx_a[ijc] = sd_d_a->getPointer( 0, ijc );
      ptr_dy_a[ijc] = sd_d_a->getPointer( 1, ijc );
      if ( NDIM > 2 ) {
         ptr_dz_a[ijc] = sd_d_a->getPointer( 2, ijc );
      }
   }
   
   vector<double*> ptr_dx_b(nc2, NULL);
   vector<double*> ptr_dy_b(nc2, NULL);
   vector<double*> ptr_dz_b(nc2, NULL);
   if ( d_with_third_phase ) {
      for(unsigned short ic=0;ic<d_ncompositions;ic++)
      for(unsigned short jc=0;jc<d_ncompositions;jc++){
         const unsigned ijc=ic+jc*d_ncompositions;
         ptr_dx_b[ijc] = sd_d_b->getPointer( 0, ijc );
         ptr_dy_b[ijc] = sd_d_b->getPointer( 1, ijc );
         if ( NDIM > 2 ) {
            ptr_dz_b[ijc] = sd_d_b->getPointer( 2, ijc );
         }
      }
   }
  
   vector<double*> ptr_dx_coeff_l(nc2);
   vector<double*> ptr_dy_coeff_l(nc2);
   vector<double*> ptr_dz_coeff_l(nc2, NULL);
   for(unsigned short ic=0;ic<d_ncompositions;ic++)
   for(unsigned short jc=0;jc<d_ncompositions;jc++){
      const unsigned ijc=ic+jc*d_ncompositions;
      ptr_dx_coeff_l[ijc] = sd_d_coeff_l->getPointer( 0, ijc );
      ptr_dy_coeff_l[ijc] = sd_d_coeff_l->getPointer( 1, ijc );
      if ( NDIM > 2 ) {
         ptr_dz_coeff_l[ijc] = sd_d_coeff_l->getPointer( 2, ijc );
      }
   }

   vector<double*> ptr_dx_coeff_a(nc2);
   vector<double*> ptr_dy_coeff_a(nc2);
   vector<double*> ptr_dz_coeff_a(nc2, NULL);
   for(unsigned short ic=0;ic<d_ncompositions;ic++)
   for(unsigned short jc=0;jc<d_ncompositions;jc++){
      const unsigned ijc=ic+jc*d_ncompositions;
      ptr_dx_coeff_a[ijc] = sd_d_coeff_a->getPointer( 0, ijc );
      ptr_dy_coeff_a[ijc] = sd_d_coeff_a->getPointer( 1, ijc );
      if ( NDIM > 2 ) {
         ptr_dz_coeff_a[ijc] = sd_d_coeff_a->getPointer( 2, ijc );
      }
   }

   vector<double*> ptr_dx_coeff_b(nc2, NULL);
   vector<double*> ptr_dy_coeff_b(nc2, NULL);
   vector<double*> ptr_dz_coeff_b(nc2, NULL);
   if ( d_with_third_phase ) {
      for(unsigned short ic=0;ic<d_ncompositions;ic++)
      for(unsigned short jc=0;jc<d_ncompositions;jc++){
         const unsigned ijc=ic+jc*d_ncompositions;
         ptr_dx_coeff_b[ijc] = sd_d_coeff_b->getPointer( 0, ijc );
         ptr_dy_coeff_b[ijc] = sd_d_coeff_b->getPointer( 1, ijc );
         if ( NDIM > 2 ) {
            ptr_dz_coeff_b[ijc] = sd_d_coeff_b->getPointer( 2, ijc );
         }
      }
   }
 
   double* ptr_phi = cd_phi->getPointer();
   double* ptr_eta = NULL;
   vector<double*> ptr_c_b;ptr_c_b.resize(d_ncompositions);
   if ( d_with_third_phase ) {
      ptr_eta = cd_eta->getPointer();
   }

   // Assuming all sd_d_* have same box
   const hier::Box& dcoeff_gbox = sd_d_l->getGhostBox();
   int imin_dcoeff = dcoeff_gbox.lower(0);
   int jmin_dcoeff = dcoeff_gbox.lower(1);
   int jp_dcoeff = dcoeff_gbox.numberCells(0);
   int kmin_dcoeff = 0;
   int kp_dcoeff = 0;
#if (NDIM == 3)
   kmin_dcoeff = dcoeff_gbox.lower(2);
   kp_dcoeff = jp_dcoeff * dcoeff_gbox.numberCells(1);
#endif

   // Assuming phi and eta have same box
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
         
   // X-side
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax+1; ii++ ) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
               (jj - jmin_dcoeff) * (jp_dcoeff + 1) +
               (kk - kmin_dcoeff) * (kp_dcoeff + dcoeff_gbox.numberCells(1) );

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - 1;

            double phi = average( ptr_phi[idx_pf], ptr_phi[idxm1_pf] );
            double hphi =
               FORT_INTERP_FUNC(
                  phi,
                  d_diff_interp_func_type.c_str() );
            
            double heta = 0.0;
            if ( d_with_third_phase ) {
               double eta =
                  average( ptr_eta[idx_pf], ptr_eta[idxm1_pf] );
               heta =
                  FORT_INTERP_FUNC(
                     eta,
                     d_diff_interp_func_type.c_str() );
            }

            for(int ic=0;ic<nc2;ic++){
               ptr_dx_l[ic][idx_dcoeff] = (1.-hphi)*ptr_dx_coeff_l[ic][idx_dcoeff];
               ptr_dx_a[ic][idx_dcoeff] = hphi*( 1.0 - heta )*ptr_dx_coeff_a[ic][idx_dcoeff];
               if ( d_with_third_phase ) {
                  ptr_dx_b[ic][idx_dcoeff] = hphi*heta*ptr_dx_coeff_b[ic][idx_dcoeff];
               }
            }
         }
      }
   }

   // Y-side
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax+1; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
               (jj - jmin_dcoeff) * jp_dcoeff +
               (kk - kmin_dcoeff) * (kp_dcoeff + jp_dcoeff);

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - jp_pf;

            double phi = 
               average( ptr_phi[idx_pf], ptr_phi[idxm1_pf] );
            double hphi =
               FORT_INTERP_FUNC(
                  phi,
                  d_diff_interp_func_type.c_str() );
            
            double heta = 0.0;
            if ( d_with_third_phase ) {
               double eta =
                  average( ptr_eta[idx_pf], ptr_eta[idxm1_pf] );
               heta =
                  FORT_INTERP_FUNC(
                     eta,
                     d_diff_interp_func_type.c_str() );
            }

            for(int ic=0;ic<nc2;ic++){
               ptr_dy_l[ic][idx_dcoeff] = (1.-hphi)*ptr_dy_coeff_l[ic][idx_dcoeff];
               ptr_dy_a[ic][idx_dcoeff] = hphi*( 1.0 - heta )*ptr_dy_coeff_a[ic][idx_dcoeff];
               if ( d_with_third_phase ) {
                  ptr_dy_b[ic][idx_dcoeff] = hphi*heta*ptr_dy_coeff_b[ic][idx_dcoeff];
               }
            }
         }
      }
   }

   if ( NDIM > 2 ) {
      // Z-side
      for ( int kk = kmin; kk <= kmax+1; kk++ ) {
         for ( int jj = jmin; jj <= jmax; jj++ ) {
            for ( int ii = imin; ii <= imax; ii++ ) {

               const int idx_dcoeff = (ii - imin_dcoeff) +
                  (jj - jmin_dcoeff) * jp_dcoeff +
                  (kk - kmin_dcoeff) * kp_dcoeff;

               const int idx_pf = (ii - imin_pf) +
                  (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;
               const int idxm1_pf = idx_pf - kp_pf;

               double phi = 
                  average( ptr_phi[idx_pf], ptr_phi[idxm1_pf] );
               double hphi =
                  FORT_INTERP_FUNC(
                     phi,
                     d_diff_interp_func_type.c_str() );

               double heta = 0.0;
               if ( d_with_third_phase ) {
                  double eta =
                     average( ptr_eta[idx_pf], ptr_eta[idxm1_pf] );
                  heta =
                     FORT_INTERP_FUNC(
                        eta,
                        d_diff_interp_func_type.c_str() );
               }

               for(int ic=0;ic<nc2;ic++){
                  ptr_dz_l[ic][idx_dcoeff] = (1.-hphi)*ptr_dz_coeff_l[ic][idx_dcoeff];
                  ptr_dz_a[ic][idx_dcoeff] = hphi*( 1.0 - heta )*ptr_dz_coeff_a[ic][idx_dcoeff];
                  if ( d_with_third_phase ) {
                     ptr_dz_b[ic][idx_dcoeff] = hphi*heta*ptr_dz_coeff_b[ic][idx_dcoeff];
                  }
               }
            }
         }
      }
   }  // if ( NDIM > 2 )
}

