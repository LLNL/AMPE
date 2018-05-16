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
#include "FuncFort.h"
#include "ConcFort.h"
#include "QuatParams.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFunctions.h"
#include "MolarVolumeStrategy.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;
using namespace std;

#include <cassert>

//=======================================================================

CALPHADFreeEnergyStrategyBinary::CALPHADFreeEnergyStrategyBinary(
   boost::shared_ptr<tbox::Database> calphad_db,
   boost::shared_ptr<tbox::Database> newton_db,
   const string& energy_interp_func_type,
   const string& conc_interp_func_type,
   const string& eta_interp_func_type,
   const string& avg_func_type,
   MolarVolumeStrategy* mvstrategy,
   const int conc_l_id,
   const int conc_a_id,
   const int conc_b_id,
   const bool with_third_phase,
   const double  phase_well_scale,
   const double eta_well_scale,
   const std::string& phase_well_func_type,
   const std::string& eta_well_func_type
   ):
      d_mv_strategy(mvstrategy),
      d_phase_well_scale(phase_well_scale),
      d_eta_well_scale(eta_well_scale),
      d_phase_well_func_type(phase_well_func_type),
      d_eta_well_func_type(eta_well_func_type)
{
   d_with_third_phase = with_third_phase;

   d_energy_interp_func_type = energy_interp_func_type;
   d_conc_interp_func_type = conc_interp_func_type;
   d_eta_interp_func_type = eta_interp_func_type;
   d_avg_func_type = avg_func_type;

   // conversion factor from [J/mol] to [pJ/(mu m)^3]
   // vm^-1 [mol/m^3] * 10e-18 [m^3/(mu m^3)] * 10e12 [pJ/J]
   //d_jpmol2pjpmumcube = 1.e-6 / d_vm;

   // R = 8.314472 J · K-1 · mol-1 
   //tbox::plog << "CALPHADFreeEnergyStrategyBinary:" << endl;
   //tbox::plog << "Molar volume L =" << vml << endl;
   //tbox::plog << "Molar volume A =" << vma << endl;
   //tbox::plog << "jpmol2pjpmumcube=" << d_jpmol2pjpmumcube << endl;
   
   d_conc_l_id = conc_l_id;
   d_conc_a_id = conc_a_id;
   d_conc_b_id = conc_b_id;

   setup(calphad_db,newton_db);
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::setup(
   boost::shared_ptr<tbox::Database> calphad_db,
   boost::shared_ptr<tbox::Database> newton_db)
{
   d_calphad_fenergy = new
      CALPHADFreeEnergyFunctionsBinary(calphad_db,newton_db,
                                 d_energy_interp_func_type,
                                 d_conc_interp_func_type,
                                 d_eta_interp_func_type,d_avg_func_type,
                                 d_with_third_phase,
                                 d_phase_well_scale,d_eta_well_scale,
                                 d_phase_well_func_type,d_eta_well_func_type);
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergyLiquid(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy, 
   const int temperature_id,
   const int fl_id,
   const bool gp ) 
{
   assert( temperature_id >= 0 );
   assert( fl_id >= 0 );
 
   assert( d_conc_l_id >= 0 );

   computeFreeEnergyPrivate(
      hierarchy,
      temperature_id,
      fl_id,
      d_conc_l_id,
      phaseL,
      gp );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeDerivFreeEnergyLiquid(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy, 
   const int temperature_id,
   const int dfl_id ) 
{
   assert( temperature_id >= 0 );
   assert( dfl_id >= 0 );

   computeDerivFreeEnergyPrivate(
      hierarchy,
      temperature_id,
      dfl_id,
      d_conc_l_id,
      phaseL );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergySolidA(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int fs_id,
   const bool gp )
{
   assert( temperature_id >= 0. );
   assert( fs_id >= 0 );

   computeFreeEnergyPrivate(
      hierarchy,
      temperature_id,
      fs_id,
      d_conc_a_id,
      phaseA,
      gp );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeDerivFreeEnergySolidA(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int dfs_id )
{
   assert( temperature_id >= 0. );
   assert( dfs_id >= 0 );

   computeDerivFreeEnergyPrivate(
      hierarchy,
      temperature_id,
      dfs_id,
      d_conc_a_id,
      phaseA );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergySolidB(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int fs_id,
   const bool gp )
{
   assert( temperature_id >= 0. );
   assert( fs_id >= 0 );

   computeFreeEnergyPrivate(
      hierarchy,
      temperature_id,
      fs_id,
      d_conc_b_id,
      phaseB,
      gp );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeDerivFreeEnergySolidB(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int dfs_id )
{
   assert( temperature_id >= 0. );
   assert( dfs_id >= 0 );

   computeDerivFreeEnergyPrivate(
      hierarchy,
      temperature_id,
      dfs_id,
      d_conc_b_id,
      phaseB );
}


//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergyLiquid(
   hier::Patch& patch, 
   const int temperature_id,
   const int fl_id,
   const bool gp ) 
{
   assert( temperature_id >= 0 );
   assert( fl_id >= 0 );
 
   assert( d_conc_l_id >= 0 );

   computeFreeEnergyPrivate(
      patch,
      temperature_id,
      fl_id,
      d_conc_l_id,
      phaseL,
      gp );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergySolidA(
   hier::Patch& patch,
   const int temperature_id,
   const int fs_id,
   const bool gp )
{
   assert( temperature_id >= 0. );
   assert( fs_id >= 0 );

   computeFreeEnergyPrivate(
      patch,
      temperature_id,
      fs_id,
      d_conc_a_id,
      phaseA,
      gp );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergySolidB(
   hier::Patch& patch,
   const int temperature_id,
   const int fs_id,
   const bool gp )
{
   assert( temperature_id >= 0. );
   assert( fs_id >= 0 );

   computeFreeEnergyPrivate(
      patch,
      temperature_id,
      fs_id,
      d_conc_b_id,
      phaseB,
      gp );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergyPrivate(
   hier::Patch& patch,
   const int temperature_id,
   const int f_id,
   const int conc_i_id,
   const PHASE_INDEX pi,
   const bool gp )
{
   assert( temperature_id >= 0 );
   assert( f_id >= 0 );
   assert( conc_i_id >= 0 );
   
   const hier::Box& pbox = patch.getBox();
 
   boost::shared_ptr< pdat::CellData<double> > temperature (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( temperature_id) ) );
 
   boost::shared_ptr< pdat::CellData<double> > f (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( f_id) ) );
   
   boost::shared_ptr< pdat::CellData<double> > c_i (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( conc_i_id) ) );
   
   computeFreeEnergyPrivatePatch(
      pbox,
      temperature,
      f,
      c_i,
      pi,
      gp );
         
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeDerivFreeEnergyPrivate(
   hier::Patch& patch,
   const int temperature_id,
   const int df_id,
   const int conc_i_id,
   const PHASE_INDEX pi )
{
   assert( temperature_id >= 0 );
   assert( df_id >= 0 );
   assert( conc_i_id >= 0 );
   
   const hier::Box& pbox = patch.getBox();
 
   boost::shared_ptr< pdat::CellData<double> > temperature (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( temperature_id) ) );
 
   boost::shared_ptr< pdat::CellData<double> > df (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( df_id) ) );
   
   boost::shared_ptr< pdat::CellData<double> > c_i (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( conc_i_id) ) );
   
   computeDerivFreeEnergyPrivatePatch(
      pbox,
      temperature,
      df,
      c_i,
      pi );
         
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergyPrivate(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int f_id,
   const int conc_i_id,
   const PHASE_INDEX pi,
   const bool gp )
{
   assert( temperature_id >= 0 );
   assert( f_id >= 0 );
   assert( conc_i_id >= 0 );
   
   for ( int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++ ) {
      boost::shared_ptr< hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );
 
      for ( hier::PatchLevel::Iterator ip(level->begin()); ip!=level->end(); ip++ ) {
         boost::shared_ptr<hier::Patch > patch = *ip;
         
         computeFreeEnergyPrivate(*patch, temperature_id, f_id, conc_i_id, pi, gp);
      }
   }
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeDerivFreeEnergyPrivate(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int df_id,
   const int conc_i_id,
   const PHASE_INDEX pi )
{
   assert( temperature_id >= 0 );
   assert( df_id >= 0 );
   assert( conc_i_id >= 0 );
   
   for ( int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++ ) {
      boost::shared_ptr< hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );
 
      for ( hier::PatchLevel::Iterator ip(level->begin()); ip!=level->end(); ip++ ) {
         boost::shared_ptr<hier::Patch > patch = *ip;
         
         computeDerivFreeEnergyPrivate(*patch, temperature_id, df_id, conc_i_id, pi);
      }
   }
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeFreeEnergyPrivatePatch(
   const hier::Box& pbox,
   boost::shared_ptr< pdat::CellData<double> > cd_temp,
   boost::shared_ptr< pdat::CellData<double> > cd_free_energy,
   boost::shared_ptr< pdat::CellData<double> > cd_conc_i,
   const PHASE_INDEX pi,
   const bool gp )
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
         
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;

            const int idx_f = (ii - imin_f) +
               (jj - jmin_f) * jp_f + (kk - kmin_f) * kp_f;

            const int idx_c_i = (ii - imin_c_i) +
               (jj - jmin_c_i) * jp_c_i + (kk - kmin_c_i) * kp_c_i;

            double t = ptr_temp[idx_temp];
            double c_i = ptr_c_i[idx_c_i];

            ptr_f[idx_f] = d_calphad_fenergy->computeFreeEnergy(t,&c_i,pi,gp);
            ptr_f[idx_f] *=d_mv_strategy->computeInvMolarVolume(t,&c_i,pi);
         }
      }
   }
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::computeDerivFreeEnergyPrivatePatch(
   const hier::Box& pbox,
   boost::shared_ptr< pdat::CellData<double> > cd_temp,
   boost::shared_ptr< pdat::CellData<double> > cd_free_energy,
   boost::shared_ptr< pdat::CellData<double> > cd_conc_i,
   const PHASE_INDEX pi )
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
         
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;

            const int idx_f = (ii - imin_f) +
               (jj - jmin_f) * jp_f + (kk - kmin_f) * kp_f;

            const int idx_c_i = (ii - imin_c_i) +
               (jj - jmin_c_i) * jp_c_i + (kk - kmin_c_i) * kp_c_i;

            double t = ptr_temp[idx_temp];
            double c_i = ptr_c_i[idx_c_i];

            d_calphad_fenergy->computeDerivFreeEnergy(t,&c_i,pi, &ptr_f[idx_f]);
            ptr_f[idx_f] *=d_mv_strategy->computeInvMolarVolume(t,&c_i,pi);
         }
      }
   }
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::addComponentRhsPhi(
   const double time,
   hier::Patch& patch,
   const int temperature_id,
   const int phase_id,
   const int eta_id,
   const int conc_id, 
   const int f_l_id,
   const int f_a_id,
   const int f_b_id,
   const int rhs_id )
{
   (void)time;

   assert( conc_id >= 0 );
   assert( phase_id >= 0 );
   assert( f_l_id >= 0 );
   assert( f_a_id >= 0 );
   assert( rhs_id >= 0 );
   assert( temperature_id >= 0 );
   if ( d_with_third_phase ) {
      assert( eta_id >= 0 );
      assert( f_b_id >= 0 );
   }
   assert( d_conc_l_id >=0 );
   assert( d_conc_a_id >=0 );

   boost::shared_ptr< pdat::CellData<double> > phase (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData(phase_id) ) );
   assert( phase );
 
   boost::shared_ptr< pdat::CellData<double> > t (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData( temperature_id) ) );
   assert( t ); 
 
   boost::shared_ptr< pdat::CellData<double> > fl (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData( f_l_id) ) );
   assert( fl );
 
   boost::shared_ptr< pdat::CellData<double> > fa (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData( f_a_id) ) );
   assert( fa );
 
   boost::shared_ptr< pdat::CellData<double> > c_l (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData(
         d_conc_l_id) ) );
   assert( c_l );
 
   boost::shared_ptr< pdat::CellData<double> > c_a (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData(
         d_conc_a_id) ) );
   assert( c_a );

   boost::shared_ptr< pdat::CellData<double> > rhs (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
         patch.getPatchData( rhs_id) ) );

   assert( rhs ); 
   assert( rhs->getGhostCellWidth() == hier::IntVector(tbox::Dimension(NDIM),0) );

   boost::shared_ptr< pdat::CellData<double> > eta;
   boost::shared_ptr< pdat::CellData<double> > fb;
   boost::shared_ptr< pdat::CellData<double> > c_b;
   if ( d_with_third_phase ) {
      eta = BOOST_CAST< pdat::CellData<double>,hier::PatchData>(
               patch.getPatchData( eta_id )); 
      assert( eta );
      fb = BOOST_CAST< pdat::CellData<double>, hier::PatchData> (
              patch.getPatchData( f_b_id ));
      assert( fb );

      c_b = BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch.getPatchData( d_conc_b_id ));
      assert( c_b );
   }
 
   const hier::Box& pbox = patch.getBox();

   addComponentRhsPhiOnPatch(
      rhs,
      t,
      phase,
      eta,
      fl,
      fa,
      fb,
      c_l,
      c_a,
      c_b,
      pbox );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::addComponentRhsPhiOnPatch(
   boost::shared_ptr< pdat::CellData<double> > cd_rhs,
   boost::shared_ptr< pdat::CellData<double> > cd_temperature,
   boost::shared_ptr< pdat::CellData<double> > cd_phi,
   boost::shared_ptr< pdat::CellData<double> > cd_eta,
   boost::shared_ptr< pdat::CellData<double> > cd_f_l,
   boost::shared_ptr< pdat::CellData<double> > cd_f_a,
   boost::shared_ptr< pdat::CellData<double> > cd_f_b,
   boost::shared_ptr< pdat::CellData<double> > cd_c_l,
   boost::shared_ptr< pdat::CellData<double> > cd_c_a,
   boost::shared_ptr< pdat::CellData<double> > cd_c_b,
   const hier::Box& pbox )
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
   
   if ( d_with_third_phase ) {
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
         
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_rhs = (ii - imin_rhs) +
               (jj - jmin_rhs) * jp_rhs + (kk - kmin_rhs) * kp_rhs;

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;

            const int idx_f_i = (ii - imin_f_i) +
               (jj - jmin_f_i) * jp_f_i + (kk - kmin_f_i) * kp_f_i;

            const int idx_c_i = (ii - imin_c_i) +
               (jj - jmin_c_i) * jp_c_i + (kk - kmin_c_i) * kp_c_i;

            double t = ptr_temp[idx_temp];
            double phi = ptr_phi[idx_pf];
            double eta = 0.0;
            double f_l = ptr_f_l[idx_f_i];
            double f_a = ptr_f_a[idx_f_i];
            double f_b = 0.0;
            double c_l = ptr_c_l[idx_c_i];
            double c_a = ptr_c_a[idx_c_i];
            double c_b = 0.0;

            double mu = computeMuA( t, c_a );
            //double mu = 0.;

            double hphi_prime =
               FORT_DERIV_INTERP_FUNC(
                  phi,
                  d_energy_interp_func_type.c_str() );

            double heta = 0.0;

            if ( d_with_third_phase ) {
               eta = ptr_eta[idx_pf];
               f_b = ptr_f_b[idx_f_i];
               c_b = ptr_c_b[idx_c_i];

               heta =
                  FORT_INTERP_FUNC(
                     eta,
                     d_eta_interp_func_type.c_str() );

            }

            ptr_rhs[idx_rhs] +=
               hphi_prime * (
                  ( f_l - ( 1.0 - heta ) * f_a - heta * f_b ) -
                  mu * ( c_l - ( 1.0 - heta ) * c_a - heta * c_b )
                  );

         }
      }
   }
}

//=======================================================================

double CALPHADFreeEnergyStrategyBinary::computeMuA(
   const double t,
   const double c )
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t,&c,phaseA,&mu);
   mu*=d_mv_strategy->computeInvMolarVolume(t,&c,phaseA);

   return mu;
}

//=======================================================================

double CALPHADFreeEnergyStrategyBinary::computeMuL(
   const double t,
   const double c )
{
   double mu;
   d_calphad_fenergy->computeDerivFreeEnergy(t,&c,phaseL,&mu);
   mu*=d_mv_strategy->computeInvMolarVolume(t,&c,phaseL);

   return mu;
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::addComponentRhsEta(
   const double time,
   hier::Patch& patch,
   const int temperature_id,
   const int phase_id,
   const int eta_id,
   const int conc_id, 
   const int f_l_id,
   const int f_a_id,
   const int f_b_id,
   const int rhs_id )
{
   (void)time;

   assert( conc_id >= 0 );
   assert( phase_id >= 0 );
   assert( f_l_id >= 0 );
   assert( f_a_id >= 0 );
   assert( rhs_id >= 0 );
   assert( temperature_id >= 0 );
   assert( eta_id >= 0 );
   assert( f_b_id >= 0 );

   boost::shared_ptr< pdat::CellData<double> > phase (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData(phase_id) ) );
   assert( phase );
 
   boost::shared_ptr< pdat::CellData<double> > eta (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( eta_id) ) );
   assert( eta );

   boost::shared_ptr< pdat::CellData<double> > t (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( temperature_id) ) );
   assert( t ); 
 
   boost::shared_ptr< pdat::CellData<double> > f_l (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( f_l_id) ) );
   assert( f_l );
 
   boost::shared_ptr< pdat::CellData<double> > f_a (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( f_a_id) ) );
   assert( f_a );
 
   boost::shared_ptr< pdat::CellData<double> > f_b (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( f_b_id) ) );
   assert( f_b );

   boost::shared_ptr< pdat::CellData<double> > c_l (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( d_conc_l_id) ) );
   assert( c_l );
 
   boost::shared_ptr< pdat::CellData<double> > c_a (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( d_conc_a_id) ) );
   assert( c_a );
 
   boost::shared_ptr< pdat::CellData<double> > c_b (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( d_conc_b_id) ) );
   assert( c_b );

   boost::shared_ptr< pdat::CellData<double> > rhs (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( rhs_id) ) );
   assert( rhs );
 
   const hier::Box& pbox = patch.getBox();
       
   addComponentRhsEtaOnPatchPrivate(
      rhs,
      t,
      phase,
      eta,
      f_l,
      f_a,
      f_b,
      c_l,
      c_a,
      c_b,
      pbox );
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::addComponentRhsEtaOnPatchPrivate(
   boost::shared_ptr< pdat::CellData<double> > cd_rhs,
   boost::shared_ptr< pdat::CellData<double> > cd_temperature,
   boost::shared_ptr< pdat::CellData<double> > cd_phi,
   boost::shared_ptr< pdat::CellData<double> > cd_eta,
   boost::shared_ptr< pdat::CellData<double> > cd_f_l,
   boost::shared_ptr< pdat::CellData<double> > cd_f_a,
   boost::shared_ptr< pdat::CellData<double> > cd_f_b,
   boost::shared_ptr< pdat::CellData<double> > cd_c_l,
   boost::shared_ptr< pdat::CellData<double> > cd_c_a,
   boost::shared_ptr< pdat::CellData<double> > cd_c_b,
   const hier::Box& pbox )
{
   double* ptr_rhs = cd_rhs->getPointer();

   const double* const ptr_temp = cd_temperature->getPointer();
   const double* const ptr_phi = cd_phi->getPointer();
   const double* const ptr_eta = cd_eta->getPointer();
   const double* const ptr_f_a = cd_f_a->getPointer();
   const double* const ptr_f_b = cd_f_b->getPointer();
   const double* const ptr_c_a = cd_c_a->getPointer();
   const double* const ptr_c_b = cd_c_b->getPointer();
   
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
         
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_rhs = (ii - imin_rhs) +
               (jj - jmin_rhs) * jp_rhs + (kk - kmin_rhs) * kp_rhs;

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;

            const int idx_f_i = (ii - imin_f_i) +
               (jj - jmin_f_i) * jp_f_i + (kk - kmin_f_i) * kp_f_i;

            const int idx_c_i = (ii - imin_c_i) +
               (jj - jmin_c_i) * jp_c_i + (kk - kmin_c_i) * kp_c_i;

            const double t = ptr_temp[idx_temp];
            const double phi = ptr_phi[idx_pf];
            const double eta = ptr_eta[idx_pf];
            const double f_a = ptr_f_a[idx_f_i];
            const double f_b = ptr_f_b[idx_f_i];
            const double c_a = ptr_c_a[idx_c_i];
            const double c_b = ptr_c_b[idx_c_i];

            const double mu = computeMuA( t, c_a );
            //const double mu = 0.;

            const double hphi =
               FORT_INTERP_FUNC(
                  phi,
                  d_energy_interp_func_type.c_str() );

            const double heta_prime =
               FORT_DERIV_INTERP_FUNC(
                  eta,
                  d_eta_interp_func_type.c_str() );

            ptr_rhs[idx_rhs] +=
               hphi * heta_prime * (
                  ( f_a - f_b ) - mu * ( c_a - c_b )
                  );

         }
      }
   }
}


//=======================================================================

void CALPHADFreeEnergyStrategyBinary::defaultComputeSecondDerivativeEnergyPhaseL(
   const double temp,
   const vector<double>& c_l,
   vector<double>& d2fdc2,
   const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(temp,&c_l[0],phaseL,d2fdc2);
   
   if( use_internal_units )
      d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temp,&c_l[0],phaseL);
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::defaultComputeSecondDerivativeEnergyPhaseA(
   const double temp,
   const vector<double>& c_a,
   vector<double>& d2fdc2,
   const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(temp,&c_a[0],phaseA,d2fdc2);
   
   if( use_internal_units )
      d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temp,&c_a[0],phaseA);
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinary::defaultComputeSecondDerivativeEnergyPhaseB(
   const double temp,
   const vector<double>& c_b,
   vector<double>& d2fdc2,
   const bool use_internal_units)
{
   d_calphad_fenergy->computeSecondDerivativeFreeEnergy(temp,&c_b[0],phaseB,d2fdc2);

   if( use_internal_units )
      d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temp,&c_b[0],phaseB);
}
