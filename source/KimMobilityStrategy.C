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
#include "KimMobilityStrategy.h"

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "QuatModel.h"

using namespace std;

KimMobilityStrategy::KimMobilityStrategy(
   QuatModel* quat_model,
   const int conc_l_id,
   const int conc_s_id,
   const int temp_id,
   const string& energy_interp_func_type,
   const string& conc_interp_func_type,
   boost::shared_ptr<tbox::Database> calphad_db,
   boost::shared_ptr<tbox::Database> newton_db,
   const unsigned ncompositions,
   const double DL, const double Q0):
      SimpleQuatMobilityStrategy(quat_model),
      d_conc_l_id(conc_l_id),
      d_conc_s_id(conc_s_id),
      d_temp_id(temp_id),
      d_ncompositions(ncompositions),
      d_DL(DL),
      d_Q0(Q0)
{
   assert( d_conc_l_id>=0 );
   assert( d_conc_s_id>=0 );
   assert( d_temp_id>=0 );
   assert( d_ncompositions>0 );

   t_compute = tbox::TimerManager::getManager()->
      getTimer("AMPE::KimMobilityStrategy::compute");

   if( ncompositions==1 ){
      d_calphad_fenergy = new
         CALPHADFreeEnergyFunctionsBinary(calphad_db,newton_db,
                                 energy_interp_func_type,
                                 conc_interp_func_type,
                                 false); // no 3rd phase
   }else{
      d_calphad_fenergy = new
         CALPHADFreeEnergyFunctionsTernary(calphad_db,newton_db,
                                 energy_interp_func_type,
                                 conc_interp_func_type);
   }

}

void KimMobilityStrategy::computePhaseMobility(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      int& phase_id,
      int& mobility_id,
      const double time,
      const CACHE_TYPE cache )
{
   (void)phase_id;
   (void)time;
   (void)cache;

   t_compute->start();

   const int maxl = hierarchy->getNumberOfLevels();

   for ( int amr_level = 0; amr_level < maxl; amr_level++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( amr_level );

      for( hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); p++ ){
         boost::shared_ptr<hier::Patch > patch = *p;

         boost::shared_ptr< pdat::CellData<double> > temperature(
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( d_temp_id) ) );
         assert( temperature );

         boost::shared_ptr< pdat::CellData<double> > concl(
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( d_conc_l_id) ) );
         assert( concl );

         boost::shared_ptr< pdat::CellData<double> > concs(
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( d_conc_s_id) ) );
         assert( concs );

         boost::shared_ptr< pdat::CellData<double> > mobility(
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData( mobility_id ) ) );
         assert( mobility );

         update(temperature, concl, concs, mobility, patch);
      }
   }

   t_compute->stop();
}

void KimMobilityStrategy::update(
   boost::shared_ptr< pdat::CellData<double> > cd_te,
   boost::shared_ptr< pdat::CellData<double> > cd_cl,
   boost::shared_ptr< pdat::CellData<double> > cd_cs,
   boost::shared_ptr< pdat::CellData<double> > cd_mob,
   boost::shared_ptr<hier::Patch > patch )
{
   assert( cd_te );
   assert( cd_cl );
   assert( cd_cs );
   assert( cd_mob );

   const hier::Box& pbox ( patch->getBox() );
   int imin[3] = {pbox.lower(0),pbox.lower(1),0};
   int imax[3] = {pbox.upper(0),pbox.upper(1),0};
#if (NDIM == 3)
   imin[2] = pbox.lower(2);
   imax[2] = pbox.upper(2);
#endif

   const hier::Box& temp_gbox = cd_te->getGhostBox();
   int min_te[3] = {temp_gbox.lower(0),temp_gbox.lower(1),0};
   int inc_j_te = temp_gbox.numberCells(0);
#if (NDIM == 3)
   min_te[2] = temp_gbox.lower(2);
   int inc_k_te = inc_j_te * temp_gbox.numberCells(1);
#else
   int inc_k_te = 0;
#endif

   //assume cs and cl have the same number of ghost values
   const hier::Box& ci_gbox = cd_cl->getGhostBox();
   int min_ci[3] = {ci_gbox.lower(0),ci_gbox.lower(1),0};
   int inc_j_ci = ci_gbox.numberCells(0);
#if (NDIM == 3)
   min_ci[2] = ci_gbox.lower(2);
   int inc_k_ci = inc_j_ci * ci_gbox.numberCells(1);
#else
   int inc_k_ci = 0;
#endif

   const hier::Box& mo_gbox = cd_mob->getGhostBox();
   int min_mo[3] = {mo_gbox.lower(0),mo_gbox.lower(1),0};
   int inc_j_mo = mo_gbox.numberCells(0);
#if (NDIM == 3)
   min_mo[2] = mo_gbox.lower(2);
   int inc_k_mo = inc_j_mo * mo_gbox.numberCells(1);
#else
   int inc_k_mo = 0;
#endif

   int idx_te = (imin[0]-min_te[0])
              + (imin[1]-min_te[1])*inc_j_te
              + (imin[2]-min_te[2])*inc_k_te;
   int idx_ci = (imin[0]-min_ci[0])
              + (imin[1]-min_ci[1])*inc_j_ci
              + (imin[2]-min_ci[2])*inc_k_ci;
   int idx_mo = (imin[0]-min_mo[0])
              + (imin[1]-min_mo[1])*inc_j_mo
              + (imin[2]-min_mo[2])*inc_k_mo;

   std::vector<double> phaseconc(2*d_ncompositions); // 2 for two phases
   double* cl=&phaseconc[0];
   double* cs=&phaseconc[d_ncompositions];

   for ( int kk = imin[2]; kk <= imax[2]; kk++ ) {
      for ( int jj = imin[1]; jj <= imax[1]; jj++ ) {
         for ( int ii = imin[0]; ii <= imax[0]; ii++ ) {

            const double temp = cd_te->getPointer()[idx_te];
            for(unsigned ic=0;ic<d_ncompositions;ic++)
               cl[ic] = cd_cl->getPointer(ic)[idx_ci];
            for(unsigned ic=0;ic<d_ncompositions;ic++)
               cs[ic] = cd_cs->getPointer(ic)[idx_ci];

            cd_mob->getPointer()[idx_mo] = evaluateMobility(temp, phaseconc);

            idx_te++;
            idx_ci++;
            idx_mo++;
         }
         idx_te+=2*cd_te->getGhostCellWidth()[0];
         idx_ci+=2*cd_cl->getGhostCellWidth()[0];
         idx_mo+=2*cd_mob->getGhostCellWidth()[0];
      }
      idx_te+=2*inc_j_te*cd_te->getGhostCellWidth()[1];
      idx_ci+=2*inc_j_ci*cd_cl->getGhostCellWidth()[1];
      idx_mo+=2*inc_j_mo*cd_mob->getGhostCellWidth()[1];
   }
}

