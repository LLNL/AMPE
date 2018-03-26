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
#include "PhaseConcentrationsStrategy.h"

PhaseConcentrationsStrategy::PhaseConcentrationsStrategy(
      const int conc_l_id,
      const int conc_a_id,
      const int conc_b_id,
      const bool with_third_phase):
   d_conc_l_id(conc_l_id),
   d_conc_a_id(conc_a_id),
   d_conc_b_id(conc_b_id),
   d_with_third_phase(with_third_phase)
{
   assert( d_conc_l_id>=0 );
   assert( d_conc_a_id>=0 );
}

void PhaseConcentrationsStrategy::computePhaseConcentrations(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int temperature_id,
   const int phase_id,
   const int eta_id,
   const int concentration_id )
{
   //tbox::pout<<"CALPHADFreeEnergyStrategy::computePhaseConcentrations()"<<endl;

   assert( temperature_id >= 0 );
   assert( phase_id >= 0 );
   assert( concentration_id >= 0 );
   assert( d_conc_l_id >= 0 );
   assert( d_conc_a_id >= 0 );
   if ( d_with_third_phase ) {
      assert( eta_id >= 0 );
      assert( d_conc_b_id >= 0 );
   }
   const int maxl = hierarchy->getNumberOfLevels();

   for ( int amr_level = 0; amr_level < maxl; amr_level++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( amr_level );

      for ( hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); p++ ) {
         boost::shared_ptr<hier::Patch > patch =
            *p;

         const hier::Box& pbox = patch->getBox();

         boost::shared_ptr< pdat::CellData<double> > temperature (
            patch->getPatchData( temperature_id ), boost::detail::dynamic_cast_tag());
         assert( temperature );
         
         boost::shared_ptr< pdat::CellData<double> > phi (
            patch->getPatchData( phase_id ), boost::detail::dynamic_cast_tag());
         assert( phi );
         
         boost::shared_ptr< pdat::CellData<double> > eta;
         if ( d_with_third_phase ) {
            eta = boost::dynamic_pointer_cast<pdat::CellData<double>,
                                              hier::PatchData>( patch->getPatchData( eta_id ));
         }

         boost::shared_ptr< pdat::CellData<double> > concentration (
            patch->getPatchData( concentration_id ), boost::detail::dynamic_cast_tag());
         assert( concentration );

         boost::shared_ptr< pdat::CellData<double> > c_l (
            patch->getPatchData( d_conc_l_id ), boost::detail::dynamic_cast_tag());
         assert( c_l );
         
         boost::shared_ptr< pdat::CellData<double> > c_a (
            patch->getPatchData( d_conc_a_id ), boost::detail::dynamic_cast_tag());
         assert( c_a );
         
         boost::shared_ptr< pdat::CellData<double> > c_b;
         if ( d_with_third_phase ) {
            c_b = boost::dynamic_pointer_cast<pdat::CellData<double>,
                                              hier::PatchData>( patch->getPatchData( d_conc_b_id )); 
            assert( c_b );
         }

          computePhaseConcentrationsOnPatch(
            temperature,
            phi,
            eta,
            concentration,
            c_l,
            c_a,
            c_b,
            patch );

      }
   }
}
