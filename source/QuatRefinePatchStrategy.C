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
#include "QuatRefinePatchStrategy.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include <boost/make_shared.hpp>

#include <cassert>

QuatRefinePatchStrategy::QuatRefinePatchStrategy(
   const string& object_name,
   boost::shared_ptr< tbox::Database > input_bc_db,
   const int phase_id,
   const int eta_id,
   const int quat_id,
   const int conc_id,
   const int temperature_id )
   : xfer::RefinePatchStrategy(tbox::Dimension(NDIM)),
     d_object_name( object_name ),
     d_phase_id( phase_id ),
     d_eta_id( eta_id ),
     d_quat_id( quat_id ),
     d_conc_id( conc_id ),
     d_temperature_id( temperature_id )
{
   assert(!object_name.empty());

   if ( d_phase_id >= 0 ) {
      d_phase_refine_strategy =
         new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM), "PhaseBcHelper" );
      d_phase_refine_strategy->setTargetDataId( d_phase_id );
      boost::shared_ptr<tbox::Database> phase_bc_db =
         input_bc_db->getDatabase( "Phase" );
      d_phase_bc_coefs = new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
         "PhaseBcCoefs", phase_bc_db );
      d_phase_refine_strategy->setCoefImplementation( d_phase_bc_coefs );
   }

   if ( d_eta_id >= 0 ) {
      d_eta_refine_strategy =
         new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM), "EtaBcHelper" );
      d_eta_refine_strategy->setTargetDataId( d_eta_id );
      boost::shared_ptr<tbox::Database> eta_bc_db =
         input_bc_db->getDatabase( "Eta" );
      d_eta_bc_coefs = new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
         "EtaBcCoefs", eta_bc_db );
      d_eta_refine_strategy->setCoefImplementation( d_eta_bc_coefs );
   }

   if ( d_quat_id >= 0 ) {
      d_quat_refine_strategy =
         new CartesianRobinBcHelperWithDepth(tbox::Dimension(NDIM), "QuatBcHelper" );
      d_quat_refine_strategy->setTargetDataId( d_quat_id );
      boost::shared_ptr<tbox::Database> quat_bc_db =
         input_bc_db->getDatabase( "Quat" );
      d_quat_bc_coefs = new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
         "QuatBcCoefs", quat_bc_db );
      d_quat_refine_strategy->setCoefImplementation( d_quat_bc_coefs );
   }

   if ( d_conc_id >= 0 ) {
      d_conc_refine_strategy =
         new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM), "ConcBcHelper" );
      d_conc_refine_strategy->setTargetDataId( d_conc_id );
      boost::shared_ptr<tbox::Database> conc_bc_db =
         input_bc_db->getDatabase( "Conc" );
      d_conc_bc_coefs = new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
         "ConcBcCoefs", conc_bc_db );
      d_conc_refine_strategy->setCoefImplementation( d_conc_bc_coefs );
   }

   if ( d_temperature_id >= 0 ) {
      d_temp_refine_strategy =
         new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM), "TemperatureBcHelper" );
      d_temp_refine_strategy->setTargetDataId( d_temperature_id );
      boost::shared_ptr<tbox::Database> temp_bc_db =
         input_bc_db->getDatabase( "Temperature" );
      d_temp_bc_coefs = new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
         "TemperatureBcCoefs", temp_bc_db );
      d_temp_refine_strategy->setCoefImplementation( d_temp_bc_coefs );
   }

}

//=======================================================================

QuatRefinePatchStrategy::~QuatRefinePatchStrategy()
{
}

//=======================================================================

void QuatRefinePatchStrategy::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill )
{
   if ( d_phase_id >= 0 ) {
      d_phase_refine_strategy->setPhysicalBoundaryConditions(
         patch,
         fill_time,
         ghost_width_to_fill );
   }

   if ( d_eta_id >= 0 ) {
      d_eta_refine_strategy->setPhysicalBoundaryConditions(
         patch,
         fill_time,
         ghost_width_to_fill );
   }

   if ( d_quat_id >= 0 ) {
      d_quat_refine_strategy->setPhysicalBoundaryConditions(
         patch,
         fill_time,
         ghost_width_to_fill );
   }

   if ( d_conc_id >= 0 ) {
      d_conc_refine_strategy->setPhysicalBoundaryConditions(
         patch,
         fill_time,
         ghost_width_to_fill );
   }

   if ( d_temperature_id >= 0 ) {
      d_temp_refine_strategy->setPhysicalBoundaryConditions(
         patch,
         fill_time,
         ghost_width_to_fill );
   }
}

//=======================================================================
//
// Print all class data members to given output stream.

void QuatRefinePatchStrategy::printClassData(ostream& os) const
{
   os << "\nQuatRefinePatchStrategy::printClassData..." << endl;
   os << "QuatRefinePatchStrategy: this = " << (QuatRefinePatchStrategy*)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_phase_id =   " << d_phase_id << endl;
   os << "d_eta_id =   " << d_eta_id << endl;
   os << "d_quat_id =   " << d_quat_id << endl;
   os << "d_conc_id =   " << d_conc_id << endl;
   os << "d_temperature_id = " << d_temperature_id << endl;
}
