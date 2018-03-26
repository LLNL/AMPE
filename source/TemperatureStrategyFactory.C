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
#include "TemperatureStrategyFactory.h"

#include "ScalarTemperatureStrategy.h"
#include "GaussianTemperatureStrategy.h"
#include "ConstantTemperatureStrategy.h"
#include "SteadyStateTemperatureCompositionSource.h"
#include "SteadyStateTemperatureGaussianSource.h"
#include "tools.h"

#include <vector>
using namespace std;
using namespace SAMRAI;

TemperatureStrategyFactory::TemperatureStrategyFactory(const int temperature_id, 
                                                       const int temperature_scratch_id,
                                                       const int conc_id, 
                                                       const int weight_id, 
                                                       const int temperature_rhs_id,
                                                       const int cp_id,
                                                       const double molar_volume, 
                                                       const bool with_concentration,
                                                       boost::shared_ptr<geom::CartesianGridGeometry > grid_geometry,
                                                       HeatCapacityStrategy* heat_capacity_strategy):
   d_temperature_id(temperature_id),
   d_temperature_scratch_id(temperature_scratch_id),
   d_conc_id(conc_id),
   d_weight_id(weight_id),
   d_temperature_rhs_id(temperature_rhs_id),
   d_cp_id(cp_id),
   d_molar_volume(molar_volume),
   d_with_concentration(with_concentration),
   d_grid_geometry(grid_geometry),
   d_heat_capacity_strategy(heat_capacity_strategy)
{
   d_temperature_bc_coefs=0;
}

double TemperatureStrategyFactory::readTemperature0(boost::shared_ptr<tbox::Database> temperature_db,
                                                    QuatModelParameters::TEMPERATURE_TYPE temperature_type)
{
   assert( temperature_db );
   
   double temperature0=-1.;
   if( temperature_type == QuatModelParameters::SCALAR || temperature_type == QuatModelParameters::GAUSSIAN ){
      if ( temperature_db->keyExists( "temperature0" ) ) {
         temperature0 = temperature_db->getDouble( "temperature0" );
         printDeprecated( "temperature0", "temperature" );
      }
      else if ( temperature_db->keyExists( "T_parameter" ) ) {
         temperature0 = temperature_db->getDouble( "T_parameter" );
         printDeprecated( "T_parameter", "temperature" );
      }
      else if ( temperature_db->keyExists( "temperature" ) ) {
         temperature0 = temperature_db->getDouble( "temperature" );
      }
      else if ( temperature_db->keyExists( "initial" ) ) {
         temperature0 = temperature_db->getDouble( "initial" );
      }
      else {
         TBOX_ERROR( "Error in TemperatureStrategyFactory: temperature not specified" );
      }
   }
   return temperature0;
}

TemperatureStrategy* TemperatureStrategyFactory::create(
   boost::shared_ptr<tbox::Database> model_db,
   boost::shared_ptr<tbox::Database> integrator_db,
   const QuatModelParameters& model_parameters)
{
   static TemperatureStrategy* strategy=0;
   QuatModelParameters::TEMPERATURE_TYPE temperature_type;
   
   boost::shared_ptr<tbox::Database> temperature_db;
   if ( model_db->keyExists( "Temperature" ) ) {
      temperature_db = model_db->getDatabase( "Temperature" );
   }else{
      temperature_db = model_db;
   }
   
   if( model_db->keyExists("BoundaryConditions") ) {
      boost::shared_ptr<tbox::Database> bc_db =
         model_db->getDatabase( "BoundaryConditions" );
      boost::shared_ptr<tbox::Database> temperature_bc_db =
         bc_db->getDatabase( "Temperature" );
      d_temperature_bc_coefs
         =new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),"TemperatureBcCoefs", temperature_bc_db );
   }
   
   // create strategy
   if ( model_parameters.isTemperatureUniform() ) { // uniform T across spatial domain
      const double temperature0=readTemperature0(temperature_db,QuatModelParameters::SCALAR);

      strategy = new ScalarTemperatureStrategy( 
            d_temperature_id,
            d_temperature_scratch_id,
            temperature0,
            temperature_db );
   }
   else if( model_parameters.isTemperatureGaussian() ){

      strategy = new GaussianTemperatureStrategy( 
            d_temperature_id,
            d_temperature_scratch_id,
            temperature_db,
            d_grid_geometry );
   }
   else if( model_parameters.with_heat_equation() ){
      
      if ( model_parameters.with_steady_temperature() ) {
         assert( d_temperature_scratch_id>=0 );
         assert( d_temperature_rhs_id>=0 );
         assert( d_temperature_bc_coefs!=0 );
         assert( d_heat_capacity_strategy );
            
         boost::shared_ptr<tbox::Database > temperature_sys_solver_database;
         if ( integrator_db->isDatabase( "TemperatureSysSolver" ) ) {
            temperature_sys_solver_database = integrator_db->getDatabase( "TemperatureSysSolver" );
         }
         if( model_parameters.isHeatSourceCompositionDependent() )
            strategy = 
               new SteadyStateTemperatureCompositionSource(d_temperature_scratch_id,
                                               d_conc_id,
                                               d_temperature_rhs_id,
                                               d_weight_id,
                                               model_parameters.thermal_diffusivity(),
                                               d_cp_id,
                                               model_parameters.T_source(),
                                               temperature_sys_solver_database,
                                               d_heat_capacity_strategy,
                                               d_temperature_bc_coefs);
         else
         {
            boost::shared_ptr<tbox::Database >
               heat_source_db = temperature_db->getDatabase( "HeatSource" );
            strategy = 
               new SteadyStateTemperatureGaussianSource(d_temperature_scratch_id,
                                               d_temperature_rhs_id,
                                               d_weight_id,
                                               model_parameters.thermal_diffusivity(),
                                               d_cp_id,
                                               heat_source_db,
                                               temperature_sys_solver_database,
                                               d_grid_geometry,
                                               d_heat_capacity_strategy,
                                               d_temperature_bc_coefs);
         }
      }else{
         strategy = new ConstantTemperatureStrategy(
                                               d_temperature_id,
                                               d_temperature_scratch_id);
      }

   }
   else {
      strategy = new ConstantTemperatureStrategy(
                                               d_temperature_id,
                                               d_temperature_scratch_id);
   }
   
   return strategy;
}
