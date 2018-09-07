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
#include "SteadyStateTemperatureGaussianSource.h"

#include "TemperatureFACSolver.h"
#include "QuatFort.h"
#include "TemperatureFACOps.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace std;

SteadyStateTemperatureGaussianSource::SteadyStateTemperatureGaussianSource(
   const int temperature_scratch_id,
   const int rhs_id,
   const int weight_id,
   const double thermal_diffusivity,
   const int  cp_id,
   boost::shared_ptr<tbox::Database> heat_source_db,
   boost::shared_ptr<tbox::Database > temperature_sys_solver_database,
   boost::shared_ptr<geom::CartesianGridGeometry > grid_geometry,
   HeatCapacityStrategy* heat_capacity_strategy,
   solv::LocationIndexRobinBcCoefs* bc_coefs ):
      SteadyStateTemperatureStrategy(temperature_scratch_id,
                                     rhs_id,
                                     weight_id,
                                     temperature_sys_solver_database,
                                     bc_coefs),
      d_thermal_diffusivity(thermal_diffusivity),
      d_cp_id(cp_id),
      d_heat_capacity_strategy(heat_capacity_strategy),
      d_grid_geometry(grid_geometry),
      d_verbose(false)
{
   assert( temperature_scratch_id>=0 );
   assert( thermal_diffusivity>0. );
   assert( thermal_diffusivity<1.e15 );
   
   tbox::plog<<"Read Gaussian Heat Source profile parameters..."<<std::endl;
   d_source_peak =
      heat_source_db->getDouble( "max_power" );
   d_source_peak *= 10.; // conversion from [mJ/cm^2] to [pJ/um^2]
   d_standard_dev =
      heat_source_db->getDouble( "standard_deviation" );  
   d_pulse_time =
      heat_source_db->getDouble( "pulse_time" );
   d_pulse_width =
      heat_source_db->getDouble( "pulse_width" );

   heat_source_db->getDoubleArray("center",&d_center0[0],NDIM);
      
   if( heat_source_db->isDouble("velocity") )
      heat_source_db->getDoubleArray("velocity",&d_velocity[0],NDIM);
   else
      for(short i=0;i<NDIM;i++)d_velocity[i]=0.;

   for(short i=0;i<NDIM;i++)d_center[i]=tbox::IEEE::getSignalingNaN();

   tbox::Dimension dim(NDIM);
   const hier::IntVector& one_vec = hier::IntVector::getOne(dim);
   const hier::IntVector periodic = grid_geometry->getPeriodicShift(one_vec);

   for ( int dd = 0; dd < NDIM; dd++ ) {
      d_periodic_flag[dd] = ( periodic[dd] != 0 );
   }
   
   if ( temperature_sys_solver_database->isBool("verbose") ) {
      d_verbose = temperature_sys_solver_database->getBool("verbose");
   }
}

void SteadyStateTemperatureGaussianSource::setCurrentTemperature(
   boost::shared_ptr<hier::PatchHierarchy > patch_hierarchy,
   const double time )
{
   //tbox::pout<<"SteadyStateTemperatureGaussianSource: solve for steady state T"<<endl;

   assert( d_temperature_sys_solver );
   assert( d_heat_capacity_strategy );
   assert( d_thermal_diffusivity>0. );
   
   
   const double factor=exp(-(time-d_pulse_time)*(time-d_pulse_time)/(2.*d_pulse_width*d_pulse_width));
   const double source_peak = factor*d_source_peak;
   
   const double* low = d_grid_geometry->getXLower();
   const double* up = d_grid_geometry->getXUpper();
   
   for(short d=0;d<NDIM;d++){
      d_center[d]=d_center0[d]+time*d_velocity[d];
      if( d_periodic_flag[d] ){
         double ll = (up[d]-low[d]);
         d_periodic_length[d]=ll;
         while( d_center[d]<low[d] )d_center[d]+=ll;
         while( d_center[d]>up[d]  )d_center[d]-=ll;
      }else{
         d_periodic_length[d]=-1.;
      }
   }

   //rescale operator and r.h.s. to have values closer to 1. in solver
   //otherwise physical values for thermal diffusivity lead to 
   //matrix elements ~O(1.e12) which cause convergence issues when r.h.s. is small
   const double scaled_thermal_diffusivity=1.e-6*d_thermal_diffusivity;
   const double scaled_source_peak        =1.e-6*source_peak;

   d_temperature_sys_solver->setOperatorCoefficients(1.,0.,-1.*scaled_thermal_diffusivity);

   d_heat_capacity_strategy->setCurrentValue(patch_hierarchy);
   
   int maxln = patch_hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++ ) {

      boost::shared_ptr<hier::PatchLevel > level = patch_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); p++) {
         boost::shared_ptr<hier::Patch > patch = *p;

         const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(patch->getPatchGeometry()) );
         const double* dx = patch_geom->getDx();
         const double* xlo = patch_geom->getXLower();
         const double* xhi = patch_geom->getXUpper();

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast =  pbox.upper();

         boost::shared_ptr< pdat::CellData<double> > rhs (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_rhs_id) ) );
         boost::shared_ptr< pdat::CellData<double> > cp (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_cp_id) ) );
            
         FORT_INITGAUSSIAN_SOURCE(dx, xlo, xhi,
            ifirst(0), ilast(0), 
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            rhs->getPointer(), rhs->getGhostCellWidth()[0],
            cp->getPointer(),  cp->getGhostCellWidth()[0],
            d_center, d_periodic_length, d_standard_dev, scaled_source_peak);
      }
   }

   this->solveSystem();
   
   if(d_verbose)
   {
      math::HierarchyCellDataOpsReal<double> mathops( patch_hierarchy );
      //tbox::pout<<"max. rhs after solve heat equation="<<mathops.max(d_rhs_id)<<endl;
      tbox::plog<<"max. T after solve heat equation="<<mathops.max(d_temperature_scratch_id)<<endl;
      tbox::plog<<"min. T after solve heat equation="<<mathops.min(d_temperature_scratch_id)<<endl;
   }
}

