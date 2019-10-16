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
#include "QuatModel.h"
#include "QuatIntegrator.h"
#include "SimpleGradStrategy.h"
#include "SimpleQuatGradStrategy.h"
#include "TemperatureFreeEnergyStrategy.h"
#include "HBSMFreeEnergyStrategy.h"
#include "KKSdiluteBinary.h"
#include "KKSdiluteEquilibriumPhaseConcentrationsStrategy.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFreeEnergyStrategyTernary.h"
#include "CALPHADFreeEnergyStrategyWithPenalty.h"
#include "KKSCompositionRHSStrategy.h"
#include "EBSCompositionRHSStrategy.h"
#include "BeckermannCompositionRHSStrategy.h"
#include "ConstantTemperatureStrategy.h"
#include "SteadyStateTemperatureStrategy.h"
#include "FuncFort.h"
#include "QuatFort.h"
#include "QuatLinearRefine.h"
#include "QuatWeightedAverage.h"
#include "MinIntCoarsen.h"
#include "ConcFort.h"
#include "ConstantMolarVolumeStrategy.h"
#include "TemperatureStrategyFactory.h"
#include "ConstantHeatCapacityStrategy.h"
#include "NKRHeatCapacityStrategy.h"
#include "PhaseFluxStrategySimple.h"
#include "PhaseFluxStrategyIsotropic.h"
#include "PhaseFluxStrategyAnisotropy.h"
#include "BiasDoubleWellUTRCFreeEnergyStrategy.h"
#include "BiasDoubleWellBeckermannFreeEnergyStrategy.h"
#include "DeltaTemperatureFreeEnergyStrategy.h"
#include "CALPHADequilibriumPhaseConcentrationsStrategy.h"
#include "HBSMequilibriumPhaseConcentrationsStrategy.h"
#include "PartitionPhaseConcentrationsStrategy.h"
#include "PhaseIndependentConcentrationsStrategy.h"
#include "AzizPartitionCoefficientStrategy.h"
#include "UniformPartitionCoefficientStrategy.h"
#include "ConstantMeltingTemperatureStrategy.h"
#include "LinearMeltingTemperatureStrategy.h"
#include "QuatIntegratorFactory.h"
#include "CompositionStrategyMobilities.h"
#include "DiffusionForConcInPhaseStrategy.h"
#include "TbasedCompositionDiffusionStrategy.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "toolsSAMRAI.h"
#include "tools.h"
#include "MobilityFactory.h"
#include "CompositionDiffusionStrategyFactory.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"

#include "PhysicalConstants.h"

#include <set>
#include <map>

#include <boost/lexical_cast.hpp>

using namespace std;

#ifdef HAVE_NETCDF4
using namespace netCDF;
#endif

const double um2tom2 = 1.e-12;

static const int NC_ERR = 2;

QuatModel::QuatModel( int ql ) :
   d_qlen( ql ),
   d_ncompositions( -1 )
{
   d_symmetry_aware = false;

   d_integrator.reset();
   d_integrator_quat_only.reset();
   d_quat_grad_strategy = nullptr;
   d_mobility_strategy = nullptr;
   d_all_refine_patch_strategy = nullptr;
   d_partition_coeff_refine_patch_strategy = nullptr;
   d_phase_conc_strategy = nullptr;
   d_partition_coeff_strategy = nullptr;

   d_temperature_strategy = nullptr;
   d_temperature_strategy_quat_only = nullptr;
   d_meltingT_strategy = nullptr;
   
   d_composition_strategy_mobilities=nullptr;
   d_composition_rhs_strategy=nullptr;
   d_free_energy_strategy_for_diffusion=nullptr;
 
   d_heat_capacity_strategy = nullptr;
   

   d_test_interval.reset();
   d_fundamental_interval.reset();
   d_scalar_diag_interval.reset();
   d_grain_extend_interval.reset();

   d_phase_id = -1;
   d_eta_id = -1;
   d_quat_id = -1;
   d_quat_relax_id = -1;
   d_conc_id = -1;
   d_conc_l_id = -1;
   d_conc_a_id = -1;
   d_conc_b_id = -1;
   d_conc_l_ref_id = -1;
   d_conc_a_ref_id = -1;
   d_conc_b_ref_id = -1;
   d_phase_scratch_id = -1;
   d_eta_scratch_id = -1;
   d_quat_scratch_id = -1;
   d_conc_scratch_id = -1;
   d_conc_l_scratch_id = -1;
   d_conc_a_scratch_id = -1;
   d_conc_b_scratch_id = -1;
   d_phase_grad_cell_id = -1;
   d_phase_grad_side_id = -1;
   d_phase_diffs_id = -1;
   d_phase_diffs_cell_id = -1;
   d_eta_grad_cell_id = -1;
   d_eta_grad_side_id = -1;
   d_eta_diffs_id = -1;
   d_quat_diffs_id = -1;
   d_quat_diffs_cell_id = -1;
   d_quat_nonsymm_diffs_cell_id = -1;
   d_quat_grad_cell_id = -1;
   d_quat_grad_side_id = -1;
   d_quat_grad_modulus_id = -1;
   d_quat_symm_rotation_id = -1;
   d_quat_symm_rotation_cell_id = -1;
   d_phase_mobility_id = -1;
   d_eta_mobility_id = -1;
   d_quat_mobility_id = -1;
   d_quat_norm_error_id = -1;
   d_weight_id = -1;
   d_work_id = -1;
   d_conc_diffusion_id = -1;
   d_conc_phase_coupling_diffusion_id = -1;
   d_conc_eta_coupling_diffusion_id = -1;
   d_conc_pfm_diffusion_l_id = -1;
   d_conc_pfm_diffusion_a_id = -1;
   d_conc_pfm_diffusion_b_id = -1;
   d_conc_diffusion_coeff_l_id = -1;
   d_conc_diffusion_coeff_a_id = -1;
   d_conc_diffusion_coeff_b_id = -1;
   d_velocity_id = -1;
   d_partition_coeff_id = -1;
   d_partition_coeff_scratch_id = -1;
   d_temperature_id = -1;
   d_temperature_scratch_id = -1;
   d_temperature_rhs_steady_id = -1;
   d_f_l_id = -1;
   d_f_a_id = -1;
   d_f_b_id = -1;
   d_cp_id = -1;
   d_conc_Mq_id = -1;

   d_number_of_grains=-1;
   d_phase_threshold=0.85;
   
   d_use_warm_start = false;

   d_conc_l_var.reset();
   d_conc_l_ref_var.reset();

   double def_val = tbox::IEEE::getSignalingNaN();

   d_time = def_val;

   tbox::RestartManager::getManager()->
      registerRestartItem( "QuatModel", this );

   d_fenergy_diag_filename = "";
   d_grains.reset();

   d_verbosity = new QuatVerbosity();
   PFModel::setVerbosity( d_verbosity );
   
   d_cafe=0;
   
   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_resetGrains_timer =
      tman->getTimer("AMPE::QuatModel::resetGrains()");

   t_phase_diffs_timer =
      tman->getTimer("AMPE::QuatModel::phaseDiffs");
}

//=======================================================================

QuatModel::~QuatModel()
{
   delete d_quat_grad_strategy;
   if( d_free_energy_strategy!=d_free_energy_strategy_for_diffusion )
      delete d_free_energy_strategy;
   if( d_free_energy_strategy_for_diffusion!=nullptr )
      delete d_free_energy_strategy_for_diffusion;
   delete d_composition_rhs_strategy;
   delete d_temperature_strategy;
   delete d_temperature_strategy_quat_only;
   if( d_heat_capacity_strategy )
      delete d_heat_capacity_strategy;
   if( d_phase_conc_strategy )delete d_phase_conc_strategy;
   if( d_partition_coeff_strategy )delete d_partition_coeff_strategy;
   d_integrator.reset();
   d_integrator_quat_only.reset();

   delete d_verbosity;
}


//=======================================================================

void QuatModel::initializeTemperature(
   boost::shared_ptr<tbox::Database> model_db,
   boost::shared_ptr<tbox::Database> integrator_db)
{
   if( d_model_parameters.with_heat_equation() ){
      if( d_model_parameters.with_concentration() ){
         tbox::pout << "QuatModel::initializeTemperature() "
                    << "Using NKR model for heat capacity..."
                    << endl;
         d_heat_capacity_strategy = new NKRHeatCapacityStrategy(
            d_model_parameters.cp(), d_cp_id, d_conc_id, d_temperature_id);
      }else{
         const std::map<short,double>& cp( d_model_parameters.cp(0) );
         std::map<short,double>::const_iterator p=cp.find(0);
         const double cpval=p->second;
         d_heat_capacity_strategy = new ConstantHeatCapacityStrategy(cpval, d_cp_id);
      }
   }

   TemperatureStrategyFactory factory(d_temperature_id,d_temperature_scratch_id,
      d_conc_id,
      d_weight_id,
      d_temperature_rhs_steady_id,
      d_cp_id,
      d_model_parameters.molar_volume_liquid(), 
      d_model_parameters.with_concentration(),
      d_grid_geometry,
      d_heat_capacity_strategy);
   
   d_temperature_strategy = factory.create(model_db,integrator_db,d_model_parameters);
   
   d_temperature_strategy_quat_only = new ConstantTemperatureStrategy(
                                               d_temperature_id,
                                               d_temperature_scratch_id);
}

//=======================================================================

void QuatModel::initializeAmr(boost::shared_ptr<tbox::Database> amr_db)
{
   d_use_warm_start = amr_db->getBoolWithDefault( "use_warm_start", false );

   if ( amr_db->isDatabase( "TaggingCriteria" ) ) {
      boost::shared_ptr<tbox::Database> tag_db =
         amr_db->getDatabase( "TaggingCriteria" );
      
      if ( tag_db->isDatabase( "Phi" ) ) {
         d_tag_phase = true;

         boost::shared_ptr<tbox::Database> p_db = tag_db->getDatabase( "Phi" );
      
         d_phase_threshold_tagged = p_db->getDouble( "threshold_tagged" );
         d_phase_threshold_untagged = p_db->getDouble( "threshold_untagged" );
      }

      if ( d_model_parameters.with_third_phase() && tag_db->isDatabase( "Eta" ) ) {
         d_tag_eta = true;
         if ( !d_model_parameters.with_third_phase() ) {
            d_tag_eta = false;
         }

         boost::shared_ptr<tbox::Database> p_db = tag_db->getDatabase( "Eta" );
      
         d_eta_threshold_tagged = p_db->getDouble( "threshold_tagged" );
         d_eta_threshold_untagged = p_db->getDouble( "threshold_untagged" );
      }

      if ( d_model_parameters.with_orientation() && tag_db->isDatabase( "Orient" ) ) {
         d_tag_quat = true;

         boost::shared_ptr<tbox::Database> q_db = tag_db->getDatabase( "Orient" );
      
         d_quat_threshold_tagged = q_db->getDouble( "threshold_tagged" );
         d_quat_threshold_untagged = q_db->getDouble( "threshold_untagged" );
      }
   }
}

//=======================================================================

void QuatModel::initializeCompositionRHSStrategy()
{
   if ( d_model_parameters.concRHSstrategyIsKKS() ){
      d_composition_rhs_strategy =
         new KKSCompositionRHSStrategy(
            d_conc_scratch_id,
            d_phase_scratch_id,
            d_conc_pfm_diffusion_id[0], // use 1x1 diffusion matrix
            d_conc_phase_coupling_diffusion_id,
            d_temperature_scratch_id,
            d_eta_scratch_id,
            d_conc_eta_coupling_diffusion_id,
            d_conc_l_scratch_id,
            d_conc_a_scratch_id,
            d_conc_b_scratch_id,
            d_model_parameters.D_liquid(),
            d_model_parameters.D_solid_A(),
            d_model_parameters.D_solid_B(),
            d_model_parameters.Q0_liquid(),
            d_model_parameters.Q0_solid_A(),
            d_model_parameters.Q0_solid_B(),
            d_model_parameters.energy_interp_func_type(),
            d_model_parameters.avg_func_type() );
   }else if ( d_model_parameters.concRHSstrategyIsEBS() ){
      
      assert( d_diffusion_for_conc_in_phase );
      
      d_composition_rhs_strategy =
         new EBSCompositionRHSStrategy(
            d_phase_scratch_id,
            d_eta_scratch_id,
            static_cast<unsigned short>(d_ncompositions),
            d_conc_l_scratch_id,
            d_conc_a_scratch_id,
            d_conc_b_scratch_id,
            d_temperature_scratch_id,
            d_conc_pfm_diffusion_l_id,
            d_conc_pfm_diffusion_a_id,
            d_conc_pfm_diffusion_b_id,
            d_conc_Mq_id,
            d_model_parameters.Q_heat_transport(),
            d_conc_pfm_diffusion_id,
            d_model_parameters.avg_func_type(),
            d_free_energy_strategy_for_diffusion,
            d_composition_strategy_mobilities,
            d_diffusion_for_conc_in_phase );
   }else if ( d_model_parameters.concRHSstrategyIsBeckermann() ){
      d_composition_rhs_strategy =
         new BeckermannCompositionRHSStrategy(
            this,
            d_conc_scratch_id,
            d_phase_scratch_id,
            d_partition_coeff_scratch_id,
            d_conc_pfm_diffusion_id[0],
            d_conc_phase_coupling_diffusion_id,
            d_model_parameters.D_liquid(),
            d_model_parameters.D_solid_A(),
            d_model_parameters.conc_interp_func_type(),
            d_model_parameters.avg_func_type() );
   }else{
      TBOX_ERROR( "Error: unknown composition RHS Strategy" );
   }
}

//=======================================================================

void QuatModel::initializeRHSandEnergyStrategies(boost::shared_ptr<tbox::MemoryDatabase>& input_db)
{
   tbox::plog<<"QuatModel::initializeRHSandEnergyStrategies()"<<endl;

   assert( d_ncompositions>=0 );
  
   const double Tref = d_model_parameters.with_rescaled_temperature() ? 
                       d_model_parameters.meltingT()/d_model_parameters.rescale_factorT() : 
                       d_model_parameters.meltingT();
 
   boost::shared_ptr<tbox::Database> model_db =
      input_db->getDatabase("ModelParameters");

   double epsilon_anisotropy = d_model_parameters.epsilon_anisotropy();
   
   if( epsilon_anisotropy>=0. )
      d_phase_flux_strategy =
         new PhaseFluxStrategyAnisotropy(d_model_parameters.epsilon_phase(), 
                                         epsilon_anisotropy, 
                                         4);
   else if( d_model_parameters.useIsotropicStencil() ) {
      d_phase_flux_strategy = new PhaseFluxStrategyIsotropic(d_model_parameters.epsilon_phase());
   } else {
      d_phase_flux_strategy = new PhaseFluxStrategySimple(d_model_parameters.epsilon_phase());
   }

   boost::shared_ptr<tbox::MemoryDatabase> calphad_db;
   boost::shared_ptr<tbox::MemoryDatabase> newton_db;
 
   if ( d_model_parameters.with_concentration() ) {
      d_conc_db = model_db->getDatabase( "ConcentrationModel" );

      if ( d_model_parameters.isConcentrationModelCALPHAD() ||
           d_model_parameters.isConcentrationModelKKSdilute() )
      {
         d_mvstrategy = new ConstantMolarVolumeStrategy(
            d_model_parameters.molar_volume_liquid(),
            d_model_parameters.molar_volume_solid_A(),
            d_model_parameters.molar_volume_solid_B());
      }

      // setup free energy strategy first since it may be needed 
      // to setup d_composition_rhs_strategy
      if ( d_model_parameters.isConcentrationModelCALPHAD() ) {
         tbox::pout << "QuatModel: "
                    << "Using CALPHAD model for concentration"
                    << endl;
         d_calphad_db=d_conc_db->getDatabase( "Calphad" );
         std::string calphad_filename = d_calphad_db->getString( "filename" );
         calphad_db.reset ( new tbox::MemoryDatabase( "calphad_db" ) );
         tbox::InputManager::getManager()->parseInputFile(
            calphad_filename, calphad_db );
         
         if ( d_conc_db->isDatabase( "NewtonSolver" ) ){
            d_newton_db = d_conc_db->getDatabase( "NewtonSolver" );
            newton_db.reset ( new tbox::MemoryDatabase( "newton_db" ) );
         }
         
         {
            if( d_ncompositions==1 ){
               d_free_energy_strategy_for_diffusion =
                  new CALPHADFreeEnergyStrategyBinary(
                     calphad_db, newton_db,
                     d_model_parameters.energy_interp_func_type(),
                     d_model_parameters.conc_interp_func_type(),
                     d_mvstrategy,
                     d_conc_l_scratch_id,
                     d_conc_a_scratch_id,
                     d_conc_b_scratch_id,
                     d_model_parameters.with_third_phase());
            }else{
               assert( d_ncompositions==2 );
               d_free_energy_strategy_for_diffusion =
                  new CALPHADFreeEnergyStrategyTernary(
                     calphad_db, newton_db,
                     d_model_parameters.energy_interp_func_type(),
                     d_model_parameters.conc_interp_func_type(),
                     d_mvstrategy,
                     d_conc_l_scratch_id,
                     d_conc_a_scratch_id);
            }
         }

         if( !calphad_db->keyExists( "PenaltyPhaseL" ) ){
             
            d_free_energy_strategy = d_free_energy_strategy_for_diffusion;

            tbox::plog << "QuatModel: "
                       << "CALPHAD with "<<d_ncompositions+1<<" species"<<endl;
            if( d_ncompositions==1 ){
            d_cafe = new CALPHADFreeEnergyFunctionsBinary(
                  calphad_db, newton_db,
                  d_model_parameters.energy_interp_func_type(),
                  d_model_parameters.conc_interp_func_type(),
                  d_model_parameters.with_third_phase());
            }else{
            d_cafe = new CALPHADFreeEnergyFunctionsTernary(
                  calphad_db, newton_db,
                  d_model_parameters.energy_interp_func_type(),
                  d_model_parameters.conc_interp_func_type());
            }
         }else{
            tbox::plog << "QuatModel: "
                       << "Adding penalty to CALPHAD energy"
                       << endl;
            
            assert( d_ncompositions==1 );
            
            d_free_energy_strategy =
               new CALPHADFreeEnergyStrategyWithPenalty(
                  calphad_db, newton_db,
                  d_model_parameters.energy_interp_func_type(),
                  d_model_parameters.conc_interp_func_type(),
                  d_mvstrategy,
                  d_conc_l_scratch_id,
                  d_conc_a_scratch_id,
                  d_conc_b_scratch_id,
                  d_ncompositions,
                  d_model_parameters.with_third_phase() );
            
            d_cafe = new CALPHADFreeEnergyFunctionsWithPenaltyBinary(
                  calphad_db, newton_db,
                  d_model_parameters.energy_interp_func_type(),
                  d_model_parameters.conc_interp_func_type(),
                  d_model_parameters.with_third_phase());
         }
      }// d_model_parameters.isConcentrationModelCALPHAD()
      else if( d_model_parameters.isConcentrationModelKKSdilute() ){
         tbox::pout << "QuatModel: "
                    << "Using KKS dilute model for concentration"
                    << endl;
         d_free_energy_strategy =
            new KKSdiluteBinary(
               d_conc_db,
               d_model_parameters.energy_interp_func_type(),
               d_model_parameters.conc_interp_func_type(),
               d_mvstrategy, 
               d_conc_l_scratch_id,
               d_conc_a_scratch_id );
         d_free_energy_strategy_for_diffusion = d_free_energy_strategy;
      }
      else if( d_model_parameters.isConcentrationModelHBSM() ){
         tbox::pout << "QuatModel: "
                    << "Using HBSM model for concentration"
                    << endl;
         d_free_energy_strategy =
            new HBSMFreeEnergyStrategy(
               d_conc_db->getDatabase( "HBSM" ),
               d_model_parameters.energy_interp_func_type(),
               d_model_parameters.molar_volume_liquid(),
               d_model_parameters.molar_volume_solid_A(),
               d_model_parameters.molar_volume_solid_B(),
               d_model_parameters.D_liquid(),
               d_model_parameters.D_solid_A(),
               d_model_parameters.D_solid_B(),
               d_model_parameters.Q0_liquid(),
               d_model_parameters.Q0_solid_A(),
               d_model_parameters.Q0_solid_B(),
               d_conc_l_scratch_id,
               d_conc_a_scratch_id,
               d_conc_b_scratch_id,
               d_model_parameters.with_third_phase() );
      } // d_model_parameters.isConcentrationModelHBSM()
      else if( d_model_parameters.with_bias_well() ){
         if( d_model_parameters.isConcentrationModelLinear() ){
            d_meltingT_strategy =
               new LinearMeltingTemperatureStrategy(Tref,
                                      d_model_parameters.average_concentration(),
                                      d_model_parameters.liquidus_slope(),
                                      d_conc_l_id,
                                      d_equilibrium_temperature_id);
         
         }else{
            d_meltingT_strategy =
               new ConstantMeltingTemperatureStrategy(Tref,
                                      d_equilibrium_temperature_id);
         }
         if( d_model_parameters.wellBiasBeckermann() ){
            d_free_energy_strategy =
               new BiasDoubleWellBeckermannFreeEnergyStrategy(
                  d_model_parameters.well_bias_alpha(),
                  d_meltingT_strategy );
         }else{
            d_free_energy_strategy =
               new BiasDoubleWellUTRCFreeEnergyStrategy(
                  d_model_parameters.well_bias_alpha(),
                  d_model_parameters.well_bias_gamma(),
                  d_meltingT_strategy );
         }
      }

      if( d_model_parameters.kks_phase_concentration() ){
         tbox::plog<<"Phase concentration determined by KKS"<<endl;
         if ( d_model_parameters.isConcentrationModelCALPHAD() ){
            d_phase_conc_strategy =
               new CALPHADequilibriumPhaseConcentrationsStrategy(
                  d_conc_l_scratch_id,
                  d_conc_a_scratch_id,
                  d_conc_b_scratch_id,
                  d_conc_l_ref_id,
                  d_conc_a_ref_id,
                  d_conc_b_ref_id,
                  d_model_parameters.energy_interp_func_type(),
                  d_model_parameters.conc_interp_func_type(),
                  d_model_parameters.with_third_phase(),
                  calphad_db, newton_db,
                  d_ncompositions );
         }else if( d_model_parameters.isConcentrationModelKKSdilute() ){
            d_phase_conc_strategy =
               new KKSdiluteEquilibriumPhaseConcentrationsStrategy(
                  d_conc_l_scratch_id,
                  d_conc_a_scratch_id,
                  d_conc_b_scratch_id,
                  d_conc_l_ref_id,
                  d_conc_a_ref_id,
                  d_conc_b_ref_id,
                  d_model_parameters.energy_interp_func_type(),
                  d_model_parameters.conc_interp_func_type(),
                  d_conc_db );
         }else{
         if ( d_model_parameters.isConcentrationModelHBSM() )
            d_phase_conc_strategy =
               new HBSMequilibriumPhaseConcentrationsStrategy(
                  d_conc_l_scratch_id,
                  d_conc_a_scratch_id,
                  d_conc_b_scratch_id,
                  d_model_parameters,
                  d_conc_db );
         }
      }else{
         if( d_model_parameters.partition_phase_concentration() ){
            d_phase_conc_strategy =
               new PartitionPhaseConcentrationsStrategy(
                  d_conc_l_scratch_id,
                  d_conc_a_scratch_id,
                  d_conc_b_scratch_id,
                  d_model_parameters.conc_interp_func_type(),
                  d_partition_coeff_id);
         }else{ // simply use cl=ca=c
            d_phase_conc_strategy =
               new PhaseIndependentConcentrationsStrategy(
                  d_conc_l_scratch_id,
                  d_conc_a_scratch_id,
                  d_conc_b_scratch_id);
         }
      }

      math::HierarchyCellDataOpsReal<double> mathops( d_patch_hierarchy );
      mathops.copyData(d_conc_l_id, d_conc_id);
      mathops.copyData(d_conc_a_id, d_conc_id);
      if( d_conc_b_id>-1 )mathops.copyData(d_conc_b_id, d_conc_id);
 
      boost::shared_ptr<tbox::Database> integrator_db =
         input_db->getDatabase("Integrator");
      if( d_model_parameters.concRHSstrategyIsEBS()
            ){
         if( d_model_parameters.conDiffusionStrategyIsCTD() ){
            d_composition_strategy_mobilities =
               new CompositionStrategyMobilities(
                  d_calphad_db,
                  (d_eta_scratch_id>-1),
                  static_cast<unsigned short>(d_ncompositions),
                  d_free_energy_strategy );
         }

         d_diffusion_for_conc_in_phase =
            CompositionDiffusionStrategyFactory::create(
               this, d_model_parameters,
               static_cast<unsigned short>(d_ncompositions),
               d_conc_l_scratch_id,
               d_conc_a_scratch_id,
               d_conc_b_scratch_id,
               d_conc_pfm_diffusion_l_id,
               d_conc_pfm_diffusion_a_id,
               d_conc_pfm_diffusion_b_id,
               d_conc_diffusion_coeff_l_id,
               d_conc_diffusion_coeff_a_id,
               d_conc_diffusion_coeff_b_id,
               d_composition_strategy_mobilities,
               d_free_energy_strategy);
      }

      initializeCompositionRHSStrategy();
 
   } // d_model_parameters.with_concentration()
   else if( d_model_parameters.with_heat_equation() ){
      if( d_model_parameters.with_bias_well() ){
         d_meltingT_strategy =
            new ConstantMeltingTemperatureStrategy(
               Tref, d_equilibrium_temperature_id);
         
         d_free_energy_strategy =
            new BiasDoubleWellUTRCFreeEnergyStrategy(
               d_model_parameters.well_bias_alpha(),
               d_model_parameters.well_bias_gamma(),
               d_meltingT_strategy );
      }else if( d_model_parameters.free_energy_type()[0]=='l' ){
         d_free_energy_strategy =
            new DeltaTemperatureFreeEnergyStrategy(
               Tref,
               d_model_parameters.latent_heat(),
               d_model_parameters.energy_interp_func_type() );
      }else
         d_free_energy_strategy =
            new TemperatureFreeEnergyStrategy(
               d_model_parameters.energy_interp_func_type(),
               d_model_parameters.eta_interp_func_type(),
               d_model_parameters.free_energy_solid_A(),
               d_model_parameters.free_energy_solid_B(),
               d_model_parameters.molar_volume_solid_A(),
               d_model_parameters.molar_volume_solid_B(),
               d_model_parameters.latent_heat(),
               Tref,
               d_model_parameters.with_third_phase() );
   
   } else { // no composition, no heat equation
      if( d_model_parameters.free_energy_type()[0]=='s' ){
         d_free_energy_strategy =
            new PhaseFreeEnergyStrategy(
               d_model_parameters.energy_interp_func_type(),
               d_model_parameters.eta_interp_func_type(),
               d_model_parameters.free_energy_liquid(),
               d_model_parameters.free_energy_solid_A(),
               d_model_parameters.free_energy_solid_B(),
            d_model_parameters.molar_volume_liquid(),
            d_model_parameters.molar_volume_solid_A(),
            d_model_parameters.molar_volume_solid_B(),
            d_model_parameters.with_third_phase() );
      }
   }

   //pure element free energy
   if( d_model_parameters.free_energy_type()[0]=='l' ){
         d_free_energy_strategy =
            new DeltaTemperatureFreeEnergyStrategy(
               Tref,
               d_model_parameters.latent_heat(),
               d_model_parameters.energy_interp_func_type());
   }
   
   if( d_model_parameters.with_Aziz_partition_coeff() ){
      
      // use d_temperature_scratch_id since that field will be uptodate when 
      // integrator needs partion function
      d_partition_coeff_strategy = new AzizPartitionCoefficientStrategy(
         d_velocity_id,
         d_temperature_scratch_id, 
         d_partition_coeff_id,
         d_cafe,
         d_model_parameters.vd(),
         d_model_parameters.keq() );
   }

   if( d_model_parameters.with_uniform_partition_coeff() ){
      
      // use d_temperature_scratch_id since that field will be uptodate when 
      // integrator needs partion function
      d_partition_coeff_strategy = new UniformPartitionCoefficientStrategy(
         d_velocity_id,
         d_temperature_scratch_id, 
         d_partition_coeff_id,
         d_model_parameters.keq() );
   }
   
}

//=======================================================================

void QuatModel::Initialize(
   boost::shared_ptr<tbox::MemoryDatabase>& input_db,
   const string& run_name,
   const bool is_from_restart,
   const string& restart_read_dirname,
   const int restore_num )
{
   tbox::plog<<"QuatModel::Initialize()"<<endl;
   d_is_from_restart=is_from_restart;

   if ( input_db->isDatabase( "FreeEnergyDiagnostics" ) ) {
      boost::shared_ptr<tbox::Database> fenergy_diag_db =
         input_db->getDatabase( "FreeEnergyDiagnostics" );
      
      d_fenergy_diag_filename
         = fenergy_diag_db->getString( "filename" );      
   }

   boost::shared_ptr<tbox::Database> model_db =
      input_db->getDatabase("ModelParameters");

   d_model_parameters.readModelParameters(model_db);

   d_ncompositions=d_model_parameters.ncompositions();

   d_tag_phase = false;
   d_tag_eta = false;
   d_tag_quat = false;

   if ( input_db->isDatabase( "Amr" ) ) {
      boost::shared_ptr<tbox::Database> amr_db = input_db->getDatabase( "Amr" );
      d_amr_enabled = amr_db->getBoolWithDefault( "enabled", true );
      
      initializeAmr(amr_db);
   }
   else {
      d_amr_enabled = false;
   }

   d_test_interval.reset( new EventInterval(
      input_db, "Test", 0.0, "step" ) );

   d_symmetry_aware = false;
   if ( input_db->isDatabase( "Symmetry" ) ) {
      boost::shared_ptr<tbox::Database> symm_db =
         input_db->getDatabase( "Symmetry" );

      if ( symm_db->keyExists( "enabled" ) ) {
         d_symmetry_aware = symm_db->getBoolWithDefault( "enabled", true );
      }
      
      d_fundamental_interval.reset( new EventInterval(
         symm_db, "Fundamental", 0.0, "step" ) );
   }

   d_scalar_diag_interval.reset( new EventInterval(
      input_db, "ScalarDiagnostics", 0.0, "step" ) );

   if ( d_scalar_diag_interval->isActive() ) {
      
      boost::shared_ptr<tbox::Database> tmp_db =
         input_db->getDatabase( "ScalarDiagnostics" );

      d_extra_energy_detail =
         tmp_db->getBoolWithDefault( "extra_energy_detail", false );
   }

   d_grain_extend_interval.reset( new EventInterval(
      input_db, "GrainExtension", 0.0, "step" ) );

   d_ncompositions = d_model_parameters.ncompositionFields();

   // Base class does setup of various common things: logfile,
   // restart, basic samrai objects, visit, ALSO VIRTUAL FUNCTIONS for
   // creating and initializing the integrator and registering
   // variables.


   d_model_parameters.readFreeEnergies(model_db);

   EventInterval tmp_interval( input_db, "Visit", 0.0, "step" );

   if ( tmp_interval.isActive() ) {
      
      boost::shared_ptr<tbox::Database> visit_db =
         input_db->getDatabase( "Visit" );

      d_model_parameters.readVisitOptions(visit_db);
   }

   d_grains.reset( new Grains(d_qlen,
      d_model_parameters.with_visit_grain_output(),
      input_db ) );

   PFModel::Initialize(
      input_db,
      run_name,
      is_from_restart,
      restart_read_dirname,
      restore_num );

   d_grains->initialize(input_db, d_all_periodic);

   // Set up Dirichlet boundary conditions
   if ( ! d_all_periodic ) {
      boost::shared_ptr<tbox::Database> bc_db =
         model_db->getDatabase( "BoundaryConditions" );

      const int phase_id = d_model_parameters.with_phase() ?
                           d_phase_scratch_id : -1; 
      double factor = d_model_parameters.with_rescaled_temperature() ?
                      1./d_model_parameters.rescale_factorT() : -1.;
      d_all_refine_patch_strategy =
         new QuatRefinePatchStrategy(
            "QuatRefinePatchStrategy",
            bc_db,
            phase_id,
            d_eta_scratch_id,
            d_quat_scratch_id,
            d_conc_scratch_id,
            d_temperature_scratch_id, factor );
  
      if( d_model_parameters.needGhosts4PartitionCoeff() ) 
         d_partition_coeff_refine_patch_strategy =
            new PartitionCoeffRefinePatchStrategy(
               "PartitionCoeffRefinePatchStrategy",
               bc_db,
               d_partition_coeff_scratch_id );
   }
   
   boost::shared_ptr<tbox::Database> integrator_db =
      input_db->getDatabase("Integrator");
      
   initializeTemperature(model_db,integrator_db);
   
   if ( ! is_from_restart ) {
      d_number_of_grains = INT_MAX;

      // Read initialization database (initial conditions)
      readInitialDatabase( input_db );

      setupHierarchy();

      // Create single level hierarchy for initial data
      setupInitialDataLevel();

      // Fill single level from initial data file
      FieldsInitializer initializer(d_grid_geometry,
                                    d_ratio_of_init_to_coarsest,
                                    d_verbosity->amrLevel() );

      int temperature_id = d_model_parameters.isTemperatureConstant() ?
                           d_temperature_id : -1;
      int qlen = (d_model_parameters.H_parameter()>=0. ) ? d_qlen : 0;
      initializer.registerFieldsIds(
         d_phase_id,
         d_eta_id,
         temperature_id,
         d_quat_id, qlen,
         d_conc_id, d_ncompositions);

      if( !d_init_c.empty() )initializer.setCvalue(d_init_c);
      if( !d_init_q.empty() )initializer.setQvalue(d_init_q);
      if( d_init_t>=0. )initializer.setTvalue(d_init_t);

      initializer.initializeLevelFromData( d_initial_level,
                                           d_init_data_filename,
                                           d_slice_index );
      //rescale initial conditions for temperature if we are solving
      //time evolution equation for T since we use reduced units
      if( d_model_parameters.with_rescaled_temperature() ){
         assert( d_model_parameters.meltingT()==d_model_parameters.meltingT() );
         math::HierarchyCellDataOpsReal<double> hopscell(d_patch_hierarchy);
 
         hopscell.scale(d_temperature_id,
                        1./d_model_parameters.rescale_factorT(),
                        d_temperature_id);
      }

      for (int ll = 0; ll < d_patch_hierarchy->getMaxNumberOfLevels(); ll++) {
         d_tag_buffer_array[ll] = 1;
      }

      if ( d_model_parameters.isTemperatureUniform() 
        || d_model_parameters.isTemperatureGaussian() 
        || d_model_parameters.isTemperatureGradient() ) {
         d_temperature_strategy->setCurrentTemperature( d_patch_hierarchy, 0.0 );
      }

   }
   else {
      const int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();

      d_patch_hierarchy->initializeHierarchy();
      //d_patch_hierarchy->getFromRestart();

      for (int ll = 0; ll < max_levels; ll++) {
         d_tag_buffer_array[ll] = 1;
      }

      for (int ll = 0; ll < max_levels; ll++) {
         boost::shared_ptr<hier::PatchLevel > null_level_ptr;

         initializeLevelData(
            d_patch_hierarchy,
            ll,
            -1.0,
            false,
            false,
            null_level_ptr,
            true );
      }

      d_gridding_algorithm->getTagAndInitializeStrategy()->
         resetHierarchyConfiguration(
            d_patch_hierarchy,
            0,
            d_patch_hierarchy->getFinestLevelNumber() );

      tbox::RestartManager::getManager()->closeRestartFile();
      
      d_integrator->setTimestep( d_previous_timestep );
   }

   d_quat_grad_strategy = new SimpleQuatGradStrategy( this );

   initializeRHSandEnergyStrategies(input_db);
   
   if ( input_db->isDatabase( "GrainDiagnostics" ) ) {
      boost::shared_ptr<tbox::Database> g_diag_db =
         input_db->getDatabase( "GrainDiagnostics" );
      
      if ( g_diag_db->keyExists( "phase_threshold" ) ) {
         d_phase_threshold = g_diag_db->getDouble( "phase_threshold" );
      }

   }

   InitializeIntegrator();
   
   copyCurrentToScratch(
      d_patch_hierarchy,
      d_time,
      d_all_refine_patch_strategy );

   
}

//=======================================================================

void QuatModel::InitializeIntegrator( void )
{
   tbox::pout<<"QuatModel::InitializeIntegrator()"<<endl;

   assert( d_phase_flux_strategy!=nullptr );
   if( d_model_parameters.with_heat_equation() ){
      assert( d_temperature_strategy );
      assert( d_heat_capacity_strategy );
   }
 
   if ( d_model_parameters.with_phase() )
   d_mobility_strategy  = MobilityFactory::create(
                   this, d_model_parameters,
                   d_conc_l_scratch_id, d_conc_a_scratch_id,
                   d_temperature_scratch_id,
                   d_ncompositions,
                   d_conc_db);
 
   d_integrator->setQuatGradStrategy( d_quat_grad_strategy );
   d_integrator->setMobilityStrategy( d_mobility_strategy );
   d_integrator->setPhaseFluxStrategy( d_phase_flux_strategy );
   if ( d_model_parameters.with_concentration() ){
      d_integrator->setCompositionDiffusionStrategy( d_diffusion_for_conc_in_phase );
      d_integrator->setCompositionRHSStrategy( d_composition_rhs_strategy );
   }
   d_integrator->setFreeEnergyStrategy( d_free_energy_strategy );
   if ( d_model_parameters.with_concentration() )
      d_integrator->setPhaseConcentrationsStrategy( d_phase_conc_strategy );
   if( d_model_parameters.with_partition_coeff() )
      d_integrator->setPartitionCoefficientStrategy( d_partition_coeff_strategy );


   if( d_model_parameters.with_orientation() ){
      d_integrator_quat_only->setQuatGradStrategy( d_quat_grad_strategy );
      d_integrator_quat_only->setMobilityStrategy( d_mobility_strategy );
      d_integrator_quat_only->setFreeEnergyStrategy( d_free_energy_strategy );   
      if ( d_model_parameters.with_concentration() )
         d_integrator_quat_only->setPhaseConcentrationsStrategy( d_phase_conc_strategy );
      d_integrator_quat_only->setTemperatureStrategy( d_temperature_strategy_quat_only );
   }

   d_integrator->setTemperatureStrategy( d_temperature_strategy );

   if( d_model_parameters.with_heat_equation() ){
      d_integrator->setHeatCapacityStrategy( d_heat_capacity_strategy );
   }

   d_integrator->setModelParameters(
      d_time,
      d_end_time,
      d_model_parameters.H_parameter(),
      d_model_parameters.epsilon_phase(),
      d_model_parameters.epsilon_eta(),
      d_model_parameters.epsilon_q(),
      d_model_parameters.quat_grad_floor(),
      d_model_parameters.quat_grad_floor_type(),
      d_model_parameters.phase_well_scale(),
      d_model_parameters.eta_well_scale(),
      d_model_parameters.orient_interp_func_type(),
      d_model_parameters.diffq_avg_func_type(),
      d_model_parameters.phase_well_func_type(),
      d_model_parameters.energy_interp_func_type(),
      d_model_parameters.conc_interp_func_type(),
      d_model_parameters.eta_well_func_type() );

   if( d_model_parameters.with_orientation() )
      d_integrator_quat_only->setModelParameters(
         d_time,
         1.e8, // end_time large enough to never be reached
         d_model_parameters.H_parameter(),
         d_model_parameters.epsilon_phase(),
         d_model_parameters.epsilon_eta(),
         d_model_parameters.epsilon_q(),
         d_model_parameters.quat_grad_floor(),
         "s",
         d_model_parameters.phase_well_scale(),
         d_model_parameters.eta_well_scale(),
         d_model_parameters.orient_interp_func_type(),
         "a", // d_avg_func_type,
         d_model_parameters.phase_well_func_type(),
         d_model_parameters.energy_interp_func_type(),
         d_model_parameters.conc_interp_func_type(),
         d_model_parameters.eta_well_func_type() );

   if ( d_model_parameters.with_concentration() ) {
      d_integrator->setConcentrationModelParameters(
         d_model_parameters.conc_mobility() );
   }

   d_integrator->initialize( d_patch_hierarchy );
   if( d_model_parameters.with_orientation() ){
      const double atol=1.e-4;
      const double rtol=atol*1.e-2;
      tbox::pout << "QuatModel --- "
                 << "set tolerance for Quaternion only integrator:"
                 << atol <<" and "<<rtol
                 << endl;
      d_integrator_quat_only->setAbsTol(atol);
      d_integrator_quat_only->setRelTol(rtol);
      d_integrator_quat_only->initialize( d_patch_hierarchy );
   }
}

//=======================================================================

void QuatModel::initializeRefineCoarsenAlgorithms()
{
   if ( d_model_parameters.with_phase() ) {
      d_phase_refine_op =
         d_grid_geometry->lookupRefineOperator(
            d_temperature_var,
            "LINEAR_REFINE" );
   }
   if ( d_model_parameters.with_third_phase() ) {
      d_eta_refine_op =
         d_grid_geometry->lookupRefineOperator(
            d_eta_var,
            "LINEAR_REFINE" );
   }      

   if ( d_model_parameters.with_orientation() ) {
      if ( d_symmetry_aware ) {
         assert( d_quat_symm_rotation_id>=0 );
         if ( d_verbosity->notSilent() ) {
            tbox::pout << "QuatModel: "
                       << "Using symmetry aware refine/coarsen operators"
                       << endl;
         }
         d_quat_refine_op.reset(
            new QuatLinearRefine( d_quat_symm_rotation_id ) );
         d_quat_coarsen_op.reset(
            new QuatWeightedAverage( true, d_quat_symm_rotation_id ) );
      }
      else {
         d_quat_refine_op =
            d_grid_geometry->lookupRefineOperator(
               d_quat_var,
               "LINEAR_REFINE" );
         d_quat_coarsen_op.reset(
            new QuatWeightedAverage(false) );
      }
   }

   if ( d_model_parameters.with_concentration() ) {
      d_conc_refine_op =
         d_grid_geometry->lookupRefineOperator(
            d_conc_var,
            "LINEAR_REFINE" );
   }

   d_curr_to_curr_refine_alg.reset( new xfer::RefineAlgorithm() );

   d_curr_to_scr_refine_alg.reset( new xfer::RefineAlgorithm() );


   // curr to curr
   if ( d_model_parameters.with_phase() ) {
      d_curr_to_curr_refine_alg->registerRefine(
         d_phase_id,          // destination
         d_phase_id,          // source
         d_phase_scratch_id,  // temporary
         d_phase_refine_op );
   }
   if ( d_model_parameters.with_third_phase() ) {
      d_curr_to_curr_refine_alg->registerRefine(
         d_eta_id,          // destination
         d_eta_id,          // source
         d_eta_scratch_id,  // temporary
         d_eta_refine_op );
   }

   assert( d_temperature_scratch_id >= 0 );
   d_curr_to_curr_refine_alg->registerRefine(
      d_temperature_id,          // destination
      d_temperature_id,          // source
      d_temperature_scratch_id,  // temporary
      d_phase_refine_op );

   if ( d_model_parameters.with_orientation() ){
      assert( d_quat_scratch_id >= 0 );
      d_curr_to_curr_refine_alg->registerRefine(
         d_quat_id,          // destination
         d_quat_id,          // source
         d_quat_scratch_id,  // temporary
         d_quat_refine_op );
   }

   if ( d_model_parameters.with_concentration() ) {
      assert( d_conc_scratch_id >= 0 );
      d_curr_to_curr_refine_alg->registerRefine(
         d_conc_id,          // destination
         d_conc_id,          // source
         d_conc_scratch_id,  // temporary
         d_conc_refine_op );
   }

   // curr to scr
   if ( d_model_parameters.with_phase() ) {
      d_curr_to_scr_refine_alg->registerRefine(
         d_phase_scratch_id,  // destination
         d_phase_id,          // source
         d_phase_scratch_id,  // temporary work space
         d_phase_refine_op );
   }
   if ( d_model_parameters.with_third_phase() ) {
      d_curr_to_scr_refine_alg->registerRefine(
         d_eta_scratch_id,  // destination
         d_eta_id,          // source
         d_eta_scratch_id,  // temporary work space
         d_eta_refine_op );
   }

   d_curr_to_scr_refine_alg->registerRefine(
      d_temperature_scratch_id,  // destination
      d_temperature_id,          // source
      d_temperature_scratch_id,  // temporary work space
      d_phase_refine_op );

   if ( d_model_parameters.with_orientation() ) {
      d_curr_to_scr_refine_alg->registerRefine(
         d_quat_scratch_id,  // destination
         d_quat_id,          // source
         d_quat_scratch_id,  // temporary work space
         d_quat_refine_op );
   }

   if ( d_model_parameters.with_concentration() ) {
      d_curr_to_scr_refine_alg->registerRefine(
         d_conc_scratch_id,          // destination
         d_conc_id,          // source
         d_conc_scratch_id,  // temporary
         d_conc_refine_op );
   }
}


//=======================================================================

void QuatModel::initializeCoarseRefineOperators()
{
   tbox::pout<<"QuatModel::InitializeOperators()"<<endl;

   assert( d_temperature_id >= 0 );
   assert( d_temperature_scratch_id >= 0 );
   assert( d_grains );
   assert( d_grid_geometry );

   initializeRefineCoarsenAlgorithms();

   d_grains->initializeRefineCoarsenAlgorithms(
      d_grid_geometry, d_quat_coarsen_op);
      
   d_integrator->initializeCoarseRefineOperators(
      d_gridding_algorithm,
      d_quat_refine_op,
      d_quat_coarsen_op );
   if( d_model_parameters.with_orientation() )
      d_integrator_quat_only->initializeCoarseRefineOperators(
         d_gridding_algorithm,
         d_quat_refine_op,
         d_quat_coarsen_op );
}

//=======================================================================

void QuatModel::copyCurrentToScratch(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const double time,
   QuatRefinePatchStrategy* patch_strategy )
{
   //tbox::plog<<"QuatModel::copyCurrentToScratch()"<<endl;
   for ( int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++ ) {
      copyCurrentToScratch( hierarchy, ln, time, patch_strategy );
   }
}

void QuatModel::copyCurrentToScratch(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int ln,
   const double time,
   QuatRefinePatchStrategy* patch_strategy )
{
   if ( patch_strategy == d_all_refine_patch_strategy ) {
      d_curr_to_scr_refine_sched[ln]->fillData( time );

   }
   else {
      boost::shared_ptr< hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );

      boost::shared_ptr<xfer::RefineSchedule> schedule(d_curr_to_scr_refine_alg->createSchedule(
         level,
         ln-1,
         hierarchy,
         patch_strategy ) );
      schedule->fillData( time );
   }
}

//=======================================================================

bool QuatModel::isSymmetryAware( void )
{
   return d_symmetry_aware;
}

//=======================================================================

void QuatModel::setupInitialDataLevel( void )
{
   assert( d_temperature_id >= 0 );

   tbox::plog << "\nsetupInitialDataLevel()..." << endl;

   PFModel::setupInitialDataLevel();
   
   assert( d_initial_level );

   if ( d_model_parameters.with_phase() ) {
      if ( !d_initial_level->checkAllocated( d_phase_id ) ) {
         d_initial_level->allocatePatchData( d_phase_id );
         d_initial_level->setTime( 0.0, d_phase_id );
      }
   }

   if ( d_model_parameters.with_third_phase() ) {
      if ( !d_initial_level->checkAllocated( d_eta_id ) ) {
         d_initial_level->allocatePatchData( d_eta_id );
         d_initial_level->setTime( 0.0, d_eta_id );
      }
   }

   if ( !d_initial_level->checkAllocated( d_temperature_id ) ) {
      d_initial_level->allocatePatchData( d_temperature_id );
      d_initial_level->setTime( 0.0, d_temperature_id );
   }

   if ( d_model_parameters.with_orientation() ) {
      if ( !d_initial_level->checkAllocated( d_quat_id ) ) {
         d_initial_level->allocatePatchData( d_quat_id );
         d_initial_level->setTime( 0.0, d_quat_id );
      }
   }

   if ( d_model_parameters.with_concentration() ) {
      if ( !d_initial_level->checkAllocated( d_conc_id ) ) {
         d_initial_level->allocatePatchData( d_conc_id );
         d_initial_level->setTime( 0.0, d_conc_id );
      }
   }
}

//=======================================================================

void QuatModel::setupHierarchy( void )
{
   PFModel::setupHierarchy();
}

//=======================================================================


void QuatModel::CreateIntegrator(
   boost::shared_ptr<tbox::Database> input_db )
{
   tbox::plog<<"QuatModel::CreateIntegrator()"<<endl;

   boost::shared_ptr<tbox::Database> model_db =
      input_db->getDatabase( "ModelParameters" );

   string time_integration = "unsplit";

   time_integration =
      model_db->getStringWithDefault( "time_integration", "unsplit" );

   if ( time_integration == "unsplit" ) {

      d_integrator.reset( QuatIntegratorFactory::create(
         "Integrator",
         d_model_parameters,
         this,
         d_grid_geometry,
         d_qlen,
         d_ncompositions,
         input_db,
         d_use_warm_start,
         d_symmetry_aware,
         d_all_periodic ) );
      if( d_model_parameters.with_orientation() )
         d_integrator_quat_only.reset(
            QuatIntegratorFactory::create(
               "quatonly",
               d_model_parameters,
               this,
               d_grid_geometry,
               d_qlen,
               0,
               input_db,
               false,
               d_symmetry_aware,
               d_all_periodic ) );

   }
   else {
      tbox::pout << time_integration << endl;
      TBOX_ERROR( "Invalid time_integration" << endl );
   }

   d_integrator->setVerbosity( d_verbosity->basicLevel() );
   if( d_model_parameters.with_orientation() )
      d_integrator_quat_only->setVerbosity( d_verbosity->basicLevel() );
}

void QuatModel::registerPhaseConcentrationVariables( )
{
   if( !d_conc_l_var )
   d_conc_l_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "conc_l",
         d_ncompositions ) );
   assert( d_conc_l_var );

   if( !d_conc_a_var )
   d_conc_a_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "conc_a",
         d_ncompositions ) );
   assert( d_conc_a_var );
   
   if ( d_model_parameters.with_third_phase() ) {
      if( !d_conc_b_var )
      d_conc_b_var.reset(
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "conc_b",
            d_ncompositions ) );
      assert( d_conc_b_var );
   }
   
   registerPhaseConcentrationVariables(d_conc_l_var,
                                       d_conc_a_var,
                                       d_conc_b_var);
}

void QuatModel::registerPhaseConcentrationVariables(
   const boost::shared_ptr< pdat::CellVariable<double> > conc_l_var,
   const boost::shared_ptr< pdat::CellVariable<double> > conc_a_var,
   const boost::shared_ptr< pdat::CellVariable<double> > conc_b_var )
{
   assert( conc_l_var );
   assert( conc_a_var );
   
   d_conc_l_var=conc_l_var;
   d_conc_a_var=conc_a_var;
   d_conc_b_var=conc_b_var;

   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::VariableContext> current =
      variable_db->getContext( "CURRENT" );
   boost::shared_ptr<hier::VariableContext> scratch =
      variable_db->getContext( "SCRATCH" );

   //we need internal composition with ghost values for EBS r.h.s.
   //in particular
   d_conc_l_id =
      variable_db->registerVariableAndContext(
         d_conc_l_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_conc_l_scratch_id =
      variable_db->registerVariableAndContext(
         d_conc_l_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
   assert( d_conc_l_id >= 0 );
   assert( d_conc_l_scratch_id >= 0 );

   d_conc_a_id =
      variable_db->registerVariableAndContext(
         d_conc_a_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_conc_a_scratch_id =
      variable_db->registerVariableAndContext(
         d_conc_a_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
   assert( d_conc_a_id >= 0 );
   assert( d_conc_a_scratch_id >= 0 );

   if ( d_model_parameters.with_third_phase() ) {
      assert( d_conc_b_var );
      d_conc_b_id =
         variable_db->registerVariableAndContext(
            d_conc_b_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      d_conc_b_scratch_id =
         variable_db->registerVariableAndContext(
            d_conc_b_var,
            scratch,
            hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
      assert( d_conc_b_id >= 0 );
      assert( d_conc_b_scratch_id >= 0 );
   }
}

void QuatModel::registerConcentrationVariables( void )
{
   tbox::plog<<"QuatModel::registerConcentrationVariables()"<<endl;
   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::VariableContext> current =
      variable_db->getContext( "CURRENT" );
   boost::shared_ptr<hier::VariableContext> scratch =
      variable_db->getContext( "SCRATCH" );
   d_conc_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "concentration", d_ncompositions ) );
   assert( d_conc_var );
   d_conc_id =
      variable_db->registerVariableAndContext(
         d_conc_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_conc_scratch_id =
      variable_db->registerVariableAndContext(
         d_conc_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
   assert( d_conc_id >= 0 );
   assert( d_conc_scratch_id >= 0 );

   {
      // no need for D, just need D_0
      d_conc_diffusion_var.reset();
   }

   for(int ic=0;ic<d_ncompositions;ic++){
      boost::shared_ptr< pdat::SideVariable<double> > conc_pfm_diffusion_var;
      conc_pfm_diffusion_var.reset(
         new pdat::SideVariable<double>(
            tbox::Dimension(NDIM), 
            "conc_pfm_diffusion"+boost::lexical_cast<std::string>(ic) ) );
      assert( conc_pfm_diffusion_var );
      d_conc_pfm_diffusion_var.push_back( conc_pfm_diffusion_var );
      d_conc_pfm_diffusion_id.push_back(
         variable_db->registerVariableAndContext(
            d_conc_pfm_diffusion_var[ic],
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) ) );
      assert( d_conc_pfm_diffusion_id[ic] >= 0 );
   }

   d_model_parameters.checkValidityConcRHSstrategy();
   if( d_model_parameters.concRHSstrategyIsKKS() 
    || d_model_parameters.concRHSstrategyIsBeckermann()){
      d_conc_phase_coupling_diffusion_var.reset(
         new pdat::SideVariable<double>(
            tbox::Dimension(NDIM), "conc_phase_coupling_diffusion" ) );
      assert( d_conc_phase_coupling_diffusion_var );
      d_conc_phase_coupling_diffusion_id =
         variable_db->registerVariableAndContext(
            d_conc_phase_coupling_diffusion_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_conc_phase_coupling_diffusion_id >= 0 );
   }else{
      d_conc_pfm_diffusion_l_var.reset(
         new pdat::SideVariable<double>(
            tbox::Dimension(NDIM), 
            "conc_pfm_diffusion_l", d_ncompositions*d_ncompositions ) );
      assert( d_conc_pfm_diffusion_l_var );
      d_conc_pfm_diffusion_l_id =
         variable_db->registerVariableAndContext(
            d_conc_pfm_diffusion_l_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_conc_pfm_diffusion_l_id >= 0 );

      d_conc_diffusion_coeff_l_var.reset(
         new pdat::SideVariable<double>(
            tbox::Dimension(NDIM),
            "conc_diffusion_coeff_l", d_ncompositions*d_ncompositions ) );
      assert( d_conc_diffusion_coeff_l_var );
      d_conc_diffusion_coeff_l_id =
         variable_db->registerVariableAndContext(
            d_conc_diffusion_coeff_l_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_conc_diffusion_coeff_l_id >= 0 );

      d_conc_pfm_diffusion_a_var.reset(
         new pdat::SideVariable<double>(tbox::Dimension(NDIM), "conc_pfm_diffusion_a",
            d_ncompositions*d_ncompositions ) );
      assert( d_conc_pfm_diffusion_a_var );
      d_conc_pfm_diffusion_a_id =
         variable_db->registerVariableAndContext(
            d_conc_pfm_diffusion_a_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_conc_pfm_diffusion_a_id >= 0 );

      d_conc_diffusion_coeff_a_var.reset(
         new pdat::SideVariable<double>(tbox::Dimension(NDIM), "conc_diffusion_coeff_a", d_ncompositions*d_ncompositions ) );
      assert( d_conc_diffusion_coeff_a_var );
      d_conc_diffusion_coeff_a_id =
         variable_db->registerVariableAndContext(
            d_conc_diffusion_coeff_a_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_conc_diffusion_coeff_a_id >= 0 );

      if( d_model_parameters.with_gradT() )
      {
         assert( d_ncompositions>0 );
         d_conc_Mq_var.reset(
            new pdat::SideVariable<double>(tbox::Dimension(NDIM), "conc_Mq", d_ncompositions ) );
         assert( d_conc_Mq_var );
         d_conc_Mq_id =
            variable_db->registerVariableAndContext(
               d_conc_Mq_var,
               current,
               hier::IntVector(tbox::Dimension(NDIM),0) );
         assert( d_conc_Mq_id >= 0 );
      }

      if ( d_model_parameters.with_third_phase() ) {
         d_conc_pfm_diffusion_b_var.reset(
            new pdat::SideVariable<double>(
               tbox::Dimension(NDIM), "conc_pfm_diffusion_b",
               d_ncompositions*d_ncompositions ) );
         assert( d_conc_pfm_diffusion_b_var );
         d_conc_pfm_diffusion_b_id =
            variable_db->registerVariableAndContext(
               d_conc_pfm_diffusion_b_var,
               current,
               hier::IntVector(tbox::Dimension(NDIM),0) );
         assert( d_conc_pfm_diffusion_b_id >= 0 );

         d_conc_diffusion_coeff_b_var.reset(
            new pdat::SideVariable<double>(
               tbox::Dimension(NDIM), "conc_diffusion_coeff_b",
               d_ncompositions*d_ncompositions ) );
         assert( d_conc_diffusion_coeff_b_var );
         d_conc_diffusion_coeff_b_id =
            variable_db->registerVariableAndContext(
               d_conc_diffusion_coeff_b_var,
               current,
               hier::IntVector(tbox::Dimension(NDIM),0) );
         assert( d_conc_diffusion_coeff_b_id >= 0 );
      }
   }

   if ( d_model_parameters.concentrationModelNeedsPhaseConcentrations() ) {
   
      registerPhaseConcentrationVariables();

   }  // if d_model_parameters.concentrationModelNeedsPhaseConcentrations()
   if ( d_model_parameters.isConcentrationModelCALPHAD() ) {
      d_conc_l_ref_var.reset(
         new pdat::CellVariable<double>(
            tbox::Dimension(NDIM), "conc_l_ref", d_ncompositions ) );
      assert( d_conc_l_ref_var );
      d_conc_a_ref_var.reset(
         new pdat::CellVariable<double>(
            tbox::Dimension(NDIM), "conc_a_ref", d_ncompositions ) );
      assert( d_conc_a_ref_var );
      d_conc_l_ref_id =
         variable_db->registerVariableAndContext(
            d_conc_l_ref_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),1) );
      d_conc_a_ref_id =
         variable_db->registerVariableAndContext(
            d_conc_a_ref_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),1) );
      assert( d_conc_l_ref_id >= 0 );
      assert( d_conc_a_ref_id >= 0 );
      if ( d_model_parameters.with_third_phase() ) {
         d_conc_b_ref_var.reset(
            new pdat::CellVariable<double>(
               tbox::Dimension(NDIM), "conc_b_ref", d_ncompositions ) );
         assert( d_conc_b_ref_var );
         d_conc_b_ref_id =
            variable_db->registerVariableAndContext(
               d_conc_b_ref_var,
               current,
               hier::IntVector(tbox::Dimension(NDIM),1) );
         assert( d_conc_b_ref_id >= 0 );
      }

   }  // if ( d_conc_model == CALPHAD )

   
   if ( d_model_parameters.with_third_phase() ) {
      d_conc_eta_coupling_diffusion_var.reset(
         new pdat::SideVariable<double>(tbox::Dimension(NDIM), "conc_eta_coupling_diffusion" ) );
      assert( d_conc_eta_coupling_diffusion_var );
      d_conc_eta_coupling_diffusion_id =
         variable_db->registerVariableAndContext(
            d_conc_eta_coupling_diffusion_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_conc_eta_coupling_diffusion_id >= 0 );
   }

   if( d_model_parameters.with_partition_coeff() ){
      d_partition_coeff_var.reset (
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "partition_coeff", 1 ));
      d_partition_coeff_id = variable_db->registerVariableAndContext(d_partition_coeff_var,
                                                           current,
                                                           hier::IntVector(tbox::Dimension(NDIM),0));
      assert( d_partition_coeff_id >= 0 );
      
      if( d_model_parameters.needGhosts4PartitionCoeff() ){
         d_partition_coeff_scratch_id = variable_db->registerVariableAndContext(d_partition_coeff_var,
                                                           scratch,
                                                           hier::IntVector(tbox::Dimension(NDIM),1));
         assert( d_partition_coeff_scratch_id >= 0 );
      }      
  }
  if ( d_model_parameters.with_velocity() ){
      d_velocity_var.reset (
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "velocity", 1 ));
      d_velocity_id = variable_db->registerVariableAndContext(d_velocity_var,
                                                           current,
                                                           hier::IntVector(tbox::Dimension(NDIM),0));
      assert( d_velocity_id >= 0 );
   }

   if( d_model_parameters.with_concentration() )
      d_integrator->RegisterConcentrationVariables(
         d_conc_var,
         d_conc_pfm_diffusion_var,
         d_conc_phase_coupling_diffusion_var,
         d_conc_eta_coupling_diffusion_var,
         d_conc_diffusion_var
      );
}

//=======================================================================

void QuatModel::registerEtaVariables( void )
{
   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::VariableContext> current =
      variable_db->getContext( "CURRENT" );
   boost::shared_ptr<hier::VariableContext> scratch =
      variable_db->getContext( "SCRATCH" );

   d_eta_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "eta" ) );
   assert( d_eta_var );
   d_eta_id =
      variable_db->registerVariableAndContext(
         d_eta_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_eta_scratch_id =
      variable_db->registerVariableAndContext(
         d_eta_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

   d_eta_diffs_var.reset(
      new pdat::SideVariable<double>(tbox::Dimension(NDIM), "eta_diffs" ) );
   assert( d_eta_diffs_var );
   d_eta_diffs_id =
      variable_db->registerVariableAndContext(
         d_eta_diffs_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

   d_eta_grad_cell_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "eta_grad_cell", NDIM ) );
   assert( d_eta_grad_cell_var );
   d_eta_grad_cell_id =
      variable_db->registerVariableAndContext(
         d_eta_grad_cell_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   d_eta_grad_side_var.reset(
      new pdat::SideVariable<double>(tbox::Dimension(NDIM), "eta_grad_side", NDIM ) );
   assert( d_eta_grad_side_var );
   d_eta_grad_side_id =
      variable_db->registerVariableAndContext(
         d_eta_grad_side_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   d_eta_mobility_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "eta_mobility" ) );
   assert( d_eta_mobility_var );
   d_eta_mobility_id =
      variable_db->registerVariableAndContext(
         d_eta_mobility_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),1) );
}

//=======================================================================

void QuatModel::registerPhaseVariables( void )
{
   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::VariableContext> current =
      variable_db->getContext( "CURRENT" );
   boost::shared_ptr<hier::VariableContext> scratch =
      variable_db->getContext( "SCRATCH" );

   d_phase_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "phase" ) );
   assert( d_phase_var );
   d_phase_id =
      variable_db->registerVariableAndContext(
         d_phase_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_phase_scratch_id =
      variable_db->registerVariableAndContext(
         d_phase_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

   d_phase_diffs_var.reset(
      new pdat::SideVariable<double>(tbox::Dimension(NDIM), "phase_diffs" ) );
   assert( d_phase_diffs_var );
   d_phase_diffs_id =
      variable_db->registerVariableAndContext(
         d_phase_diffs_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

   if ( d_model_parameters.with_extra_visit_output() ) {
      d_phase_diffs_cell_var.reset(
         new pdat::CellVariable<double>(
            tbox::Dimension(NDIM), "phase_diffs_cell", NDIM ) );
      assert( d_phase_diffs_cell_var );
      d_phase_diffs_cell_id =
         variable_db->registerVariableAndContext(
            d_phase_diffs_cell_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
   }      

   d_phase_grad_cell_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "phase_grad_cell", NDIM ) );
   assert( d_phase_grad_cell_var );
   d_phase_grad_cell_id =
      variable_db->registerVariableAndContext(
         d_phase_grad_cell_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   d_phase_grad_side_var.reset(
      new pdat::SideVariable<double>(tbox::Dimension(NDIM), "phase_grad_side", NDIM ) );
   assert( d_phase_grad_side_var );
   d_phase_grad_side_id =
      variable_db->registerVariableAndContext(
         d_phase_grad_side_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   d_phase_mobility_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "phase_mobility" ) );
   assert( d_phase_mobility_var );
   d_phase_mobility_id =
      variable_db->registerVariableAndContext(
         d_phase_mobility_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),1) );
}

//=======================================================================

void QuatModel::registerOrientationVariables( void )
{
   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::VariableContext> current =
      variable_db->getContext( "CURRENT" );
   boost::shared_ptr<hier::VariableContext> scratch =
      variable_db->getContext( "SCRATCH" );

   const int symm_depth = isSymmetryAware() ? d_qlen * 2 : d_qlen * 1;

   d_quat_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat", d_qlen ) );
   assert( d_quat_var );
   d_quat_id =
      variable_db->registerVariableAndContext(
         d_quat_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_quat_scratch_id =
      variable_db->registerVariableAndContext(
         d_quat_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

   d_quat_relax_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_relax", d_qlen ) );
   assert( d_quat_relax_var );
   d_quat_relax_id =
      variable_db->registerVariableAndContext(
         d_quat_relax_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_quat_relax_scratch_id =
      variable_db->registerVariableAndContext(
         d_quat_relax_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

   d_quat_diffs_var.reset(
      new pdat::SideVariable<double>(tbox::Dimension(NDIM), "quat_diffs", symm_depth ) );
   assert( d_quat_diffs_var );
   d_quat_diffs_id =
      variable_db->registerVariableAndContext(
         d_quat_diffs_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

   d_quat_grad_cell_var.reset(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_grad_cell", NDIM*d_qlen ) );
   assert( d_quat_grad_cell_var );
   d_quat_grad_cell_id =
      variable_db->registerVariableAndContext(
         d_quat_grad_cell_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   d_quat_grad_side_var.reset(
      new pdat::SideVariable<double>(tbox::Dimension(NDIM), "quat_grad_side", NDIM*d_qlen ) );
   assert( d_quat_grad_side_var );
   d_quat_grad_side_id =
      variable_db->registerVariableAndContext(
         d_quat_grad_side_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   if ( d_model_parameters.with_extra_visit_output() ) {
      d_quat_diffs_cell_var.reset(
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_diffs_cell", d_qlen * NDIM ) );
      assert( d_quat_diffs_cell_var );
      d_quat_diffs_cell_id =
         variable_db->registerVariableAndContext(
            d_quat_diffs_cell_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );

      d_quat_nonsymm_diffs_cell_var.reset(
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_nonsymm_diffs_cell", d_qlen * NDIM ) );
      assert( d_quat_nonsymm_diffs_cell_var );
      d_quat_nonsymm_diffs_cell_id =
         variable_db->registerVariableAndContext(
            d_quat_nonsymm_diffs_cell_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );

      d_quat_norm_error_var.reset(
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_norm_error" ) );
      assert( d_quat_norm_error_var );
      d_quat_norm_error_id =
         variable_db->registerVariableAndContext(
            d_quat_norm_error_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );

      if ( d_symmetry_aware ) {
         d_quat_symm_rotation_cell_var.reset(
            new pdat::CellVariable<int>(tbox::Dimension(NDIM), "quat_symm_rotation_cell", NDIM ) );
         assert( d_quat_symm_rotation_cell_var );
         d_quat_symm_rotation_cell_id =
            variable_db->registerVariableAndContext(
               d_quat_symm_rotation_cell_var,
               current,
               hier::IntVector(tbox::Dimension(NDIM),0) );
      }
   }

   if ( d_symmetry_aware ) {
      d_quat_symm_rotation_var.reset(
         new pdat::SideVariable<int>(tbox::Dimension(NDIM), "quat_symm_rotation" ) );
      assert( d_quat_symm_rotation_var );
      d_quat_symm_rotation_id =
         variable_db->registerVariableAndContext(
            d_quat_symm_rotation_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
   }

   d_quat_grad_modulus_var.reset (
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_grad_modulus" ));
   assert( d_quat_grad_modulus_var );
   d_quat_grad_modulus_id =
      variable_db->registerVariableAndContext(
         d_quat_grad_modulus_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   // we have ghost values for quat_mobility so that values are defined
   // at physical boundaries when needed to compute sqrt_mobility
   d_quat_mobility_var.reset (
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_mobility" ));
   assert( d_quat_mobility_var );
   d_quat_mobility_id =
      variable_db->registerVariableAndContext(
         d_quat_mobility_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),1) );

   d_quat_diffusion_var.reset (
      new pdat::SideVariable<double>(tbox::Dimension(NDIM), "quat_diffusion" ));
   assert( d_quat_diffusion_var );
   d_quat_diffusion_id =
      variable_db->registerVariableAndContext(
         d_quat_diffusion_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
}

//=======================================================================

void QuatModel::registerPatchDataForRestart( void )
{
   if ( d_model_parameters.with_phase() ) {
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_phase_id );
   }
   if ( d_model_parameters.with_third_phase() ) {
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_eta_id );
   }

   hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_temperature_id );

   if ( d_model_parameters.with_orientation() ){
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_quat_id );
   }
   if ( d_model_parameters.with_concentration() ){
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_conc_id );
      if ( d_model_parameters.concentrationModelNeedsPhaseConcentrations() )
      {
         hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_conc_l_scratch_id );
         hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_conc_a_scratch_id );
         if ( d_model_parameters.with_third_phase() )
            hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart( d_conc_b_scratch_id );
      }
   }
}

//=======================================================================

void QuatModel::RegisterVariables( void )
{
   tbox::pout<<"QuatModel::RegisterVariables()"<<endl;
   
   assert( d_grains );
   
   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::VariableContext> current =
      variable_db->getContext( "CURRENT" );
   boost::shared_ptr<hier::VariableContext> scratch =
      variable_db->getContext( "SCRATCH" );

   assert( current );
   assert( scratch );

   if ( d_model_parameters.with_phase() ) {
      registerPhaseVariables();
   }
   if ( d_model_parameters.with_third_phase() ) {
      registerEtaVariables();
   }

   d_temperature_var.reset (
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "temperature" ));
   assert( d_temperature_var );
   d_temperature_id =
      variable_db->registerVariableAndContext(
         d_temperature_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   d_temperature_scratch_id =
      variable_db->registerVariableAndContext(
         d_temperature_var,
         scratch,
         hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
   if( d_model_parameters.with_steady_temperature() ){
      d_temperature_rhs_steady_var.reset (
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "temperature_rhs_", 1 ));
      d_temperature_rhs_steady_id =
         variable_db->registerVariableAndContext(
            d_temperature_rhs_steady_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
      assert( d_temperature_rhs_steady_id >= 0 );
   }
   
   if( d_model_parameters.with_bias_well() ){
      d_equilibrium_temperature_var.reset(
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "equilibrium_temperature", 1 ));
      
      d_equilibrium_temperature_id =
         variable_db->registerVariableAndContext(
            d_equilibrium_temperature_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_equilibrium_temperature_id >= 0 );
   }

   if ( d_model_parameters.with_orientation() ) {
      registerOrientationVariables();
   }

   d_weight_var.reset (
      new pdat::CellVariable<double>(tbox::Dimension(NDIM),"weight"));
   assert( d_weight_var );
   d_weight_id =
      variable_db->registerVariableAndContext(
         d_weight_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   d_work_var.reset (
      new pdat::CellVariable<double>(tbox::Dimension(NDIM),"work"));
   assert( d_work_var );
   d_work_id =
      variable_db->registerVariableAndContext(
         d_work_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );

   if( d_model_parameters.with_heat_equation() ){
      d_cp_var.reset (
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "cp" ));
      assert( d_cp_var );
      //ghost value needed in particular for mobility of off-diagonal
      //phase-temperature block in preconditioner
      d_cp_id =
         variable_db->registerVariableAndContext(
            d_cp_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),1) );
      assert( d_cp_id >= 0 );
   }
   
   d_integrator->RegisterVariables(
      d_phase_var,
      d_eta_var,
      d_quat_var,
      d_quat_grad_cell_var,
      d_quat_grad_side_var,
      d_quat_grad_modulus_var,
      d_phase_mobility_var,
      d_eta_mobility_var,
      d_quat_mobility_var,
      d_quat_diffusion_var,
      d_quat_diffs_var,
      d_quat_symm_rotation_var,
      d_weight_var,
      d_temperature_var,
      d_cp_var );
   
   if ( d_model_parameters.with_orientation() ){
      d_integrator_quat_only->RegisterVariables(
         d_phase_var,
         d_eta_var,
         d_quat_relax_var,
         d_quat_grad_cell_var,
         d_quat_grad_side_var,
         d_quat_grad_modulus_var,
         d_phase_mobility_var,
         d_eta_mobility_var,
         d_quat_mobility_var,
         d_quat_diffusion_var,
         d_quat_diffs_var,
         d_quat_symm_rotation_var,
         d_weight_var,
         d_temperature_var,
         d_cp_var );
   }

   //
   // free energy variables
   //
   d_f_l_var.reset (
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "f_l" ));
   assert( d_f_l_var );
   d_f_l_id =
      variable_db->registerVariableAndContext(
         d_f_l_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   assert( d_f_l_id >= 0 );

   d_f_a_var.reset (
      new pdat::CellVariable<double>(tbox::Dimension(NDIM), "f_a" ));
   assert( d_f_a_var );
   d_f_a_id =
      variable_db->registerVariableAndContext(
         d_f_a_var,
         current,
         hier::IntVector(tbox::Dimension(NDIM),0) );
   assert( d_f_a_id >= 0 );

   if ( d_model_parameters.with_third_phase() ) {
      d_f_b_var.reset (
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "f_b" ));
      assert( d_f_b_var );
      d_f_b_id =
         variable_db->registerVariableAndContext(
            d_f_b_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( d_f_b_id >= 0 );
   }

   d_integrator->RegisterFreeEnergyVariables(
      d_f_l_var,
      d_f_a_var,
      d_f_b_var );
   if ( d_model_parameters.with_orientation() )
      d_integrator_quat_only->RegisterFreeEnergyVariables(
         d_f_l_var,
         d_f_a_var,
         d_f_b_var );

   // now create and register concentration variables
   if ( d_model_parameters.with_concentration() ){
      registerConcentrationVariables();
   }

   registerPatchDataForRestart();

   if ( d_model_parameters.with_visit_energy_output() ) {
      d_energy_diag_var.reset (
         new pdat::CellVariable<double>(tbox::Dimension(NDIM), "energy_diag" ));
      assert( d_energy_diag_var );
      d_energy_diag_id =
         variable_db->registerVariableAndContext(
            d_energy_diag_var,
            current,
            hier::IntVector(tbox::Dimension(NDIM),0) );
   }

   d_grains->registerVariables();

   d_integrator->setupBC();

   if ( d_model_parameters.with_orientation() ){
      d_integrator_quat_only->setupBC();
   }
}

//=======================================================================

void QuatModel::RegisterWithVisit( void )
{
   assert( d_visit_data_writer );

   if ( d_model_parameters.with_phase() ) {
      d_visit_data_writer->registerPlotQuantity(
         "phase", "SCALAR", d_phase_id, 0 );
   }
   if ( d_model_parameters.with_third_phase() ) {
      assert( d_eta_id>=0 );
      d_visit_data_writer->registerPlotQuantity(
         "eta", "SCALAR", d_eta_id, 0 );
   }

   if ( !d_model_parameters.isTemperatureUniform() ) {  // if SCALAR, then constant in space
      assert( d_temperature_id >= 0 );
      d_visit_data_writer->registerPlotQuantity(
         "temperature", "SCALAR", d_temperature_id, 0 );
   }

   if ( d_model_parameters.with_orientation()
     && d_model_parameters.evolveQuat() ) {
      assert( d_quat_id >= 0 );
      for ( int n = 0; n < d_qlen; n++ ) {
         string visit_name("q" + tbox::Utilities::intToString(n, 1));

         d_visit_data_writer->registerPlotQuantity(
            visit_name, "SCALAR", d_quat_id, n );
      }
   }

   if ( d_model_parameters.with_concentration() ) {
      assert( d_conc_id >= 0 );
      for ( int n = 0; n < d_ncompositions; n++ ) {
         string visit_name("concentration" + tbox::Utilities::intToString(n, 1));
         d_visit_data_writer->registerPlotQuantity(
            visit_name, "SCALAR", d_conc_id, n );

         if ( d_model_parameters.concentrationModelNeedsPhaseConcentrations()
          && d_model_parameters.with_extra_visit_output() ){
            assert( d_conc_l_scratch_id>=0 );
            assert( d_conc_a_scratch_id>=0 );
            string visit_namel("conc_l" + tbox::Utilities::intToString(n, 1));
            d_visit_data_writer->registerPlotQuantity(
               visit_namel, "SCALAR", d_conc_l_scratch_id, n );
            string visit_namea("conc_a" + tbox::Utilities::intToString(n, 1));
            d_visit_data_writer->registerPlotQuantity(
              visit_namea , "SCALAR", d_conc_a_scratch_id, n );
            if ( d_model_parameters.with_third_phase() ){
               assert( d_conc_b_id>=0 );
               string visit_nameb("conc_b" + tbox::Utilities::intToString(n, 1));
               d_visit_data_writer->registerPlotQuantity(
                  visit_nameb, "SCALAR", d_conc_b_scratch_id, n );
            }
         }
      }
   }
   
   if( d_model_parameters.with_velocity() ){
      assert( d_velocity_id>=0 );
      d_visit_data_writer->registerPlotQuantity(
         "velocity", "SCALAR", d_velocity_id, 0);
   }

   if ( d_model_parameters.with_extra_visit_output() ) {
      if ( d_phase_diffs_cell_id >= 0 ) {
         for ( int d = 0; d < NDIM; d++ ) {
            string visit_name( "phase_diffs_cell"
                             + tbox::Utilities::intToString(d, 1) );
            d_visit_data_writer->registerPlotQuantity(
               visit_name, "SCALAR", d_phase_diffs_cell_id, d );
         }
      }

      if ( d_phase_grad_cell_id >= 0 ) {
         for ( int d=0; d<NDIM; d++ ) {
            string visit_name("phase_grad_cell"
                            + tbox::Utilities::intToString(d, 1));
            d_visit_data_writer->registerPlotQuantity(
               visit_name, "SCALAR", d_phase_grad_cell_id, d );
         }
      }

      if ( d_quat_grad_cell_id >= 0 ) {
         for ( int d=0; d<NDIM; d++ ) {
            for ( int n = 0; n < d_qlen; n++ ) {
               string visit_name(
                  "quat_grad_cell_d" +
                  tbox::Utilities::intToString(d, 1) + "_q" +
                  tbox::Utilities::intToString(n, 1) );

               d_visit_data_writer->registerPlotQuantity(
                  visit_name, "SCALAR", d_quat_grad_cell_id, d*d_qlen + n );
            }
         }
      }

      if ( d_quat_diffs_cell_id >= 0 ) {
         for ( int d = 0; d < NDIM; d++ ) {
            for ( int n = 0; n < d_qlen; n++ ) {
               string visit_symm_name(
                  "quat_diffs_cell_d"
                  + tbox::Utilities::intToString(d, 1) 
                  + "_q"+tbox::Utilities::intToString(n, 1) );

               d_visit_data_writer->registerPlotQuantity(
                  visit_symm_name, "SCALAR",
                  d_quat_diffs_cell_id, d*d_qlen + n );

               string visit_nonsymm_name(
                  "quat_nonsymm_diffs_cell_d"
                  + tbox::Utilities::intToString(d, 1)
                  + "_q"+tbox::Utilities::intToString(n, 1) );

               d_visit_data_writer->registerPlotQuantity(
                  visit_nonsymm_name, "SCALAR",
                  d_quat_nonsymm_diffs_cell_id, d*d_qlen + n );
            }
         }
      }

      if ( d_quat_grad_modulus_id >= 0 ) {
         d_visit_data_writer->registerPlotQuantity(
            "quat_grad_modulus", "SCALAR", d_quat_grad_modulus_id, 0 );
      }

      if ( d_quat_norm_error_id >= 0 ) {
         d_visit_data_writer->registerPlotQuantity(
            "quat_norm_error", "SCALAR", d_quat_norm_error_id, 0 );
      }

      if ( d_phase_mobility_id >= 0 ) {
         d_visit_data_writer->registerPlotQuantity(
            "phase_mobility", "SCALAR", d_phase_mobility_id, 0 );
      }

      if ( d_eta_mobility_id >= 0 ) {
         d_visit_data_writer->registerPlotQuantity(
            "eta_mobility", "SCALAR", d_eta_mobility_id, 0 );
      }

      if ( d_quat_mobility_id >= 0 ) {
         d_visit_data_writer->registerPlotQuantity(
            "quat_mobility", "SCALAR", d_quat_mobility_id, 0 );
      }

      if ( d_model_parameters.with_orientation() && d_symmetry_aware ) {

         for ( int d = 0; d < NDIM; d++ ) {
            string visit_name(
               "quat_symm_rotation_cell" +
               tbox::Utilities::intToString(d, 1) );

            d_visit_data_writer->registerPlotQuantity(
               visit_name, "SCALAR", d_quat_symm_rotation_cell_id, d );
         }
      }
      
      if ( d_model_parameters.with_partition_coeff() ){
         d_visit_data_writer->registerPlotQuantity(
            "partition_coeff", "SCALAR", d_partition_coeff_id, 0);
      }

      if ( d_cp_id>=0 ){
         d_visit_data_writer->registerPlotQuantity(
            "cp", "SCALAR", d_cp_id, 0);
      }
   }  // if ( d_model_parameters.with_extra_visit_output() )

   if ( d_model_parameters.with_visit_grain_output()
     && d_grain_diag_interval->isActive() ) {
      d_grains->registerWithVisit(d_visit_data_writer);
   }

   if ( d_model_parameters.with_visit_energy_output() ) {
      d_visit_data_writer->registerPlotQuantity(
         "energy", "SCALAR", d_energy_diag_id, 0 );
   }

   d_integrator->RegisterWithVisit( d_visit_data_writer );

}

//=======================================================================

void QuatModel::Run( void )
{
   if ( d_model_parameters.with_orientation() && d_symmetry_aware )
      computeSymmetryRotations( d_patch_hierarchy, d_time );
   
   return PFModel::Run();
}

//=======================================================================

void QuatModel::extendGrainOrientation( void )
{
   assert( d_grains );
   assert( d_quat_scratch_id>=0 );
   assert( d_quat_id>=0 );
   assert( d_phase_id>=0 );

   // Fill ghosts of original quat data
   copyCurrentToScratch(
      d_patch_hierarchy,
      d_time,
      d_all_refine_patch_strategy );

   d_grains->extendGrainOrientation(
      d_patch_hierarchy, d_time,
      d_quat_scratch_id, d_phase_id,
      d_quat_id );

   copyCurrentToScratch(
      d_patch_hierarchy,
      d_time,
      d_all_refine_patch_strategy );
}

//=======================================================================

bool QuatModel::resetGrains( void )
{
   assert( d_grains );

   // test if number of grains has changed first...
   int old_number_of_grains = d_number_of_grains;
   d_number_of_grains = d_grains->getNumberOfGrains();
   
   // if the number of grains has not decreased, don't do anything
   if( old_number_of_grains <= d_number_of_grains 
      && !d_grain_extend_interval->hasIntervalPassed(d_cycle, d_time ) )
   {
      return false;
   }
   
   t_resetGrains_timer->start();

   tbox::pout<<"Old number of grains: "<<old_number_of_grains<<endl;
   tbox::pout<<"New number of grains: "<<d_number_of_grains<<endl;
   
   double total_energy, phase_energy;
   double well_energy, free_energy, eta_energy;
   double original_orient_energy;
   double original_qint_energy  ;

   evaluateEnergy(
      d_patch_hierarchy,
      d_time,
      total_energy,
      phase_energy,
      eta_energy,
      original_orient_energy,
      original_qint_energy,
      well_energy,
      free_energy
   );

   if ( d_extra_energy_detail ) {
      tbox::pout << setprecision(8);
      tbox::pout << "  Total energy     = " << total_energy << endl;
      tbox::pout << "    phi energy     = " << phase_energy << endl;
      tbox::pout << "    orient energy  = " << original_orient_energy << endl;
      tbox::pout << "    qint energy    = " << original_qint_energy << endl;
      tbox::pout << "    well energy    = " << well_energy << endl;
      tbox::pout << "    free energy    = " << free_energy << endl;
      if ( d_model_parameters.with_third_phase() ) {
         tbox::pout << "    eta energy     = " << eta_energy << endl;
      }
   }

   extendGrainOrientation();

   if ( d_symmetry_aware )
   {
      if( d_fundamental_interval->isActive() ){
         makeQuatFundamental( d_patch_hierarchy, d_time );
      }

      computeSymmetryRotations( d_patch_hierarchy, d_time );
   }

   double orient_energy=10.e6;
   double qint_energy  =10.e6;
   for(int it=0;it<10;it++){
      smoothQuat( d_patch_hierarchy, d_time );
      evaluateEnergy(
         d_patch_hierarchy,
         d_time,
         total_energy,
         phase_energy,
         eta_energy,
         orient_energy,
         qint_energy,
         well_energy,
         free_energy
      );

      tbox::pout<< setprecision(12);
      tbox::pout<<"Smooth out quaternions, orient energy  = " << orient_energy
                << ", qint energy    = " << qint_energy 
                << ", total quat energy    = " << orient_energy+qint_energy 
                << endl;
   }
   
   math::HierarchyCellDataOpsReal<double> cellops( d_patch_hierarchy );
   cellops.copyData( d_quat_relax_id, d_quat_id, false );
   
   d_integrator_quat_only->initialize( d_patch_hierarchy );

   //double time=d_time;
   double time=0.;
   
   tbox::pout << "Relax quaternions..." << endl;
   for(int it=1;it<200;it++){
      
      const double dt
         = d_integrator_quat_only->Advance(d_patch_hierarchy);
      time += dt;
      
      double dqe=qint_energy;
      double doe=orient_energy;
      
      // diagnostics
      cellops.copyData( d_quat_id, d_quat_relax_id, false );
      evaluateEnergy(
         d_patch_hierarchy,
         d_time,
         total_energy,
         phase_energy,
         eta_energy,
         orient_energy,
         qint_energy,
         well_energy,
         free_energy
      );

      dqe-=qint_energy;
      dqe=fabs(dqe);
      doe-=orient_energy;
      doe=fabs(doe);
      
      tbox::pout<< setprecision(12);
      tbox::pout<<"Smooth out quaternions with dt="<<dt
                << ", orient energy  = " << orient_energy
                << ", qint energy    = " << qint_energy 
                << ", total quat energy    = " << orient_energy+qint_energy 
                << endl;
   
      const double tol =0.1;
      if( dqe <tol*phase_energy*dt
       && qint_energy  < original_qint_energy+tol )
      {
         tbox::pout<<"Quaternions converged after "<<it<<" iterations..."<<endl;
         break;
      }
   }

   cellops.copyData( d_quat_id, d_quat_relax_id, false );
   
   if ( d_extra_energy_detail ) {
      tbox::pout << setprecision(8);
      tbox::pout << "  Total energy     = " << total_energy << endl;
      tbox::pout << "    phi energy     = " << phase_energy << endl;
      tbox::pout << "    orient energy  = " << orient_energy << endl;
      tbox::pout << "    qint energy    = " << qint_energy << endl;
      tbox::pout << "    well energy    = " << well_energy << endl;
      tbox::pout << "    free energy    = " << free_energy << endl;
      if ( d_model_parameters.with_third_phase() ) {
         tbox::pout << "    eta energy     = " << eta_energy << endl;
      }
   }

   t_resetGrains_timer->stop();
   
   return true;
}

//=======================================================================

double QuatModel::Advance( void )
{
   bool update_solution = false;

   if ( d_model_parameters.with_concentration() 
     && d_model_parameters.isConcentrationModelCALPHAD())
      resetRefPhaseConcentrations();

   if ( d_model_parameters.with_orientation() ) {
      if (  d_grain_extend_interval->isActive() ){
         update_solution = resetGrains();
      }

      if ( d_symmetry_aware && !update_solution ){
         if ( d_fundamental_interval->hasIntervalPassed( d_cycle, d_time ) ) {
            makeQuatFundamental( d_patch_hierarchy, d_time );
            update_solution = true;
         
            if( d_symmetry_aware )
               computeSymmetryRotations( d_patch_hierarchy, d_time );
         }
      }
   }

   if ( d_test_interval->hasIntervalPassed( d_cycle, d_time ) ) {
      // do some test, e.g.,
      // update_solution = true;
   }

   if ( update_solution ) {
      d_integrator->updateSolution(
         d_patch_hierarchy,
         0,
         d_patch_hierarchy->getFinestLevelNumber() );
   }

   const double dt =
      d_integrator->Advance(d_patch_hierarchy);

   d_time += dt;
   
   if( d_model_parameters.with_heat_equation() )
      d_heat_capacity_strategy->setCurrentValue( d_patch_hierarchy );

   
   return dt;
}

//-----------------------------------------------------------------------

void QuatModel::postAdvanceDiagnostics( void )
{
   if ( d_scalar_diag_interval->hasIntervalPassed(
           d_cycle,
           d_time ) ) {

      printScalarDiagnostics();

   }
}

//-----------------------------------------------------------------------

void QuatModel::preRunDiagnostics( void )
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   math::HierarchyCellDataOpsReal<double> mathops( d_patch_hierarchy );

   double surface_e=d_model_parameters.surfaceEnergy();
   tbox::pout<<"Surface energy (J/m^2): "<<surface_e<<endl;
   double width=d_model_parameters.interfacialWidth();
   tbox::pout<<"Interfacial width (um): "<<width<<endl;

   if ( d_model_parameters.with_concentration() )
      d_composition_rhs_strategy->printDiagnostics( d_patch_hierarchy );

   const double temperature =
      d_temperature_strategy->getCurrentMinTemperature(
         d_patch_hierarchy, d_time );

   if ( d_fenergy_diag_filename != ""  ) {
      
      // pre-run diagnostics
      if ( d_model_parameters.with_concentration() 
        && d_model_parameters.isConcentrationModelCALPHAD() ){
         assert( temperature>0. );
         if( mpi.getRank()==0 ){
            ofstream ffile("FvsT.dat", ios::out);
            assert( d_free_energy_strategy_for_diffusion );
            d_free_energy_strategy_for_diffusion->preRunDiagnostics(ffile);
         }

         if( mpi.getRank()==0 && d_cafe!=0 ){
            //energy vs. composition for phi=0 and phi=1
            d_cafe->printEnergyVsComposition(temperature);
         }

         // compute equilibrium composition for pair L,A
         double phi_min = mathops.min( d_phase_id );
         //double phi_max = mathops.max( d_phase_id );
         
         double ceq[4]; // 2 phases x 2 compositions max.
         bool found_ceq=false;
         if( phi_min<0.1 ){
            found_ceq=computeCeq(temperature,phaseL, phaseA,&ceq[0]);

            // compute equilibrium composition for pair L,B
            if( d_model_parameters.with_third_phase() ){
               found_ceq=computeCeq(temperature,phaseL, phaseB,&ceq[0]);
            }
         }
         
         if( d_model_parameters.with_third_phase() )
         {
            found_ceq=computeCeq(temperature,phaseA, phaseB,&ceq[0]);
         }
         
         if( d_cafe!=0 && found_ceq )
         {
            if( phi_min<0.1 )
               d_cafe->energyVsPhiAndC(
                  temperature, &ceq[0], found_ceq, 
                  d_model_parameters.phase_well_scale(),
                  d_model_parameters.phase_well_func_type(),
                  false);
            if( d_model_parameters.with_third_phase() ){
               d_cafe->energyVsPhiAndC(
                  temperature, &ceq[0], found_ceq,
                  d_model_parameters.phase_well_scale(),
                  d_model_parameters.phase_well_func_type(),
                  true);
            }
         }
         mpi.Barrier();
         
         if( found_ceq && !d_is_from_restart)
         {
            setRefPhaseConcentrationsToEquilibrium(ceq);
            if( d_model_parameters.initPhaseConcAtEq() )
               setPhaseConcentrationsToEquilibrium(ceq);
         }
         
      } // with_concentration

   }
         
   mpi.Barrier();
   
   if ( d_scalar_diag_interval->includeInitial( d_time ) ) {
      printScalarDiagnostics();
   }

   if( d_model_parameters.with_concentration() )
      preRunDiagnosticsMobilityInPhases(temperature);
}

//-----------------------------------------------------------------------

bool QuatModel::computeCeq(const double temperature, 
                           const PHASE_INDEX pi0, const PHASE_INDEX pi1,
                           double* ceq)const
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   math::HierarchyCellDataOpsReal<double> mathops( d_patch_hierarchy );

   double cmin = mathops.min( d_conc_id );
   double cmax = mathops.max( d_conc_id );
   double dc=cmax-cmin;
   cmin = max( 0.25*cmin,    cmin-0.2*dc );
   cmax = min( 1.-0.25*(1.-cmax), cmax+0.2*dc );
   if( cmax-cmin <= 1.e-8 ){
      cmax = min(1.,cmin + 0.1);
      cmin = max(0.,cmin - 0.1);
   }
   tbox::pout<<"QuatModel::computeCeq(): "<<endl
             <<"Try to estimate equilibrium concentrations between cmin="
             <<cmin <<" and cmax="<<cmax<<"..."<<endl;
   tbox::pout<<"T="<<temperature<<endl;

   double ceq_init0=cmin;
   double ceq_init1=cmax;
   double lceq[4]={ceq_init0,ceq_init1,ceq_init0,ceq_init1};

   if( d_model_parameters.knownInitCinPhase() ){
      const unsigned int offset=d_ncompositions;
      for(int ic=0;ic<d_ncompositions;ic++){
         lceq[       ic]=d_model_parameters.getInitCphaseL(ic);
         lceq[offset+ic]=d_model_parameters.getInitCphaseA(ic);
      }
   }
 
   // compute equilibrium concentrations
   bool found_ceq = false;
   if( mpi.getRank()==0 ) // do it on PE0 only to avoid error message prints from all PEs
   {
      found_ceq =
         d_cafe->computeCeqT(temperature,pi0,pi1,&lceq[0], 50, true);
      if( lceq[0]>1. )found_ceq = false;
      if( lceq[0]<0. )found_ceq = false;
      if( lceq[1]>1. )found_ceq = false;
      if( lceq[1]<0. )found_ceq = false;
      
      if( !found_ceq )
      {
         tbox::pout<<"Try again with different initial conditions..."<<endl;
         lceq[0]=ceq_init1;
         lceq[1]=ceq_init0;
         found_ceq =
            d_cafe->computeCeqT(temperature,pi0,pi1,&lceq[0],50,true);
         if( lceq[0]>1. )found_ceq = false;
         if( lceq[0]<0. )found_ceq = false;
         if( lceq[1]>1. )found_ceq = false;
         if( lceq[1]<0. )found_ceq = false;
      }
      
      if( found_ceq ){
         tbox::plog<<"Found equilibrium concentrations: "
                   <<lceq[0]<<", "<<lceq[1]<<"..."<<endl;
         if( d_ncompositions>1 )
            tbox::plog<<"                                  "
                      <<lceq[2]<<", "<<lceq[3]<<"..."<<endl;

         if( d_model_parameters.isConcentrationModelCALPHAD() )
         {
            vector<double> d2fdc2(1);
            CALPHADFreeEnergyFunctionsBinary* cafe=
               dynamic_cast<CALPHADFreeEnergyFunctionsBinary*>(d_cafe);
            cafe->computeSecondDerivativeFreeEnergy(
                    temperature,&lceq[0],pi0,d2fdc2);
            for(vector<double>::const_iterator it=d2fdc2.begin();
                                               it!=d2fdc2.end();
                                             ++it)
            tbox::plog<<"d2fdc2="<<*it<<endl;
            cafe->computeSecondDerivativeFreeEnergy(
                    temperature,&lceq[0],pi1,d2fdc2);
            for(vector<double>::const_iterator it=d2fdc2.begin();
                                               it!=d2fdc2.end();
                                             ++it)
            tbox::plog<<"d2fdc2="<<*it<<endl;
         }
      }else{
         tbox::plog<<"ERROR: Equilibrium concentrations not found... "<<endl;
      }
   
      if( !found_ceq )mpi.abort();
   }
   
   mpi.Bcast( &lceq[0], 4, MPI_DOUBLE, 0);
   
   int flag = (int)found_ceq;
   mpi.Bcast( &flag, 1, MPI_INT, 0);
   found_ceq=(bool)flag;
   
   ceq[0]=lceq[0];
   ceq[1]=lceq[1];
   ceq[2]=lceq[2];
   ceq[3]=lceq[3];
 
   return found_ceq;
}


//-----------------------------------------------------------------------

void QuatModel::preRunDiagnosticsMobilityInPhases( const double temperature )
{  
   if ( d_model_parameters.isConcentrationModelCALPHAD() )
   {
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

      assert( d_calphad_db );
      string calphad_filename = d_calphad_db->getString( "filename" );
 
      boost::shared_ptr<tbox::MemoryDatabase> calphad_db
         ( new tbox::MemoryDatabase( "calphad_db" ) );
      tbox::InputManager::getManager()->parseInputFile(
         calphad_filename, calphad_db );
      
      if( calphad_db->isDatabase("MobilityParameters") ){
         boost::shared_ptr<tbox::Database> mobility_db
            ( calphad_db->getDatabase( "MobilityParameters" ) );
         boost::shared_ptr<tbox::Database> species0_db
            ( mobility_db->getDatabase( "Species0" ) );
         boost::shared_ptr<tbox::Database> species1_db
            ( mobility_db->getDatabase( "Species1" ) );

         CALPHADMobility calphad_mobility0_phaseL("MobilitySpecies0");   
         calphad_mobility0_phaseL.initialize(species0_db->getDatabase( "PhaseL" ));
       
         CALPHADMobility calphad_mobility1_phaseL("MobilitySpecies1");
         calphad_mobility1_phaseL.initialize(species1_db->getDatabase( "PhaseL" ));
       
         CALPHADMobility calphad_mobility0_phaseA("MobilitySpecies0");
         calphad_mobility0_phaseA.initialize(species0_db->getDatabase( "PhaseA" ));
       
         CALPHADMobility calphad_mobility1_phaseA("MobilitySpecies1");
         calphad_mobility1_phaseA.initialize(species1_db->getDatabase( "PhaseA" ));
       
         const double tempmin=temperature*0.5;
         const double tempmax=temperature*2.;

         if( mpi.getRank()==0 )
         {
         ofstream tfile("D.dat", ios::out);
         tfile<<"#Diffusion in liquid phase for species 0 [m^2/s] vs. 10000./T"<<endl;
         calphad_mobility0_phaseL.printDiffusionVsTemperature(
            tempmin, tempmax,
            tfile );

         tfile<<endl<<"#Diffusion in liquid phase for species 1 [m^2/s] vs. 10000./T"<<endl;
         calphad_mobility0_phaseL.printDiffusionVsTemperature(
            tempmin, tempmax,
            tfile );

         tfile<<endl<<"#Diffusion in phase A for species 0 [m^2/s] vs. 10000./T"<<endl;
         calphad_mobility0_phaseA.printDiffusionVsTemperature(
            tempmin, tempmax,
            tfile );

         tfile<<endl<<"#Diffusion in phase A for species 1 [m^2/s] vs. 10000./T"<<endl;
         calphad_mobility1_phaseA.printDiffusionVsTemperature(
            tempmin, tempmax,
            tfile );
         }
      }
   }
}

//-----------------------------------------------------------------------

void QuatModel::postRunDiagnostics( void )
{
   d_integrator->printSolverTotals();
}

//-----------------------------------------------------------------------

void QuatModel::writeRestartFile( void )
{
   PFModel::writeRestartFile();
}

//=======================================================================

void QuatModel::printScalarDiagnostics( void )
{
   double total_energy, phase_energy, orient_energy, qint_energy;
   double well_energy, free_energy, eta_energy;

   evaluateEnergy(
      d_patch_hierarchy,
      d_time,
      total_energy,
      phase_energy,
      eta_energy,
      orient_energy,
      qint_energy,
      well_energy,
      free_energy,
      d_model_parameters.grand_potential() );

   if( d_model_parameters.with_heat_equation() ){
      double thermal_energy = computeThermalEnergy(d_patch_hierarchy);
      tbox::pout << "Thermal energy [pJ]= "<< thermal_energy << endl;
   }

   if ( ! d_time_info_interval->eventOccurredAtTime( d_time ) ) {
       tbox::pout << "cycle # " << d_cycle
                  << " : t = " << d_time << endl;
   }

   if ( d_extra_energy_detail ) {
      tbox::pout << setprecision(8);
      tbox::pout << "  Total energy [pJ]    = " << total_energy << endl;
      tbox::pout << "    phi energy [pJ]    = " << phase_energy << endl;
      tbox::pout << "    orient energy [pJ] = " << orient_energy << endl;
      tbox::pout << "    qint energy [pJ]   = " << qint_energy << endl;
      tbox::pout << "    well energy [pJ]   = " << well_energy << endl;
      tbox::pout << "    free energy [pJ]   = " << free_energy << endl;
      if( d_model_parameters.grand_potential() )
      {
         tbox::pout << "  GP [pJ] = " << total_energy << endl;
      }
      if ( d_model_parameters.with_third_phase() ) {
         tbox::pout << "    eta energy [pJ]    = " << eta_energy << endl;
      }
   }
   else {
      tbox::pout << "  Total energy [pJ] = " << total_energy << endl;
   }

   if ( d_temperature_strategy!=nullptr ) {
      double t = d_temperature_strategy->getCurrentMinTemperature( d_patch_hierarchy, d_time );
      tbox::pout << "  Min. Temperature = " << t << endl;
      t = d_temperature_strategy->getCurrentMaxTemperature( d_patch_hierarchy, d_time );
      tbox::pout << "  Max. Temperature = " << t << endl;
      t = d_temperature_strategy->getCurrentAverageTemperature( d_patch_hierarchy, d_time );
      tbox::pout << "  Average Temperature = " << t << endl;
   }

   // computes volume of physical domain
   const double* low = d_grid_geometry->getXLower();
   const double* up = d_grid_geometry->getXUpper();
   double vol = 1.;
   for(int d=0;d<NDIM;d++)
      vol *= (up[d]-low[d]);

   double vphi=vol;
   if ( d_model_parameters.with_phase() ){   
      vphi = evaluateVolumeSolid( d_patch_hierarchy );
      tbox::pout << "  Volume fraction of solid phase = " << vphi/vol << endl;    
   }

   if ( d_model_parameters.with_third_phase() ) {
      const double vphi_eta = evaluateVolumeEta( d_patch_hierarchy );
      tbox::pout << "  Volume fraction of eta phase = " << vphi_eta/vol << endl;    
   }

   if ( d_model_parameters.with_concentration() ){
      assert( d_work_id != -1 );

      math::HierarchyCellDataOpsReal<double> mathops( d_patch_hierarchy );

      for(int ic=0;ic<d_ncompositions;ic++){
         const double c0V0 = 
            evaluateIntegralConcentration(d_patch_hierarchy,ic);
         tbox::pout << "  Integral concentration "<<ic<<"= " << c0V0 << endl;

         copyDepthCellData(d_patch_hierarchy, d_work_id, 0,
                                              d_conc_id, ic);

         double cmax = mathops.max( d_work_id );
         tbox::pout << "  Max. concentration "<<ic<<"= " << cmax << endl;

         // average concentration
         const double c0 = c0V0 / vol;

         // now computes coring factor according to HBSM formula
         const double cphi = 
            evaluateIntegralPhaseConcentration( d_patch_hierarchy, ic );

         const double cex = ( cphi - c0*vphi ) / c0V0;
         tbox::pout << "  Cex (HBSM) for component "<<ic<<" = " << cex << endl;
      }
   }
}

//=======================================================================

void QuatModel::findAndNumberGrains( void )
{
   tbox::pout << "findAndNumberGrains" << endl;
   assert( d_grains );

   d_grains->findAndNumberGrains(d_patch_hierarchy,
      d_phase_id, d_weight_id, d_time);
}

//=======================================================================

void QuatModel::computeGrainDiagnostics( void )
{
   tbox::pout << "Computing grain diagnostics" << endl;

   findAndNumberGrains();
   d_grains->computeGrainVolumes(d_patch_hierarchy, d_weight_id);

   if ( d_model_parameters.with_concentration() ) {
      d_grains->computeGrainConcentrations(
         d_patch_hierarchy, d_time,
         d_conc_id, d_weight_id);
   }
}

//=======================================================================

void QuatModel::computeMinMaxQModulus(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy )
{
   if ( ! d_model_parameters.with_orientation() )return;
   if ( d_qlen == 1 ) return;
   
   double mx = 0.;
   double mn = 2.0;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++ ) {

      boost::shared_ptr<hier::PatchLevel > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
         boost::shared_ptr<hier::Patch > patch = *p;

         boost::shared_ptr< pdat::CellData<double> > y (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_quat_id) ) );
         const hier::Box& pbox = patch->getBox();

         boost::shared_ptr< pdat::CellData<double> > q_norm_err;
         if ( d_model_parameters.with_extra_visit_output() ) {
            q_norm_err = boost::dynamic_pointer_cast<pdat::CellData<double>,hier::PatchData>
                         ( patch->getPatchData(d_quat_norm_error_id) );
         }

         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i!= iend; ++i) {
            pdat::CellIndex cell = *i;
            double qnorm2 = 0.;
            for ( int q = 0 ; q < d_qlen ; q++ ) {
               qnorm2 += (*y)(cell,q) * (*y)(cell,q);
            }
            double qnorm = sqrt( qnorm2 );
            if ( qnorm > mx ) mx = qnorm;
            if ( qnorm < mn ) mn = qnorm;

            if ( d_model_parameters.with_extra_visit_output() ) {
               (*q_norm_err)(cell) = qnorm - 1.;
            }
         }
      }
   }

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   mpi.AllReduce(&mx,1,MPI_MAX);
   mpi.AllReduce(&mn,1,MPI_MIN);

   tbox::pout << endl;
   tbox::pout << "Q magnitude maximum - 1: " << mx - 1. << endl;
   tbox::pout << "Q magnitude minimum - 1: " << mn - 1. << endl;
}

//=======================================================================

void QuatModel::Regrid(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy )
{
   d_regrid_refine_alg.reset( new xfer::RefineAlgorithm() );

   d_regrid_refine_alg->registerRefine(
      d_phase_id,          // destination
      d_phase_id,          // source
      d_phase_scratch_id,  // temporary
      d_phase_refine_op );

   d_regrid_refine_alg->registerRefine(
      d_temperature_id,          // destination
      d_temperature_id,          // source
      d_temperature_scratch_id,  // temporary
      d_phase_refine_op );

   if ( d_model_parameters.with_third_phase() ) {
      d_regrid_refine_alg->registerRefine(
         d_eta_id,          // destination
         d_eta_id,          // source
         d_eta_scratch_id,  // temporary
         d_eta_refine_op );
   }

   if ( d_model_parameters.with_orientation() ){
      assert( d_quat_scratch_id >= 0 );
      d_regrid_refine_alg->registerRefine(
         d_quat_id,          // destination
         d_quat_id,          // source
         d_quat_scratch_id,  // temporary
         d_quat_refine_op );
      d_regrid_refine_alg->registerRefine(
         d_quat_relax_id,          // destination
         d_quat_relax_id,          // source
         d_quat_relax_scratch_id,  // temporary
         d_quat_refine_op );
   }

   if ( d_model_parameters.with_heat_equation() ) {
      assert( d_temperature_scratch_id >= 0 );
      d_regrid_refine_alg->registerRefine(
         d_temperature_id,          // destination
         d_temperature_id,          // source
         d_temperature_scratch_id,  // temporary
         d_phase_refine_op );
   }

   if ( d_model_parameters.with_concentration() ) {
      assert( d_conc_scratch_id >= 0 );
      d_regrid_refine_alg->registerRefine(
         d_conc_id,          // destination
         d_conc_id,          // source
         d_conc_scratch_id,  // temporary
         d_conc_refine_op );
   }

   if ( d_use_warm_start ) {

      set< int > cpodes_id_set;
      set< int > phase_id_set;
      set< int > eta_id_set;
      set< int > orient_id_set;
      set< int > conc_id_set;
      set< int > temp_id_set;

      d_integrator->getCPODESIdsRequiringRegrid(
         cpodes_id_set, phase_id_set, eta_id_set, orient_id_set, conc_id_set, temp_id_set );

      set< int >::iterator it;

      static map< int, int > id_map;

      static bool first_time = true;

      if ( first_time ) {
         hier::VariableDatabase* variable_db =
            hier::VariableDatabase::getDatabase();

         boost::shared_ptr< hier::VariableContext > scratch =
            variable_db->getContext( "SCRATCH" );

         for ( it = cpodes_id_set.begin(); it != cpodes_id_set.end(); it++ ) {

            std::ostringstream name;
            name << "warm_start_tmp_" << *it;

            boost::shared_ptr< pdat::CellVariable<double> > var
               (new pdat::CellVariable<double>(tbox::Dimension(NDIM), name.str() ));

            int id = variable_db->registerVariableAndContext(
               var, scratch, hier::IntVector( tbox::Dimension(NDIM),NGHOSTS ) );

            id_map[*it] = id;
         }

         first_time = false;
      }

      for ( it = phase_id_set.begin(); it != phase_id_set.end(); it++ ) {
         d_regrid_refine_alg->registerRefine(
            *it, *it, id_map[*it], d_phase_refine_op );
      }

      for ( it = eta_id_set.begin(); it != eta_id_set.end(); it++ ) {
         d_regrid_refine_alg->registerRefine(
            *it, *it, id_map[*it], d_phase_refine_op );
      }

      for ( it = orient_id_set.begin(); it != orient_id_set.end(); it++ ) {
         d_regrid_refine_alg->registerRefine(
            *it, *it, id_map[*it], d_quat_refine_op );
      }

      for ( it = conc_id_set.begin(); it != conc_id_set.end(); it++ ) {
         d_regrid_refine_alg->registerRefine(
            *it, *it, id_map[*it], d_conc_refine_op );
      }

      for ( it = temp_id_set.begin(); it != temp_id_set.end(); it++ ) {
         d_regrid_refine_alg->registerRefine(
            *it, *it, id_map[*it], d_phase_refine_op );
      }

   }  // if ( d_use_warm_start )

   PFModel::Regrid( hierarchy );
}

//=======================================================================
//
// Methods inherited from Serializable
//

void QuatModel::putToRestart(
   const boost::shared_ptr<tbox::Database>& db )const
{
   PFModel::putToRestart( db );
}

//=======================================================================
//
// Method inherited from StandardTagAndInitStrategy
//

void QuatModel::initializeLevelData(
   /*! Hierarchy to initialize */
   const boost::shared_ptr<hier::PatchHierarchy >& hierarchy,
   /*! Level to initialize */
   const int level_number,
   const double time,
   const bool can_be_refined,
   /*! Whether level is being introduced for the first time */
   const bool initial_time,
   /*! Level to copy data from */
   const boost::shared_ptr<hier::PatchLevel >& old_level,
   const bool allocate_data )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   if (old_level) {
      TBOX_ASSERT(level_number == old_level->getLevelNumber());
   }
   TBOX_ASSERT((hierarchy->getPatchLevel(level_number)));
#endif

   tbox::pout<<"QuatModel::initializeLevelData()"<<endl;
   
   assert( d_curr_to_curr_refine_alg );

   // Note that this method is pure virtual in PFModel and MUST be implemented here.
   // However, we might also have some default behavior in PFModel, so it should be called.
   PFModel::initializeLevelData(
      hierarchy, level_number, time, can_be_refined, initial_time, old_level,
      allocate_data );

   boost::shared_ptr<hier::PatchLevel > level (
      hierarchy->getPatchLevel( level_number ) );
   assert( level );

   AllocateLocalPatchData( level, time, allocate_data );

   d_integrator->initializeLevelData(
      hierarchy, level_number, time, can_be_refined, initial_time, old_level,
      allocate_data );
   if( d_model_parameters.with_orientation() )
      d_integrator_quat_only->initializeLevelData(
         hierarchy, level_number, time, can_be_refined, initial_time, old_level,
         allocate_data );

   d_grains->initializeLevelData(
      hierarchy, level_number, time, can_be_refined, initial_time, old_level,
      allocate_data );


   if ( initial_time ) {

      if ( level_number == d_level_of_init_data ) {

         d_initial_level = d_patch_hierarchy->getPatchLevel( d_level_of_init_data );
         assert( d_initial_level );

         boost::shared_ptr<xfer::RefineSchedule > schedule (
            d_curr_to_curr_refine_alg->createSchedule(
               level, d_initial_level, d_all_refine_patch_strategy ) );

         schedule->fillData( time, false );

      }
      else if ( level_number > d_level_of_init_data ) {

         boost::shared_ptr<xfer::RefineSchedule > schedule(
            d_curr_to_curr_refine_alg->createSchedule(
               level, old_level, level_number-1,
               hierarchy, d_all_refine_patch_strategy ) );

         schedule->fillData( time, false );
      }
      else {
         xfer::CoarsenAlgorithm coarsen_alg(tbox::Dimension(NDIM));

         coarsen_alg.registerCoarsen(
            d_phase_id, d_phase_id,
            d_grid_geometry->lookupCoarsenOperator(
               d_phase_var,
               "CONSERVATIVE_COARSEN") );

         if ( d_model_parameters.with_third_phase() ) {
            coarsen_alg.registerCoarsen(
               d_eta_id, d_eta_id,
               d_grid_geometry->lookupCoarsenOperator(
                  d_eta_var,
                  "CONSERVATIVE_COARSEN") );
         }

         if ( d_model_parameters.with_heat_equation() ){
            coarsen_alg.registerCoarsen(
               d_temperature_id, d_temperature_id,
               d_grid_geometry->lookupCoarsenOperator(
                  d_temperature_var,
                  "CONSERVATIVE_COARSEN") );
         }

         if ( d_model_parameters.with_concentration() ) {
            coarsen_alg.registerCoarsen(
               d_conc_id, d_conc_id,
               d_grid_geometry->lookupCoarsenOperator(
                  d_conc_var,
                  "CONSERVATIVE_COARSEN") );
         }

         if ( d_model_parameters.with_orientation() ){
            assert( d_quat_coarsen_op );
   
            coarsen_alg.registerCoarsen(
               d_quat_id, d_quat_id,
               d_quat_coarsen_op );
         }
         boost::shared_ptr<xfer::CoarsenSchedule > schedule =
            coarsen_alg.createSchedule(
               level, d_initial_level, nullptr );

         schedule->coarsenData();
      }

   }
   else {  // if ( initial_time )
      
      if ( ( level_number > 0 ) || old_level ) {
         // move solution from old_level to level
         d_regrid_refine_alg->createSchedule(
            level,
            old_level,
            level_number-1,
            hierarchy,
            d_all_refine_patch_strategy )
            ->fillData( time );
      }
   }

   if ( ! d_model_parameters.with_concentration() ){
      for ( hier::PatchLevel::Iterator ip(level->begin()); ip!=level->end(); ip++ ) {
         boost::shared_ptr<hier::Patch > patch = *ip;

         if ( d_model_parameters.free_energy_type()[0]=='s' ){
            boost::shared_ptr< pdat::CellData<double> > fl (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_f_l_id) ) );
            assert( fl );
            fl->fillAll( d_model_parameters.free_energy_liquid() );
         }
         boost::shared_ptr< pdat::CellData<double> > fa (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_f_a_id) ) );
         assert( fa );
         fa->fillAll( d_model_parameters.free_energy_solid_A() );

         if ( d_model_parameters.with_third_phase() ) {
            boost::shared_ptr< pdat::CellData<double> > fb (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_f_b_id) ) );
            assert( fb );
            fb->fillAll( d_model_parameters.free_energy_solid_B() );
         }
      }
   }
}

//=======================================================================

void QuatModel::AllocateQuatLocalPatchData(
   const boost::shared_ptr< hier::PatchLevel > level,
   const double time,
   const bool zero_data )
{
   assert( d_quat_id >= 0 );

   if ( !level->checkAllocated( d_quat_id ) )
      level->allocatePatchData( d_quat_id, time );

   if ( !level->checkAllocated( d_quat_relax_id ) )
      level->allocatePatchData( d_quat_relax_id, time );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_quat_scratch_id, level, time, zero_data );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_quat_relax_scratch_id, level, time, zero_data );

   AllocateAndZeroData< pdat::SideData<double> >(
      d_quat_diffusion_id, level, time, zero_data );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_quat_mobility_id, level, time, zero_data );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_quat_grad_cell_id, level, time, zero_data );

   AllocateAndZeroData< pdat::SideData<double> >(
      d_quat_grad_side_id, level, time, zero_data );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_quat_grad_modulus_id, level, time, zero_data );

   AllocateAndZeroData< pdat::SideData<double> >(
      d_quat_diffs_id, level, time, zero_data );

   if ( d_symmetry_aware ) {
      assert( d_quat_symm_rotation_id>=0 );
      AllocateAndZeroData< pdat::SideData<int> >(
         d_quat_symm_rotation_id, level, time, zero_data );
   }

   if ( d_model_parameters.with_extra_visit_output() ) {
      AllocateAndZeroData< pdat::CellData<double> >(
         d_quat_diffs_cell_id, level, time, zero_data );
      AllocateAndZeroData< pdat::CellData<double> >(
         d_quat_nonsymm_diffs_cell_id, level, time, zero_data );

      AllocateAndZeroData< pdat::CellData<double> >(
         d_quat_norm_error_id, level, time, zero_data );

      if ( d_symmetry_aware ) {
         assert( d_quat_symm_rotation_cell_id>=0 );
         AllocateAndZeroData< pdat::CellData<int> >(
            d_quat_symm_rotation_cell_id, level, time, zero_data );
      }
   }

   if ( d_model_parameters.with_visit_energy_output() ) {
      AllocateAndZeroData< pdat::CellData<double> >(
         d_energy_diag_id, level, time, zero_data );
   }
}
//=======================================================================

void QuatModel::AllocateLocalPatchData(
   const boost::shared_ptr< hier::PatchLevel > level,
   const double time,
   const bool zero_data )
{
   assert( d_temperature_scratch_id>=0 );
   
   if( d_phase_id>=0 ){
      if ( !level->checkAllocated( d_phase_id ) ) {
         level->allocatePatchData( d_phase_id, time );
      }

      AllocateAndZeroData< pdat::CellData<double> >(
         d_phase_scratch_id, level, time, zero_data );
   }

   if ( d_model_parameters.with_third_phase() ) {
      if ( !level->checkAllocated( d_eta_id ) ) {
         level->allocatePatchData( d_eta_id, time );
      }

      AllocateAndZeroData< pdat::CellData<double> >(
         d_eta_scratch_id, level, time, zero_data );
   }

   if ( !level->checkAllocated( d_temperature_id ) ) {
      level->allocatePatchData( d_temperature_id, time );
   }

   AllocateAndZeroData< pdat::CellData<double> >(
      d_temperature_scratch_id, level, time, zero_data );
   if( d_model_parameters.with_steady_temperature() ){
      assert( d_temperature_rhs_steady_id >= 0 );
      AllocateAndZeroData< pdat::CellData<double> >(
         d_temperature_rhs_steady_id, level, time, zero_data );
   }
   if( d_model_parameters.with_bias_well() ){
      assert( d_equilibrium_temperature_id >= 0 );
      AllocateAndZeroData< pdat::CellData<double> >(
         d_equilibrium_temperature_id, level, time, zero_data );
   }
   
   if( d_model_parameters.with_heat_equation() ){
      assert( d_cp_id >= 0 );
   
      AllocateAndZeroData< pdat::CellData<double> >(
         d_cp_id, level, time, zero_data );
   }
   
   if ( d_model_parameters.with_concentration() ){
      if ( !level->checkAllocated( d_conc_id ) ) {
         AllocateAndZeroData< pdat::CellData<double> >(
            d_conc_id, level, time, zero_data );
      }
      AllocateAndZeroData< pdat::CellData<double> >(
         d_conc_scratch_id, level, time, zero_data );
      for(int ic=0;ic<d_ncompositions;ic++){
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_pfm_diffusion_id[ic], level, time, zero_data );
      }
      if( d_model_parameters.concRHSstrategyIsKKS()
       || d_model_parameters.concRHSstrategyIsBeckermann() ){
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_phase_coupling_diffusion_id, level, time, zero_data );
      }else{
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_pfm_diffusion_l_id, level, time, zero_data );
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_pfm_diffusion_a_id, level, time, zero_data );
         if ( d_model_parameters.with_third_phase() ) {
            AllocateAndZeroData< pdat::SideData<double> >(
               d_conc_pfm_diffusion_b_id, level, time, zero_data );
         }
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_diffusion_coeff_l_id, level, time, zero_data );
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_diffusion_coeff_a_id, level, time, zero_data );
         if ( d_model_parameters.with_third_phase() ) {
            AllocateAndZeroData< pdat::SideData<double> >(
               d_conc_diffusion_coeff_b_id, level, time, zero_data );
         }
      }

      if ( d_model_parameters.concentrationModelNeedsPhaseConcentrations() ) {
         if ( !level->checkAllocated( d_conc_l_id ) ) {
            AllocateAndZeroData< pdat::CellData<double> >(
               d_conc_l_id, level, time, zero_data );
            AllocateAndZeroData< pdat::CellData<double> >(
               d_conc_a_id, level, time, zero_data );
            if ( d_model_parameters.with_third_phase() )
               AllocateAndZeroData< pdat::CellData<double> >(
                  d_conc_b_id, level, time, zero_data );
         }
         if ( d_model_parameters.isConcentrationModelCALPHAD() )
         if ( !level->checkAllocated( d_conc_l_ref_id ) ) {
            AllocateAndZeroData< pdat::CellData<double> >(
               d_conc_l_ref_id, level, time, zero_data );
            AllocateAndZeroData< pdat::CellData<double> >(
               d_conc_a_ref_id, level, time, zero_data );
            if ( d_model_parameters.with_third_phase() )
               AllocateAndZeroData< pdat::CellData<double> >(
                  d_conc_b_ref_id, level, time, zero_data );
         }
         AllocateAndZeroData< pdat::CellData<double> >(
            d_conc_l_scratch_id, level, time, zero_data );
         AllocateAndZeroData< pdat::CellData<double> >(
            d_conc_a_scratch_id, level, time, zero_data );
         if ( d_model_parameters.with_third_phase() )
            AllocateAndZeroData< pdat::CellData<double> >(
               d_conc_b_scratch_id, level, time, zero_data );
      }

      if ( d_model_parameters.with_third_phase() ) {
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_eta_coupling_diffusion_id, level, time, zero_data );
      }
      
      if ( d_model_parameters.with_gradT() ) {
         AllocateAndZeroData< pdat::SideData<double> >(
            d_conc_Mq_id, level, time, zero_data );
      }
      if ( d_model_parameters.with_partition_coeff() ){
         AllocateAndZeroData< pdat::CellData<double> >(
            d_partition_coeff_id, level, time, zero_data );
         if( d_model_parameters.needGhosts4PartitionCoeff() )
            AllocateAndZeroData< pdat::CellData<double> >(
               d_partition_coeff_scratch_id, level, time, zero_data );
      }
      
   } // d_model_parameters.with_concentration()
   
   if ( d_model_parameters.with_velocity() ){
      AllocateAndZeroData< pdat::CellData<double> >(
         d_velocity_id, level, time, zero_data );
   }

   if ( d_model_parameters.with_phase() ) {
      AllocateAndZeroData< pdat::SideData<double> >(
         d_phase_diffs_id, level, time, zero_data );
      if ( d_model_parameters.with_extra_visit_output() ) {
         AllocateAndZeroData< pdat::CellData<double> >(
            d_phase_diffs_cell_id, level, time, zero_data );
      }

      AllocateAndZeroData< pdat::CellData<double> >(
         d_phase_grad_cell_id, level, time, zero_data );

      AllocateAndZeroData< pdat::CellData<double> >(
         d_phase_mobility_id, level, time, zero_data );
   }

   if ( d_model_parameters.with_third_phase() ) {
      AllocateAndZeroData< pdat::SideData<double> >(
         d_eta_diffs_id, level, time, zero_data );

      AllocateAndZeroData< pdat::CellData<double> >(
         d_eta_grad_cell_id, level, time, zero_data );

      AllocateAndZeroData< pdat::CellData<double> >(
         d_eta_mobility_id, level, time, zero_data );
   }

   // AllocateAndZeroData< pdat::SideData<double> >(
   //    d_phase_grad_side_id, level, time, zero_data );

   if ( d_model_parameters.with_orientation() ) {
      AllocateQuatLocalPatchData( level, time, zero_data );
   }

   AllocateAndZeroData< pdat::CellData<double> >(
      d_weight_id, level, time, zero_data );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_work_id, level, time, zero_data );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_f_l_id, level, time, zero_data );

   AllocateAndZeroData< pdat::CellData<double> >(
      d_f_a_id, level, time, zero_data );

   if ( d_model_parameters.with_third_phase() ) {
      AllocateAndZeroData< pdat::CellData<double> >(
         d_f_b_id, level, time, zero_data );
   }
}
//-----------------------------------------------------------------------

void QuatModel::DeallocateIntermediateLocalPatchData(
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy )
{
   for (int ln0 = 0; ln0 <= d_patch_hierarchy->getFinestLevelNumber(); ln0++) {
      boost::shared_ptr<hier::PatchLevel > level (
         hierarchy->getPatchLevel( ln0 ) );

      DeallocateIntermediateLocalPatchData( level );
   }
}

void QuatModel::DeallocateIntermediateLocalPatchData(
   const boost::shared_ptr< hier::PatchLevel > level )
{
   if ( d_model_parameters.with_concentration() ){
      if ( d_model_parameters.isConcentrationModelCALPHAD() )
      if ( level->checkAllocated( d_conc_l_ref_id ) ) {
         level->deallocatePatchData( d_conc_l_ref_id );
         level->deallocatePatchData( d_conc_a_ref_id );
         if ( d_model_parameters.with_third_phase() )
            level->deallocatePatchData( d_conc_b_ref_id );
      }
   }
   
   if ( d_model_parameters.with_orientation() ) {
      level->deallocatePatchData( d_quat_grad_cell_id );
      level->deallocatePatchData( d_quat_grad_side_id );
      level->deallocatePatchData( d_quat_grad_modulus_id );
      level->deallocatePatchData( d_quat_diffs_id );
      level->deallocatePatchData( d_quat_relax_id );
      level->deallocatePatchData( d_quat_relax_scratch_id );
   
      if ( d_symmetry_aware ) {
         level->deallocatePatchData( d_quat_symm_rotation_id );
         if ( d_model_parameters.with_extra_visit_output() ) {
            level->deallocatePatchData( d_quat_symm_rotation_cell_id );
         }
      }
   }
}

//-----------------------------------------------------------------------

template <typename T>
void QuatModel::AllocateAndZeroData(
   const int data_id,
   const boost::shared_ptr< hier::PatchLevel > level,
   const double time,
   const bool zero_data )
{
   assert( data_id>=0 );
   
   if ( !level->checkAllocated( data_id ) )
      level->allocatePatchData( data_id, time );

   if ( zero_data ) {
      for ( hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p ) {
         boost::shared_ptr<T> data (
            BOOST_CAST< T, hier::PatchData>(p->getPatchData( data_id ) ) );
         data->fillAll( 0 );
      }
   }
}

//=======================================================================
//
// Method inherited from StandardTagAndInitStrategy

void QuatModel::resetHierarchyConfiguration(
   const boost::shared_ptr<hier::PatchHierarchy >& hierarchy,
   const int coarsest_level,
   const int finest_level )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (coarsest_level <= finest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));
   for (int ln0 = 0; ln0 <= finest_level; ln0++) {
      TBOX_ASSERT((hierarchy->getPatchLevel(ln0)));
   }
#endif
   tbox::pout<<"QuatModel::resetHierarchyConfiguration()"<<endl;
   
   int nlev = hierarchy->getNumberOfLevels();

   // Note that this method is pure virtual in PFModel and MUST be implemented here.
   // However, we might have some default behavior in PFModel, so it should be called.
   PFModel::resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level );

   /*
    * Check for overlapping boxes and abort if we find any.  For the time
    * being, we need to avoid this situation since the post-1.13 versions
    * of the Hypre SMG solver do not tolerate overlapping boxes.
    */
   for ( int ln = coarsest_level; ln <= finest_level; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );
      hier::BoxContainer bl( level->getBoxes() );
      if ( bl.boxesIntersect() ) {
         for ( hier::BoxContainer::iterator bli=bl.begin();
               bli != bl.end();
               ++bli ) {
            tbox::pout << *bli << endl;
         }
         stringstream message;
         message << "Boxes intersect on level " << ln;
         tbox::Utilities::abort(message.str(), __FILE__, __LINE__);
      }
   }


   d_integrator->resetHierarchyConfiguration(
      hierarchy, coarsest_level, finest_level );
   if( d_model_parameters.with_orientation() )
      d_integrator_quat_only->resetHierarchyConfiguration(
         hierarchy, coarsest_level, finest_level );

   d_curr_to_scr_refine_sched.resize( nlev );
   
   for ( int ln = coarsest_level; ln <= finest_level; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );

      boost::shared_ptr<xfer::RefineSchedule > schedule (
         d_curr_to_scr_refine_alg->createSchedule(
            level,
            ln-1,
            hierarchy,
            d_all_refine_patch_strategy ) );
      d_curr_to_scr_refine_sched[ln] = schedule ;
   }

   if ( d_grain_diag_interval->isActive() 
     || d_grain_extend_interval->isActive() ) {
      d_grains->resetHierarchyConfiguration(
         hierarchy, coarsest_level, finest_level);
      
   }

   computeVectorWeights( d_patch_hierarchy, -1, -1 );
}

//=======================================================================
//
// Method inherited from StandardTagAndInitStrategy

void QuatModel::applyGradientDetector(
   boost::shared_ptr<hier::PatchHierarchy >& hierarchy,
   int level_number,
   double time,
   int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too )
{
   boost::shared_ptr<hier::PatchLevel >
      level = hierarchy->getPatchLevel(level_number);

   copyCurrentToScratch(
      hierarchy,
      level_number,
      time,
      d_all_refine_patch_strategy );

   if ( d_tag_phase ) {
      computePhaseDiffs(
         level,
         d_phase_scratch_id,
         d_phase_diffs_id,
         time );

      computePhaseGradCell(
         level,
         d_phase_diffs_id,
         d_phase_grad_cell_id,
         time );
   }

   if ( d_tag_eta ) {
      computeEtaDiffs(
         level,
         d_eta_scratch_id,
         d_eta_diffs_id,
         time );

      computeEtaGradCell(
         level,
         d_eta_diffs_id,
         d_eta_grad_cell_id,
         time );
   }

   if ( d_tag_quat ) {
      computeQuatDiffs(
         level,
         d_quat_scratch_id,
         d_quat_diffs_id,
         time );

      computeQuatGradCell(
         level,
         d_quat_diffs_id,
         d_quat_grad_cell_id,
         time );
   }

   for (hier::PatchLevel::Iterator ip(level->begin()); ip!=level->end(); ++ip) {
      boost::shared_ptr<hier::Patch > patch = *ip;

      tagGradientDetectorCells(
         *patch,
         time,
         initial_time,
         tag_index,
         uses_richardson_extrapolation_too );
   }
}

void QuatModel::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_richardson_extrapolation_too)
{
   (void) regrid_time;
   (void) initial_error;
   (void) uses_richardson_extrapolation_too;

   boost::shared_ptr< pdat::CellData<int> > tags(
      BOOST_CAST< pdat::CellData<int>, hier::PatchData>( patch.getPatchData( tag_index) ) );
   assert( tags );
   tags->fillAll( 0 );

   const hier::Index& patch_lower = patch.getBox().lower();
   const hier::Index& patch_upper = patch.getBox().upper();

   const hier::Box& tag_ghost_box = tags->getGhostBox();
   const hier::Index& tag_gbox_lower = tag_ghost_box.lower();
   const hier::Index& tag_gbox_upper = tag_ghost_box.upper();

   boost::shared_ptr< geom::CartesianPatchGeometry > patch_geom ( 
      BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(patch.getPatchGeometry()) );
   TBOX_ASSERT(patch_geom);

   const double * dx = patch_geom->getDx();

   if ( d_tag_phase ) {
      boost::shared_ptr< pdat::CellData<double> > phase_grad (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( d_phase_grad_cell_id) ) );
      assert( phase_grad );

      const hier::Box& gbox = phase_grad->getGhostBox();
      const hier::Index& gbox_lower = gbox.lower();
      const hier::Index& gbox_upper = gbox.upper();

      FORT_TAG_GRAD(
         patch_lower(0), patch_upper(0),
         patch_lower(1), patch_upper(1),
#if (NDIM==3)
         patch_lower(2), patch_upper(2),
#endif
         phase_grad->getPointer(0),
         phase_grad->getPointer(1),
#if (NDIM==3)
         phase_grad->getPointer(2),
#endif
         gbox_lower(0), gbox_upper(0),
         gbox_lower(1), gbox_upper(1),
#if (NDIM==3)
         gbox_lower(2), gbox_upper(2),
#endif
         dx,
         tags->getPointer(0),
         tag_gbox_lower(0), tag_gbox_upper(0),
         tag_gbox_lower(1), tag_gbox_upper(1),
#if (NDIM==3)
         tag_gbox_lower(2), tag_gbox_upper(2),
#endif
         true,
         d_phase_threshold_untagged );
   }

   if ( d_tag_eta ) {
      boost::shared_ptr< pdat::CellData<double> > eta_grad (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( d_eta_grad_cell_id) ) );
      assert( eta_grad );

      const hier::Box& gbox = eta_grad->getGhostBox();
      const hier::Index& gbox_lower = gbox.lower();
      const hier::Index& gbox_upper = gbox.upper();

      FORT_TAG_GRAD(
         patch_lower(0), patch_upper(0),
         patch_lower(1), patch_upper(1),
#if (NDIM==3)
         patch_lower(2), patch_upper(2),
#endif
         eta_grad->getPointer(0),
         eta_grad->getPointer(1),
#if (NDIM==3)
         eta_grad->getPointer(2),
#endif
         gbox_lower(0), gbox_upper(0),
         gbox_lower(1), gbox_upper(1),
#if (NDIM==3)
         gbox_lower(2), gbox_upper(2),
#endif
         dx,
         tags->getPointer(0),
         tag_gbox_lower(0), tag_gbox_upper(0),
         tag_gbox_lower(1), tag_gbox_upper(1),
#if (NDIM==3)
         tag_gbox_lower(2), tag_gbox_upper(2),
#endif
         true,
         d_eta_threshold_untagged );
   }

   if ( d_tag_quat ) {

      boost::shared_ptr< pdat::CellData<double> > quat_grad (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch.getPatchData( d_quat_grad_cell_id) ) );
      assert( quat_grad );

      const hier::Box& gbox = quat_grad->getGhostBox();
      const hier::Index& gbox_lower = gbox.lower();
      const hier::Index& gbox_upper = gbox.upper();

      FORT_TAG_QUAT_GRAD(
         patch_lower(0), patch_upper(0),
         patch_lower(1), patch_upper(1),
#if (NDIM==3)
         patch_lower(2), patch_upper(2),
#endif
         quat_grad->getPointer( 0 * d_qlen ),
         quat_grad->getPointer( 1 * d_qlen ),
#if (NDIM==3)
         quat_grad->getPointer( 2 * d_qlen ),
#endif
         gbox_lower(0), gbox_upper(0),
         gbox_lower(1), gbox_upper(1),
#if (NDIM==3)
         gbox_lower(2), gbox_upper(2),
#endif
         d_qlen,
         dx,
         tags->getPointer(0),
         tag_gbox_lower(0), tag_gbox_upper(0),
         tag_gbox_lower(1), tag_gbox_upper(1),
#if (NDIM==3)
         tag_gbox_lower(2), tag_gbox_upper(2),
#endif
         true,
         d_quat_threshold_untagged );
   }
}
//=======================================================================

// Read initialization database

void QuatModel::readInitialDatabase(
   boost::shared_ptr<tbox::Database> input_db )
{
   PFModel::readInitialDatabase( input_db );

   hier::BoxContainer boxes = d_grid_geometry->getPhysicalDomain();
   assert( boxes.size() == 1 );
}

//=======================================================================

void QuatModel::WriteInitialConditionsFile( void )
{
   tbox::plog<<"Write initial conditions file..."<<endl;

   // get new PatchLevel with uniform mesh at level "d_initial_conditions_level"
   boost::shared_ptr<hier::PatchLevel > flattened_level = 
      FlattenHierarchy(
         d_patch_hierarchy,
         d_initial_conditions_level,
         d_time );

   // get size of uniform mesh to write
   hier::BoxContainer boxes = d_grid_geometry->getPhysicalDomain();
   assert( boxes.size() == 1 );

   hier::Box bf (boxes.front());
   bf.refine( flattened_level->getRatioToLevelZero() );

   int nx_prob = bf.numberCells(0);
   int ny_prob = bf.numberCells(1);
   int nz_prob = 1;
#if (NDIM > 2)   
   nz_prob = bf.numberCells(2);
#endif

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   const int npp = mpi.getSize();

   for ( int pp = 0; pp < npp; pp++ ) {

      //tbox::plog<<"pp="<<pp<<endl;
      if ( mpi.getRank() == pp ) {

#ifdef HAVE_NETCDF3
         NcFile* f;
         NcVar* nc_phase;
         NcVar* nc_eta=nullptr;
         NcVar** nc_conc =new NcVar*[d_ncompositions];
         NcVar** nc_qcomp=new NcVar*[d_qlen];
         NcVar* nc_temp=nullptr;
#endif
#ifdef HAVE_NETCDF4
         NcFile* f;
         NcVar nc_phase;
         NcVar nc_eta;
         NcVar* nc_conc =new NcVar[d_ncompositions];
         NcVar* nc_qcomp=new NcVar[d_qlen];
         NcVar nc_temp;
#endif

         if ( pp == 0 ) {
#ifdef HAVE_NETCDF3
            f = new NcFile( d_initial_conditions_file_name.c_str(), NcFile::Replace );
            if ( ! f->is_valid() ) {
               TBOX_ERROR("Cannot open file " << d_initial_conditions_file_name << endl);
            }
#endif
#ifdef HAVE_NETCDF4
            f=new NcFile( d_initial_conditions_file_name, NcFile::replace );
            if( f->isNull()) {
               TBOX_ERROR("Cannot open file " << d_initial_conditions_file_name << endl);
            }else{
               clog<<"Open/replace file "<<d_initial_conditions_file_name<<endl;
            }
#endif
 
#ifdef HAVE_NETCDF3
            NcDim* nc_nx = f->add_dim( "x", nx_prob );
            NcDim* nc_ny = f->add_dim( "y", ny_prob );
            NcDim* nc_nz = f->add_dim( "z", nz_prob );
            f->add_dim( "qlen", d_qlen );

            nc_phase = f->add_var( "phase", ncFloat, nc_nz, nc_ny, nc_nx );

            if ( d_model_parameters.with_third_phase() ) {
               nc_eta = f->add_var( "eta", ncFloat, nc_nz, nc_ny, nc_nx );
            }

            if ( d_model_parameters.with_orientation() ) {
               for ( int ii = 0; ii < d_qlen; ii++ ) {
                  std::ostringstream o;
                  o << "quat" << ii+1;
                  nc_qcomp[ii] =
                     f->add_var( o.str().c_str(), ncFloat, nc_nz, nc_ny, nc_nx );
               }
            }

            if ( d_model_parameters.with_concentration() ) {
               for ( int ii = 0; ii < d_ncompositions; ii++ ) {
                  std::ostringstream o;
                  o << "concentration";
                  if( d_ncompositions>1 )o << ii+1;
                  nc_conc[ii] =
                     f->add_var( o.str().c_str(), ncFloat, nc_nz, nc_ny, nc_nx );
               }
            }

            nc_temp = f->add_var( "temperature", ncFloat, nc_nz, nc_ny, nc_nx );

#endif
#ifdef HAVE_NETCDF4
            //cout<<"add variables from PE 0..."<<endl;
            NcDim nc_nx = f->addDim( "x", nx_prob );
            NcDim nc_ny = f->addDim( "y", ny_prob );
            NcDim nc_nz = f->addDim( "z", nz_prob );
            //f->addDim( "qlen", d_qlen );

            vector<NcDim> dims;
            dims.push_back(nc_nz);
            dims.push_back(nc_ny);
            dims.push_back(nc_nx);
            nc_phase = f->addVar("phase", ncFloat, dims);
            if( nc_phase.isNull() ){
               TBOX_ERROR( "Could add variable 'phase'" << endl );
            }
            if ( d_model_parameters.with_third_phase() ) {
               nc_eta = f->addVar( "eta", ncFloat, dims);
            }

            if ( d_model_parameters.with_orientation() ) {
               for ( int ii = 0; ii < d_qlen; ii++ ) {
                  std::ostringstream o;
                  o << "quat" << ii+1;
                  nc_qcomp[ii] =
                     f->addVar( o.str(), ncFloat, dims);
               }
            }

            if ( d_model_parameters.with_concentration() ) {
               for ( int ii = 0; ii < d_ncompositions; ii++ ) {
                  std::ostringstream o;
                  o << "concentration";
                  if( d_ncompositions>1 )o << ii+1;
                  nc_conc[ii] =
                     f->addVar( o.str(), ncFloat, dims);
               }
            }

            nc_temp = f->addVar( "temperature", ncFloat, dims);
            //cout<<"variables added on PE 0..."<<endl;
#endif
         }
         else { // pp!=0
#ifdef HAVE_NETCDF3
            f = new NcFile( d_initial_conditions_file_name.c_str(), NcFile::Write );
            if ( ! f->is_valid() ) {
               TBOX_ERROR("Cannot open file " << d_initial_conditions_file_name << endl);
            }
#endif
#ifdef HAVE_NETCDF4
            f=new NcFile( d_initial_conditions_file_name, NcFile::write );
            if ( f->isNull() ) {
               TBOX_ERROR("Cannot open file " << d_initial_conditions_file_name << endl);
            //}else{
            //   clog<<"Open/write file "<<d_initial_conditions_file_name<<endl;
            }
#endif

#ifdef HAVE_NETCDF3
            nc_phase = f->get_var( "phase" );

            if ( d_model_parameters.with_third_phase() ) {
               nc_eta = f->get_var( "eta" );
            }

            if ( d_model_parameters.with_orientation() ) {
               for ( int ii = 0; ii < d_qlen; ii++ ) {
                  std::ostringstream o;
                  o << "quat" << ii+1;
                  nc_qcomp[ii] = f->get_var( o.str().c_str() );
               }
            }

            if ( d_model_parameters.with_concentration() ) {
               for ( int ii = 0; ii < d_ncompositions; ii++ ) {
                  std::ostringstream o;
                  o << "concentration";
                  if( d_ncompositions>1 )o << ii+1;
                  nc_conc[ii] = f->get_var( o.str().c_str() );
               }
            }

            nc_temp = f->get_var( "temperature" );
#endif
#ifdef HAVE_NETCDF4
            //clog<<"add variables from PE >0..."<<endl;
            nc_phase = f->getVar( "phase" );

            if ( d_model_parameters.with_third_phase() ) {
               nc_eta = f->getVar( "eta" );
            }

            if ( d_model_parameters.with_orientation() ) {
               for ( int ii = 0; ii < d_qlen; ii++ ) {
                  std::ostringstream o;
                  o << "quat" << ii+1;
                  nc_qcomp[ii] = f->getVar( o.str() );
               }
            }

            if ( d_model_parameters.with_concentration() ) {
               for ( int ii = 0; ii < d_ncompositions; ii++ ) {
                  std::ostringstream o;
                  o << "concentration";
                  if( d_ncompositions>1 )o << ii+1;
                  nc_conc[ii] = f->getVar( o.str() );
               }
            }

            nc_temp = f->getVar( "temperature" );
#endif
         } // pp==0 or not        

#ifdef HAVE_NETCDF3
         if ( nc_phase == nullptr ) {
            TBOX_ERROR("Could not create variable 'phase'" << endl);
         }

         if ( d_model_parameters.with_third_phase() && nc_eta == nullptr ) {
            TBOX_ERROR("Could not create variable 'eta'" << endl);
         }

         if ( d_model_parameters.with_orientation() ) {
            for ( int ii = 0; ii < d_qlen; ii++ ) {
               std::ostringstream o;
               o << "quat" << ii+1;
               if ( nc_qcomp[ii] == nullptr ) {
                  TBOX_ERROR( "Could not create variable "<< o.str() << endl );
               }
            }
         }

         if ( d_model_parameters.with_concentration() ) {
            for ( int ii = 0; ii < d_ncompositions; ii++ ) {
               std::ostringstream o;
               o << "concentration";
               if( d_ncompositions>1 )o << ii+1;
               if ( nc_conc[ii] == nullptr ) {
                  TBOX_ERROR( "Could not create variable "<< o.str() << endl );
               }
            }
         }

         if ( nc_temp == nullptr ) {
            TBOX_ERROR( "Could not create variable 'temperature'" << endl );
         }
#endif

#ifdef HAVE_NETCDF4
         if ( nc_phase.isNull() ){
            TBOX_ERROR("Could not create variable 'phase'" << endl);
         }

         if ( d_model_parameters.with_third_phase() && nc_eta.isNull() ){
            TBOX_ERROR("Could not create variable 'eta'" << endl);
         }

         if ( d_model_parameters.with_orientation() ) {
            for ( int ii = 0; ii < d_qlen; ii++ ) {
               std::ostringstream o;
               o << "quat" << ii+1;
               if ( nc_qcomp[ii].isNull() ) {
                  TBOX_ERROR( "Could not create variable "<< o.str() << endl );
               }
            }
         }

         if ( d_model_parameters.with_concentration() ) {
            for ( int ii = 0; ii < d_ncompositions; ii++ ) {
               std::ostringstream o;
               o << "concentration";
               if( d_ncompositions>1 )o << ii+1;
               if ( nc_conc[ii].isNull() ) {
                  TBOX_ERROR( "Could not create variable "<< o.str() << endl );
               }
            }
         }

         if ( nc_temp.isNull() ) {
            TBOX_ERROR( "Could not create variable 'temperature'" << endl );
         }
#endif

         //cout<<"Write data into variable objects..."<<endl;
         for ( hier::PatchLevel::Iterator p(flattened_level->begin()); 
                                          p != flattened_level->end(); 
                                        ++p ) {
            boost::shared_ptr<hier::Patch > patch = *p;

            boost::shared_ptr< pdat::CellData<double> > phase_data (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_phase_id) ) );
            assert( phase_data );

            const hier::Box& this_b = patch->getBox();
            int nx = this_b.numberCells( 0 );
            int ny = this_b.numberCells( 1 );
            int nz = 1;
#if (NDIM > 2)
            nz = this_b.numberCells( 2 );
#endif
            int lowz=0;
#if (NDIM > 2)
            lowz=this_b.lower(2);
#endif

#ifdef HAVE_NETCDF3
            nc_phase->set_cur( lowz, this_b.lower(1), this_b.lower(0) );
            nc_phase->put( phase_data->getPointer(), nz, ny, nx );
#endif
#ifdef HAVE_NETCDF4
            std::vector<size_t> startp(3);
            startp[0]=lowz;
            startp[1]=this_b.lower(1);
            startp[2]=this_b.lower(0);
            std::vector<size_t> countp(3);
            countp[0]=nz;
            countp[1]=ny;
            countp[2]=nx;

            //cout<<"Write data into variable 'phase'"<<endl;
            //cout<<"nx="<<countp[0]<<", ny="<<countp[1]<<", nz="<<countp[2]<<endl;
            nc_phase.putVar(startp, countp, phase_data->getPointer());
            //cout<<"Data written into variable 'phase'"<<endl;
#endif
            if ( d_model_parameters.with_third_phase() ) {
               boost::shared_ptr< pdat::CellData<double> > eta_data (
                  BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_eta_id) ) );
               assert( eta_data );

#ifdef HAVE_NETCDF3
               nc_eta->set_cur( lowz, this_b.lower(1), this_b.lower(0) );
               nc_eta->put( eta_data->getPointer(), nz, ny, nx );
#endif
#ifdef HAVE_NETCDF4
               nc_eta.putVar(startp, countp, eta_data->getPointer());
#endif
            }

            boost::shared_ptr< pdat::CellData<double> > temp_data (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_temperature_id) ) );
            assert( temp_data );

#ifdef HAVE_NETCDF3
            nc_temp->set_cur( lowz, this_b.lower(1), this_b.lower(0) );
            nc_temp->put( temp_data->getPointer(), nz, ny, nx );
#endif
#ifdef HAVE_NETCDF4
            //cout<<"Write data into variable 'temperature'"<<endl;
            nc_temp.putVar(startp, countp, temp_data->getPointer());
#endif

            if ( d_model_parameters.with_orientation() ){
               boost::shared_ptr< pdat::CellData<double> > quat_data (
                  BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_quat_id) ) );
               assert( quat_data );

               for ( int dd = 0; dd < d_qlen; dd++ ) {
#ifdef HAVE_NETCDF3
                  nc_qcomp[dd]->set_cur( lowz, this_b.lower(1), this_b.lower(0) );
                  nc_qcomp[dd]->put( quat_data->getPointer( dd ), nz, ny, nx );
#endif
#ifdef HAVE_NETCDF4
                  nc_qcomp[dd].putVar(startp, countp, quat_data->getPointer( dd ));
#endif               
               }
            }

            if ( d_model_parameters.with_concentration() ){
               boost::shared_ptr< pdat::CellData<double> > conc_data (
                  BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
                     patch->getPatchData( d_conc_id) ) );
               assert( conc_data );

               for ( int dd = 0; dd < d_ncompositions; dd++ ) {
#ifdef HAVE_NETCDF3
                  nc_conc[dd]->set_cur( lowz, this_b.lower(1), this_b.lower(0) );
                  nc_conc[dd]->put( conc_data->getPointer( dd ), nz, ny, nx );
#endif
#ifdef HAVE_NETCDF4
                  nc_conc[dd].putVar(startp, countp, conc_data->getPointer( dd ));
#endif
               }
            }
            
         }
         //cout<<"Close file..."<<endl;
#ifdef HAVE_NETCDF3
         f->close();
#endif
#ifdef HAVE_NETCDF4
         delete f;
#endif
         delete[] nc_qcomp;
         delete[] nc_conc;
      }
      mpi.Barrier();
   }  // for ( int pp ...
   flattened_level.reset();
}

//=======================================================================

boost::shared_ptr<hier::PatchLevel >
QuatModel::FlattenHierarchy(
   const boost::shared_ptr< hier::PatchHierarchy > src_hierarchy,
   const int level_number,
   const double time )
{
   assert( level_number>=0 );
   assert( level_number<10 );

   // It will be assumed that the levels of src_hierarchy are
   // synchronized at this point.

   hier::IntVector ratio(tbox::Dimension(NDIM),1);
   for (int l = level_number; l > 0; l--)
      ratio *= d_patch_hierarchy->getRatioToCoarserLevel(l);
   d_ratio_of_init_to_coarsest = ratio;

   // Compute physical domain box array describing the index space of the physical domain managed 
   // by this geometry object. 
   // If any entry of ratio vector is negative, the index space is coarsened with respect to the 
   // physical domain description. Otherwise, the index space is refined.    
   
   // get boxes corresponding to level 0 of hierarchy, then refine them to "level_number"
   hier::BoxLevel layer0(ratio, d_grid_geometry);
   hier::BoxContainer boxes;
   boost::shared_ptr<hier::PatchLevel> zero_level =
      src_hierarchy->getPatchLevel( 0 );
   //iterate over patches
   for ( hier::PatchLevel::Iterator p(zero_level->begin()); p != zero_level->end(); ++p ) {
      const hier::Box& box = (*p)->getBox();
      boxes.pushBack(box);
   }
   
   //hier::BoxContainer boxes ( zero_level->getBoxes() );
   boxes.refine(ratio);   
   //boxes.print(cout);
   
   hier::BoxContainer::const_iterator boxes_itr=boxes.begin();
   for (int ib=0; ib < boxes.size(); ib++, boxes_itr++) {
      layer0.addBox(hier::Box(*boxes_itr, hier::LocalId(ib), layer0.getMPI().getRank()));
   }

   boost::shared_ptr<hier::PatchLevel> src_level =
      src_hierarchy->getPatchLevel( level_number );
   boost::shared_ptr<hier::PatchLevel > flattened_level (
      new hier::PatchLevel( layer0, d_grid_geometry, 
                            hier::VariableDatabase::getDatabase()->getPatchDescriptor()) );

   // allocate data on newly created level
   flattened_level->allocatePatchData( d_phase_id );
   flattened_level->allocatePatchData( d_phase_scratch_id );
   flattened_level->setTime( time, d_phase_id );
   flattened_level->setTime( time, d_phase_scratch_id );

   if ( d_model_parameters.with_third_phase() ) {
      flattened_level->allocatePatchData( d_eta_id );
      flattened_level->allocatePatchData( d_eta_scratch_id );
      flattened_level->setTime( time, d_eta_id );
      flattened_level->setTime( time, d_eta_scratch_id );
   }

   flattened_level->allocatePatchData( d_temperature_id );
   flattened_level->allocatePatchData( d_temperature_scratch_id );
   flattened_level->setTime( time, d_temperature_id );
   flattened_level->setTime( time, d_temperature_scratch_id );

   if ( d_model_parameters.with_orientation() ) {
      flattened_level->allocatePatchData( d_quat_id );
      flattened_level->allocatePatchData( d_quat_scratch_id );
      flattened_level->setTime( time, d_quat_id );
      flattened_level->setTime( time, d_quat_scratch_id );
   }
   if ( d_model_parameters.with_concentration() ) {
      flattened_level->allocatePatchData( d_conc_id );
      flattened_level->allocatePatchData( d_conc_scratch_id );
      flattened_level->setTime( time, d_conc_id );
      flattened_level->setTime( time, d_conc_scratch_id );
   }

   // fill new level with data
   boost::shared_ptr<xfer::RefineSchedule > schedule (
      d_curr_to_curr_refine_alg->createSchedule(
         flattened_level,
         src_level,
         level_number - 1,
         src_hierarchy,
         d_all_refine_patch_strategy ));

   schedule->fillData( time );

   schedule.reset();
   return flattened_level;
}

//=======================================================================

// phase_id is CellData with NGHOSTS
// phase_diffs_id is SideData with NGHOSTS

void QuatModel::computePhaseDiffs(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& phase_id,
   int& phase_diffs_id,
   const double time,
   const CACHE_TYPE cache )
{
   t_phase_diffs_timer->start();

   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computePhaseDiffs( patch_level, phase_id, phase_diffs_id, time );
   }

   t_phase_diffs_timer->stop();
}

void QuatModel::computePhaseDiffs(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& phase_id,
   int& phase_diffs_id,
   const double time )
{
   if ( phase_id < 0 ) phase_id = d_phase_scratch_id;
   if ( phase_diffs_id < 0 ) phase_diffs_id = d_phase_diffs_id;

   computeVarDiffs( level, phase_id, phase_diffs_id, time );

   if ( d_model_parameters.with_extra_visit_output() ) {
      
      for ( hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p ) {
         boost::shared_ptr<hier::Patch > patch = *p;
         const hier::Box& box = patch->getBox();

         boost::shared_ptr< pdat::SideData<double> > diff_data (
            BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( phase_diffs_id) ) );
         assert( diff_data );
         assert( diff_data->getGhostCellWidth() ==
                 hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

         boost::shared_ptr< pdat::CellData<double> > cell_diffs_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_phase_diffs_cell_id) ) );
         assert( cell_diffs_data );

         pdat::CellIterator iend(pdat::CellGeometry::end(box));
         for ( pdat::CellIterator i(pdat::CellGeometry::begin(box)); i!=iend; ++i ) {
            const pdat::CellIndex ccell = *i;
            const pdat::SideIndex xside(
               ccell, pdat::SideIndex::X, pdat::SideIndex::Lower );
            const pdat::SideIndex yside(
               ccell, pdat::SideIndex::Y, pdat::SideIndex::Lower );
#if (NDIM == 3)
            const pdat::SideIndex zside(
               ccell, pdat::SideIndex::Z, pdat::SideIndex::Lower );
#endif
            (*cell_diffs_data)(ccell,0) = (*diff_data)(xside);
            (*cell_diffs_data)(ccell,1) = (*diff_data)(yside);
#if (NDIM == 3)
            (*cell_diffs_data)(ccell,2) = (*diff_data)(zside);
#endif
         }
      }
   }
}

//=======================================================================

// eta_id is CellData with NGHOSTS
// eta_diffs_id is SideData with NGHOSTS

void QuatModel::computeEtaDiffs(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& eta_id,
   int& eta_diffs_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeEtaDiffs( patch_level, eta_id, eta_diffs_id, time );
   }
}

void QuatModel::computeEtaDiffs(
   const boost::shared_ptr< hier::PatchLevel > patch_level,
   int& eta_id,
   int& eta_diffs_id,
   const double time )
{
   if ( eta_id < 0 ) eta_id = d_eta_scratch_id;
   if ( eta_diffs_id < 0 ) eta_diffs_id = d_eta_diffs_id;

   computeVarDiffs( patch_level, eta_id, eta_diffs_id, time );
}

//=======================================================================

// var_id is CellData with NGHOSTS
// diff_id is SideData with NGHOSTS

void QuatModel::computeVarDiffs(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& var_id,
   int& diffs_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeVarDiffs( patch_level, var_id, diffs_id, time );
   }
}

void QuatModel::smoothQuat(
   const boost::shared_ptr< hier::PatchLevel > level )
{
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& box = patch->getBox();
      const hier::Index& ifirst = box.lower();
      const hier::Index& ilast = box.upper();

      boost::shared_ptr< pdat::CellData<double> > quat (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_quat_id) ) );
      assert( quat );
      assert( quat->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      boost::shared_ptr< pdat::CellData<double> > quat_scratch (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_quat_scratch_id) ) );
      assert( quat_scratch );
      assert( quat_scratch->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::SideData<double> > diff_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( d_quat_diffs_id) ) );
      assert( diff_data );
      assert( diff_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::CellData<double> > phase (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_phase_id) ) );
      assert( phase );
      assert( phase->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      FORT_SMOOTHQUAT(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         d_qlen,
         quat_scratch->getPointer(), NGHOSTS,
         quat->getPointer(),         0,
         diff_data->getPointer( 0, 0 ),
         diff_data->getPointer( 1, 0 ),
#if (NDIM == 3)
         diff_data->getPointer( 2, 0 ),
#endif
         diff_data->getGhostCellWidth()[0],
         phase->getPointer(), 0,
         d_phase_threshold );
   }
}

void QuatModel::smoothQuat(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const double time )
{
   //tbox::pout<<"smoothQuat..."<<endl;
      
   // Fill ghosts of original quat data
   copyCurrentToScratch(
      hierarchy,
      time,
      d_all_refine_patch_strategy );
   
   int maxln = hierarchy->getFinestLevelNumber();
   for ( int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatDiffs(
         patch_level,
         d_quat_scratch_id,
         d_quat_diffs_id,
         time );
      
      smoothQuat(patch_level);
   }
}

void QuatModel::computeVarDiffs(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& var_id,
   int& diffs_id,
   const double time )
{
   (void)time;

   if ( var_id < 0 ) var_id = d_phase_scratch_id;
   if ( diffs_id < 0 ) diffs_id = d_phase_diffs_id;
   assert( var_id>=0 );
   assert( diffs_id>=0 );

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& box = patch->getBox();
      const hier::Index& ifirst = box.lower();
      const hier::Index& ilast = box.upper();

      boost::shared_ptr< pdat::CellData<double> > var_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( var_id) ) );
      assert( var_data );
      assert( var_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::SideData<double> > diff_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( diffs_id) ) );
      assert( diff_data );
      assert( diff_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      const hier::Box & var_gbox = var_data->getGhostBox();
      const hier::Index& qlower = var_gbox.lower();
      const hier::Index& qupper = var_gbox.upper();

      const hier::Box & diff_gbox = diff_data->getGhostBox();
      const hier::Index& dlower = diff_gbox.lower();
      const hier::Index& dupper = diff_gbox.upper();

      FORT_DIFFS(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         var_data->getPointer(),
         qlower[0], qupper[0],
         qlower[1], qupper[1],
#if (NDIM == 3)
         qlower[2], qupper[2],
#endif
         diff_data->getPointer(0),
         diff_data->getPointer(1),
#if (NDIM == 3)
         diff_data->getPointer(2),
#endif
         dlower[0], dupper[0],
         dlower[1], dupper[1]
#if (NDIM == 3)
         , dlower[2], dupper[2]
#endif
         );
   }
}

//=======================================================================

// Computes phase gradients at cell centers.

// phase_diffs_id is SideData with NGHOSTS
// phase_grad_id is CellData with no ghosts with depth NDIM

void QuatModel::computePhaseGradCell(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& phase_diffs_id,
   int& phase_grad_cell_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for ( int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computePhaseGradCell(
         patch_level, phase_diffs_id, phase_grad_cell_id, time );
   }
}

void QuatModel::computePhaseGradCell(
   const boost::shared_ptr< hier::PatchLevel > patch_level,
   int& phase_diffs_id,
   int& phase_grad_cell_id,
   const double time )
{
   if ( phase_diffs_id < 0 ) phase_diffs_id = d_phase_diffs_id;
   if ( phase_grad_cell_id < 0 ) phase_grad_cell_id = d_phase_grad_cell_id;

   computeVarGradCell( patch_level, phase_diffs_id, phase_grad_cell_id, time );
}

//=======================================================================

// Computes eta gradients at cell centers.

// eta_diffs_id is SideData with NGHOSTS
// eta_grad_id is CellData with no ghosts with depth NDIM

void QuatModel::computeEtaGradCell(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& eta_diffs_id,
   int& eta_grad_cell_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for ( int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeEtaGradCell(
         patch_level, eta_diffs_id, eta_grad_cell_id, time );
   }
}

void QuatModel::computeEtaGradCell(
   const boost::shared_ptr< hier::PatchLevel > patch_level,
   int& eta_diffs_id,
   int& eta_grad_cell_id,
   const double time )
{
   if ( eta_diffs_id < 0 ) eta_diffs_id = d_eta_diffs_id;
   if ( eta_grad_cell_id < 0 ) eta_grad_cell_id = d_eta_grad_cell_id;

   computeVarGradCell( patch_level, eta_diffs_id, eta_grad_cell_id, time );
}

//=======================================================================

// Computes gradients at cell centers.

// diff_id is SideData with NGHOSTS
// grad_id is CellData with no ghosts with depth NDIM

void QuatModel::computeVarGradCell(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& diffs_id,
   int& grad_cell_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeVarGradCell(
         patch_level, diffs_id, grad_cell_id, time );
   }
}

void QuatModel::computeVarGradCell(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& diffs_id,
   int& grad_cell_id,
   const double time )
{
   (void)time;

   assert ( diffs_id >= 0 );
   assert ( grad_cell_id >= 0 );

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::SideData<double> > diff_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( diffs_id) ) );
      assert( diff_data );
      assert( diff_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::CellData<double> > grad_cell_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( grad_cell_id) ) );
      assert( grad_cell_data );
      assert( grad_cell_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      const hier::Box & diff_gbox = diff_data->getGhostBox();
      const hier::Index& d_lower = diff_gbox.lower();
      const hier::Index& d_upper = diff_gbox.upper();

      const hier::Box & grad_gbox = grad_cell_data->getGhostBox();
      const hier::Index& g_lower = grad_gbox.lower();
      const hier::Index& g_upper = grad_gbox.upper();

      boost::shared_ptr< geom::CartesianPatchGeometry > patch_geom (
         BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(patch->getPatchGeometry()) );
      TBOX_ASSERT(patch_geom);

      const double * dx = patch_geom->getDx();

      assert( grad_cell_data->getDepth() == NDIM );

      FORT_GRAD_CELL(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         diff_data->getPointer(0),
         diff_data->getPointer(1),
#if (NDIM == 3)
         diff_data->getPointer(2),
#endif
         d_lower[0], d_upper[0],
         d_lower[1], d_upper[1],
#if (NDIM == 3)
         d_lower[2], d_upper[2],
#endif
         dx,
         grad_cell_data->getPointer(0),
         grad_cell_data->getPointer(1),
#if (NDIM == 3)
         grad_cell_data->getPointer(2),
#endif
         g_lower[0], g_upper[0],
         g_lower[1], g_upper[1]
#if (NDIM == 3)
         , g_lower[2], g_upper[2]
#endif
         );
   }
}

//-----------------------------------------------------------------------

// Computes gradients at sides.

// diff_id is SideData with NGHOSTS
// grad_id is SideData with no ghosts, depth of NDIM

void QuatModel::computeVarGradSide(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& diffs_id,
   int& grad_side_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeVarGradSide(
         patch_level, diffs_id, grad_side_id, time );
   }
}

void QuatModel::computeVarGradSide(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& diffs_id,
   int& grad_side_id,
   const double time )
{
   (void)time;

   if ( diffs_id < 0 ) diffs_id = d_phase_diffs_id;
   if ( grad_side_id < 0 ) grad_side_id = d_phase_grad_side_id;

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::SideData<double> > diff_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(
             patch->getPatchData( diffs_id) ) );
      assert( diff_data );
      assert( diff_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::SideData<double> > grad_side_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(
             patch->getPatchData( grad_side_id) ) );
      assert( grad_side_data );
      assert( grad_side_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      const hier::Box & diff_gbox = diff_data->getGhostBox();
      const hier::Index& d_lower = diff_gbox.lower();
      const hier::Index& d_upper = diff_gbox.upper();

      const hier::Box & grad_gbox = grad_side_data->getGhostBox();
      const hier::Index& g_lower = grad_gbox.lower();
      const hier::Index& g_upper = grad_gbox.upper();

      boost::shared_ptr< geom::CartesianPatchGeometry > patch_geom ( 
         BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(
             patch->getPatchGeometry()) );
      TBOX_ASSERT(patch_geom);

      const double * dx = patch_geom->getDx();

      // there is a gradient component for each dimension x,y,z
      assert( grad_side_data->getDepth() == NDIM );

      FORT_GRAD_SIDE(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         diff_data->getPointer(0),
         diff_data->getPointer(1),
#if (NDIM == 3)
         diff_data->getPointer(2),
#endif
         d_lower[0], d_upper[0],
         d_lower[1], d_upper[1],
#if (NDIM == 3)
         d_lower[2], d_upper[2],
#endif
         dx,
         grad_side_data->getPointer(0,0), //side 0, depth 0 (x component)
         grad_side_data->getPointer(0,1), //side 0, depth 1 (y component)
#if (NDIM == 3)
         grad_side_data->getPointer(0,2),
#endif
         grad_side_data->getPointer(1,0),
         grad_side_data->getPointer(1,1),
#if (NDIM == 3)
         grad_side_data->getPointer(1,2),
#endif
#if (NDIM == 3)
         grad_side_data->getPointer(2,0),
         grad_side_data->getPointer(2,1),
         grad_side_data->getPointer(2,2),
#endif
         g_lower[0], g_upper[0],
         g_lower[1], g_upper[1]
#if (NDIM == 3)
         , g_lower[2], g_upper[2]
#endif
         );
   }
}

//=======================================================================

// Computes differences at cell sides. This will always compute the
// "normal" diffs, and will also compute the symmetry-aware diffs
// if symmetry is on.

// quat_id is CellData with NGHOSTS
// diff_id is SideData with NGHOSTS

void QuatModel::computeQuatDiffs(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& quat_id,
   int& quat_diffs_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatDiffs( patch_level, quat_id, quat_diffs_id, time );
   }
   
   assert( quat_diffs_id >= 0 );
}

void QuatModel::computeQuatDiffs(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& quat_id,
   int& quat_diffs_id,
   const double time )
{
   (void)time;

   if ( quat_id < 0 ) quat_id = d_quat_scratch_id;
   if ( quat_diffs_id < 0 ) quat_diffs_id = d_quat_diffs_id;

   assert( quat_id >= 0 );
   assert( quat_diffs_id >= 0 );

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& box = patch->getBox();
      const hier::Index& ifirst = box.lower();
      const hier::Index& ilast = box.upper();

      boost::shared_ptr< pdat::CellData<double> > quat_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
             patch->getPatchData( quat_id) ) );
      assert( quat_data );
      assert( quat_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
      assert( quat_data->getDepth() == d_qlen );

      boost::shared_ptr< pdat::SideData<double> > diff_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(
             patch->getPatchData( quat_diffs_id) ) );
      assert( diff_data );
      assert( diff_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      // If symmetry is on, "symmetric diffs" are stored first, with
      // "nonsymmetric diffs" stored offset by d_qlen.
      const int symm_depth_offset    = 0;
      const int nonsymm_depth_offset = isSymmetryAware() ? d_qlen : 0;

      FORT_QUATDIFFS(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         d_qlen,
         quat_data->getPointer(),
         quat_data->getGhostCellWidth()[0],
         diff_data->getPointer( 0, nonsymm_depth_offset ),
         diff_data->getPointer( 1, nonsymm_depth_offset ),
#if (NDIM == 3)
         diff_data->getPointer( 2, nonsymm_depth_offset ),
#endif
         diff_data->getGhostCellWidth()[0]
         );

      if ( d_symmetry_aware ) {

         assert( d_quat_symm_rotation_id>=0 );
         boost::shared_ptr< pdat::SideData<int> > rotation_index (
            BOOST_CAST< pdat::SideData<int>, hier::PatchData>(
                patch->getPatchData( d_quat_symm_rotation_id) ) );
         assert( rotation_index );
         assert( rotation_index->getGhostCellWidth() ==
                 hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

         FORT_QUATDIFFS_SYMM(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen,
            quat_data->getPointer(),
            quat_data->getGhostCellWidth()[0],
            diff_data->getPointer( 0, symm_depth_offset ),
            diff_data->getPointer( 1, symm_depth_offset ),
#if (NDIM == 3)
            diff_data->getPointer( 2, symm_depth_offset ),
#endif
            diff_data->getGhostCellWidth()[0],
            rotation_index->getPointer( 0 ),
            rotation_index->getPointer( 1 ),
#if (NDIM == 3)
            rotation_index->getPointer( 2 ),
#endif
            rotation_index->getGhostCellWidth()[0]
            );
      }

      if ( d_model_parameters.with_extra_visit_output() ) {
         boost::shared_ptr< pdat::CellData<double> > cell_diffs_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData( d_quat_diffs_cell_id) ) );
         assert( cell_diffs_data );

         boost::shared_ptr< pdat::CellData<double> > nonsymm_cell_diffs_data (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData( d_quat_nonsymm_diffs_cell_id) ) );
         assert( nonsymm_cell_diffs_data );

         pdat::CellIterator iend(pdat::CellGeometry::end(box));
         for ( pdat::CellIterator i(pdat::CellGeometry::begin(box)); i!=iend; ++i ) {
            const pdat::CellIndex ccell = *i;
            const pdat::SideIndex xside(
               ccell, pdat::SideIndex::X, pdat::SideIndex::Lower );
            const pdat::SideIndex yside(
               ccell, pdat::SideIndex::Y, pdat::SideIndex::Lower );
#if (NDIM == 3)
            const pdat::SideIndex zside(
               ccell, pdat::SideIndex::Z, pdat::SideIndex::Lower );
#endif
            for ( int q = 0; q < d_qlen; q++ ) {
               int d = symm_depth_offset + q;
               (*cell_diffs_data)(ccell,q+0*d_qlen) = (*diff_data)(xside,d);
               (*cell_diffs_data)(ccell,q+1*d_qlen) = (*diff_data)(yside,d);
#if (NDIM == 3)
               (*cell_diffs_data)(ccell,q+2*d_qlen) = (*diff_data)(zside,d);
#endif
               d = nonsymm_depth_offset + q;
               (*nonsymm_cell_diffs_data)(ccell,q+0*d_qlen) = (*diff_data)(xside,d);
               (*nonsymm_cell_diffs_data)(ccell,q+1*d_qlen) = (*diff_data)(yside,d);
#if (NDIM == 3)
               (*nonsymm_cell_diffs_data)(ccell,q+2*d_qlen) = (*diff_data)(zside,d);
#endif
            }
         }
      }
   }
}

//-----------------------------------------------------------------------

// Computes gradients at cell centers. This will always compute the
// symmetry-aware gradient if symmetry is on.

// diff_id is SideData with NGHOSTS, depth of NDIM*d_qlen
// grad_id is CellData with no ghosts

void QuatModel::computeQuatGradCell(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& quat_diffs_id,
   int& grad_cell_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatGradCell(
         patch_level, quat_diffs_id, grad_cell_id, time );
   }
}

void QuatModel::computeQuatGradCell(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& quat_diffs_id,
   int& grad_cell_id,
   const double time )
{
   (void)time;
   
   if ( quat_diffs_id < 0 ) quat_diffs_id = d_quat_diffs_id;
   if ( grad_cell_id < 0 ) grad_cell_id = d_quat_grad_cell_id;
   assert( quat_diffs_id>=0 );
   assert( grad_cell_id>=0 );

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::SideData<double> > diff_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( quat_diffs_id) ) );
      assert( diff_data );
      assert( diff_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::CellData<double> > grad_cell_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( grad_cell_id) ) );
      assert( grad_cell_data );
      assert( grad_cell_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( grad_cell_data->getDepth() == NDIM * d_qlen );

      boost::shared_ptr< geom::CartesianPatchGeometry > patch_geom ( 
         BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(patch->getPatchGeometry()) );
      TBOX_ASSERT(patch_geom);

      const double * dx = patch_geom->getDx();

      if ( d_symmetry_aware ) {

         boost::shared_ptr< pdat::SideData<int> > rotation_index (
            BOOST_CAST< pdat::SideData<int>, hier::PatchData>(patch->getPatchData( d_quat_symm_rotation_id) ) );
         assert( rotation_index );
         assert( rotation_index->getGhostCellWidth() ==
                 hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

         FORT_QUATGRAD_CELL_SYMM(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen, dx,
            diff_data->getPointer( 0, 0 ),
            diff_data->getPointer( 1, 0 ),
#if (NDIM == 3)
            diff_data->getPointer( 2, 0 ),
#endif
            NGHOSTS,
            grad_cell_data->getPointer( 0 * d_qlen ),
            grad_cell_data->getPointer( 1 * d_qlen ),
#if (NDIM == 3)
            grad_cell_data->getPointer( 2 * d_qlen ),
#endif
            0,
            rotation_index->getPointer( 0 ),
            rotation_index->getPointer( 1 ),
#if (NDIM == 3)
            rotation_index->getPointer( 2 ),
#endif
            NGHOSTS
            );

      }
      else {

         FORT_QUATGRAD_CELL(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen, dx,
            diff_data->getPointer( 0, 0 ),
            diff_data->getPointer( 1, 0 ),
#if (NDIM == 3)
            diff_data->getPointer( 2, 0 ),
#endif
            NGHOSTS,
            grad_cell_data->getPointer( 0 * d_qlen ),
            grad_cell_data->getPointer( 1 * d_qlen ),
#if (NDIM == 3)
            grad_cell_data->getPointer( 2 * d_qlen ),
#endif
            0
            );

      }
   }
}

//-----------------------------------------------------------------------

// Computes gradients at side.

// diff_id is SideData with NGHOSTS, depth of NDIM*d_qlen
// grad_id is SideData with no ghosts

void QuatModel::computeQuatGradSide(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& quat_diffs_id,
   int& grad_side_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatGradSide(
         patch_level, quat_diffs_id, grad_side_id, time );
   }
}

/*
 * Compute gradients, using differences
 */
void QuatModel::computeQuatGradSide(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& quat_diffs_id,
   int& grad_side_id,
   const double time )
{
   (void)time;

   //tbox::pout<<"computeQuatGradSide()"<<endl;
   
   if ( quat_diffs_id < 0 ) quat_diffs_id = d_quat_diffs_id;
   if ( grad_side_id < 0 ) grad_side_id = d_quat_grad_side_id;
   assert( grad_side_id>=0 );

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::SideData<double> > diff_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( quat_diffs_id) ) );
      assert( diff_data );
      assert( diff_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::SideData<double> > grad_side_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( grad_side_id) ) );
      assert( grad_side_data );
      assert( grad_side_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );
      assert( grad_side_data->getDepth() == NDIM * d_qlen );

      boost::shared_ptr< geom::CartesianPatchGeometry > patch_geom ( 
         BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(patch->getPatchGeometry()) );
      TBOX_ASSERT(patch_geom);

      const double * dx = patch_geom->getDx();

      if ( d_symmetry_aware ) {

         boost::shared_ptr< pdat::SideData<int> > rotation_index (
            BOOST_CAST< pdat::SideData<int>, hier::PatchData>(patch->getPatchData( d_quat_symm_rotation_id) ) );
         assert( rotation_index );
         assert( rotation_index->getGhostCellWidth() ==
                 hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

         if( !d_model_parameters.useIsotropicStencil() ){
         
         //tbox::pout<<"compute quat grad at sides..."<<endl;
         
         FORT_QUATGRAD_SIDE_SYMM(
            ifirst(0), ilast(0), 
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen, dx,
            diff_data->getPointer( 0, 0 ),
            diff_data->getPointer( 1, 0 ),
#if (NDIM == 3)
            diff_data->getPointer( 2, 0 ),
#endif
            NGHOSTS,
            grad_side_data->getPointer( 0, 0 * d_qlen ), // grad_x, xside
            grad_side_data->getPointer( 0, 1 * d_qlen ), // grad_y, xside
#if (NDIM == 3)
            grad_side_data->getPointer( 0, 2 * d_qlen ), // grad_z, xside
#endif
            grad_side_data->getPointer( 1, 0 * d_qlen ), // grad_x, yside
            grad_side_data->getPointer( 1, 1 * d_qlen ), // grad_y, yside
#if (NDIM == 3)
            grad_side_data->getPointer( 1, 2 * d_qlen ), // grad_z, yside
#endif
#if (NDIM == 3)
            grad_side_data->getPointer( 2, 0 * d_qlen ),
            grad_side_data->getPointer( 2, 1 * d_qlen ),
            grad_side_data->getPointer( 2, 2 * d_qlen ),
#endif
            0,
            rotation_index->getPointer(0),
            rotation_index->getPointer(1),
#if (NDIM == 3)
            rotation_index->getPointer(2),
#endif
            NGHOSTS
            );

         }else{
      
         //tbox::pout<<"compute quat grad at sides using wide stencil..."<<endl;
         
         FORT_QUATGRAD_SIDE_SYMM_ISOTROPIC(
            ifirst(0), ilast(0), 
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen, dx,
            diff_data->getPointer( 0, 0 ),
            diff_data->getPointer( 1, 0 ),
#if (NDIM == 3)
            diff_data->getPointer( 2, 0 ),
#endif
            NGHOSTS,
            grad_side_data->getPointer( 0, 0 * d_qlen ), // grad_x, xside
            grad_side_data->getPointer( 0, 1 * d_qlen ), // grad_y, xside
#if (NDIM == 3)
            grad_side_data->getPointer( 0, 2 * d_qlen ), // grad_z, xside
#endif
            grad_side_data->getPointer( 1, 0 * d_qlen ), // grad_x, yside
            grad_side_data->getPointer( 1, 1 * d_qlen ), // grad_y, yside
#if (NDIM == 3)
            grad_side_data->getPointer( 1, 2 * d_qlen ), // grad_z, yside
#endif
#if (NDIM == 3)
            grad_side_data->getPointer( 2, 0 * d_qlen ),
            grad_side_data->getPointer( 2, 1 * d_qlen ),
            grad_side_data->getPointer( 2, 2 * d_qlen ),
#endif
            0,
            rotation_index->getPointer(0),
            rotation_index->getPointer(1),
#if (NDIM == 3)
            rotation_index->getPointer(2),
#endif
            NGHOSTS
            );
         }
      }
      else {

         if( !d_model_parameters.useIsotropicStencil() ){
         FORT_QUATGRAD_SIDE(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen, dx,
            diff_data->getPointer( 0, 0 ),
            diff_data->getPointer( 1, 0 ),
#if (NDIM == 3)
            diff_data->getPointer( 2, 0 ),
#endif
            NGHOSTS,
            // output
            grad_side_data->getPointer( 0, 0 * d_qlen ), // grad_x, xside
            grad_side_data->getPointer( 0, 1 * d_qlen ), // grad_y, xside
#if (NDIM == 3)
            grad_side_data->getPointer( 0, 2 * d_qlen ), // grad_z, xside
#endif
            grad_side_data->getPointer( 1, 0 * d_qlen ), // grad_x, yside
            grad_side_data->getPointer( 1, 1 * d_qlen ), // grad_y, yside
#if (NDIM == 3)
            grad_side_data->getPointer( 1, 2 * d_qlen ), // grad_z, yside
#endif
#if (NDIM == 3)
            grad_side_data->getPointer( 2, 0 * d_qlen ),
            grad_side_data->getPointer( 2, 1 * d_qlen ),
            grad_side_data->getPointer( 2, 2 * d_qlen ),
#endif
            0
            );
         }else{
         FORT_QUATGRAD_SIDE_ISOTROPIC(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen, dx,
            diff_data->getPointer( 0, 0 ),
            diff_data->getPointer( 1, 0 ),
#if (NDIM == 3)
            diff_data->getPointer( 2, 0 ),
#endif
            NGHOSTS,
            // output
            grad_side_data->getPointer( 0, 0 * d_qlen ), // grad_x, xside
            grad_side_data->getPointer( 0, 1 * d_qlen ), // grad_y, xside
#if (NDIM == 3)
            grad_side_data->getPointer( 0, 2 * d_qlen ), // grad_z, xside
#endif
            grad_side_data->getPointer( 1, 0 * d_qlen ), // grad_x, yside
            grad_side_data->getPointer( 1, 1 * d_qlen ), // grad_y, yside
#if (NDIM == 3)
            grad_side_data->getPointer( 1, 2 * d_qlen ), // grad_z, yside
#endif
#if (NDIM == 3)
            grad_side_data->getPointer( 2, 0 * d_qlen ),
            grad_side_data->getPointer( 2, 1 * d_qlen ),
            grad_side_data->getPointer( 2, 2 * d_qlen ),
#endif
            0
            );
         
         }
      }
   }
}

//-----------------------------------------------------------------------

// Computes gradients at cell centers.

// grad_cell_id is CellData with no ghosts, depth of NDIM*d_qlen
// grad_modulus_id is CellData with no ghosts, depth of 1

void QuatModel::computeQuatGradModulus(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& grad_cell_id,
   int& grad_modulus_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache != FORCE ) return;
   old_time = time;

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatGradModulus(
         patch_level,
         grad_cell_id,
         grad_modulus_id,
         time );
   }
}

void QuatModel::computeQuatGradModulusFromSides(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& grad_side_id,
   int& grad_modulus_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache != FORCE ) return;
   old_time = time;

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatGradModulusFromSides(
         patch_level,
         grad_side_id,
         grad_modulus_id,
         time );
   }
}

void QuatModel::computeQuatGradModulus(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& grad_cell_id,
   int& grad_modulus_id,
   const double time )
{
   (void)time;

   if ( grad_cell_id < 0 ) grad_cell_id = d_quat_grad_cell_id;
   if ( grad_modulus_id < 0 )
      grad_modulus_id = d_quat_grad_modulus_id;

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::CellData<double> > grad_cell_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( grad_cell_id) ) );
      assert( grad_cell_data );
      assert( grad_cell_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      boost::shared_ptr< pdat::CellData<double> > grad_modulus_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( grad_modulus_id) ) );
      assert( grad_modulus_data );
      assert( grad_modulus_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      const hier::Box & grad_gbox = grad_cell_data->getGhostBox();
      const hier::Index& g_lower = grad_gbox.lower();
      const hier::Index& g_upper = grad_gbox.upper();

      const hier::Box & mod_gbox = grad_modulus_data->getGhostBox();
      const hier::Index& m_lower = mod_gbox.lower();
      const hier::Index& m_upper = mod_gbox.upper();

      assert( grad_cell_data->getDepth() == NDIM * d_qlen );
      assert( grad_modulus_data->getDepth() == 1 );

      FORT_QUATGRAD_MODULUS(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         d_qlen,
         grad_cell_data->getPointer( 0 * d_qlen ),
         grad_cell_data->getPointer( 1 * d_qlen ),
#if (NDIM == 3)
         grad_cell_data->getPointer( 2 * d_qlen ),
#endif
         g_lower[0], g_upper[0],
         g_lower[1], g_upper[1],
#if (NDIM == 3)
         g_lower[2], g_upper[2],
#endif
         grad_modulus_data->getPointer(),
         m_lower[0], m_upper[0],
         m_lower[1], m_upper[1],
#if (NDIM == 3)
         m_lower[2], m_upper[2],
#endif
         d_model_parameters.quat_grad_floor_type().c_str(),
         d_model_parameters.quat_grad_floor()
         );

   }
}

void QuatModel::computeQuatGradModulusFromSides(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& grad_side_id,
   int& grad_modulus_id,
   const double time )
{
   (void)time;

   if ( grad_side_id < 0 ) grad_side_id = d_quat_grad_side_id;
   assert( grad_side_id>=0 );
   if ( grad_modulus_id < 0 )
      grad_modulus_id = d_quat_grad_modulus_id;

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::SideData<double> > grad_side_data (
         BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch->getPatchData( grad_side_id) ) );
      assert( grad_side_data );
      assert( grad_side_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      boost::shared_ptr< pdat::CellData<double> > grad_modulus_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( grad_modulus_id) ) );
      assert( grad_modulus_data );
      assert( grad_modulus_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      const hier::Box & grad_gbox = grad_side_data->getGhostBox();
      const hier::Index& g_lower = grad_gbox.lower();
      const hier::Index& g_upper = grad_gbox.upper();

      const hier::Box & mod_gbox = grad_modulus_data->getGhostBox();
      const hier::Index& m_lower = mod_gbox.lower();
      const hier::Index& m_upper = mod_gbox.upper();

      assert( grad_side_data->getDepth() == NDIM * d_qlen );
      assert( grad_modulus_data->getDepth() == 1 );

      FORT_QUATGRAD_MODULUS_FROM_SIDES_COMPACT(
//      FORT_QUATGRAD_MODULUS_FROM_SIDES(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         d_qlen,
         grad_side_data->getPointer( 0 ),
         grad_side_data->getPointer( 1 ),
#if (NDIM == 3)
         grad_side_data->getPointer( 2 ),
#endif
         g_lower[0], g_upper[0],
         g_lower[1], g_upper[1],
#if (NDIM == 3)
         g_lower[2], g_upper[2],
#endif
         grad_modulus_data->getPointer(),
         m_lower[0], m_upper[0],
         m_lower[1], m_upper[1],
#if (NDIM == 3)
         m_lower[2], m_upper[2],
#endif
         d_model_parameters.quat_grad_floor_type().c_str(),
         d_model_parameters.quat_grad_floor()
         );

   }
}

//=======================================================================

void QuatModel::normalizeQuat(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int quat_id )
{
   if ( quat_id < 0 ) return;
   if ( d_qlen == 1 ) return;
   
   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      normalizeQuat( patch_level, quat_id );
   }
}

void QuatModel::normalizeQuat(
   const boost::shared_ptr< hier::PatchLevel > level,
   const int quat_id )
{
   assert( quat_id >= 0 );
   if ( d_qlen == 1 ) return;
   
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;

      boost::shared_ptr< pdat::CellData<double> > quat (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( quat_id) ) );
      assert( quat );
      const hier::Box& gbox = quat->getGhostBox();

      // issue:mew: Potentially replace this loop with fortran kernel.

      pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
      for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox)); i!=iend; ++i) {
         pdat::CellIndex cell = *i;
         double qnorm2 = 0.;
         for ( int q = 0; q < d_qlen; q++ ) {
            qnorm2 += (*quat)(cell,q) * (*quat)(cell,q);
         }
         const double invqnorm = 1. / sqrt(qnorm2);
         for ( int q = 0; q < d_qlen; q++ ) {
            (*quat)(cell,q) *= invqnorm;
         }
      }
   }
}

//=======================================================================

void QuatModel::applyPolynomial(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int src_cell_data_id,
   const int dst_cell_data_id )

{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      applyPolynomial( patch_level, src_cell_data_id, dst_cell_data_id );
   }
}

void QuatModel::applyPolynomial(
   const boost::shared_ptr< hier::PatchLevel > level,
   const int src_cell_data_id,
   const int dst_cell_data_id )

{
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;

      boost::shared_ptr< pdat::CellData<double> > sdata (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( src_cell_data_id) ) );
      assert( sdata );
      const hier::Box& gbox = sdata->getGhostBox();

      boost::shared_ptr< pdat::CellData<double> > ddata (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( dst_cell_data_id) ) );
      assert( ddata );

      // issue:mew: Potentially replace this loop with fortran kernel.
      const char interp = energyInterpChar(d_model_parameters.energy_interp_func_type());
      pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
      for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox)); i!=iend; ++i) {
         pdat::CellIndex cell = *i;
         const double phi = (*sdata)(cell);
         const double hphi =
             FORT_INTERP_FUNC( phi, &interp );
         (*ddata)(cell)=hphi;
      }
   }
}

//=======================================================================

// Computes phase mobility at cell centers.

// phase_id is CellData with NGHOSTS
// mobility_id is CellData with no ghosts

void QuatModel::computeUniformPhaseMobility(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& phase_id,
   int& mobility_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeUniformPhaseMobility(
         patch_level, phase_id, mobility_id, time );
   }
}

void QuatModel::computeUniformPhaseMobility(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& phase_id,
   int& mobility_id,
   const double time )
{
   (void)time;

   if ( phase_id < 0 ) phase_id = d_phase_scratch_id;
   if ( mobility_id < 0 ) mobility_id = d_phase_mobility_id;

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {

      boost::shared_ptr<hier::Patch > patch = *p;

      boost::shared_ptr< pdat::CellData<double> > mobility_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( mobility_id) ) );
      assert( mobility_data );
      assert( mobility_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),1) );

      mobility_data->fillAll( d_model_parameters.phase_mobility() );

   }
}

//-----------------------------------------------------------------------

void QuatModel::computePhaseMobilityPatch(
   const hier::Box& pbox,
   boost::shared_ptr< pdat::CellData<double> > cd_temp,
   boost::shared_ptr< pdat::CellData<double> > cd_mobility )
{
   double* ptr_temp = cd_temp->getPointer();
   double* ptr_m = cd_mobility->getPointer();
   
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

   const hier::Box& m_gbox = cd_mobility->getGhostBox();
   int imin_m = m_gbox.lower(0);
   int jmin_m = m_gbox.lower(1);
   int jp_m = m_gbox.numberCells(0);
   int kmin_m = 0;
   int kp_m = 0;
#if (NDIM == 3)
   kmin_m = m_gbox.lower(2);
   kp_m = jp_m * m_gbox.numberCells(1);
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

            const int idx_m = (ii - imin_m) +
               (jj - jmin_m) * jp_m + (kk - kmin_m) * kp_m;

            double t = ptr_temp[idx_temp];

            double Rt = t * gas_constant_R_JpKpmol;

            ptr_m[idx_m] =
               d_model_parameters.phase_mobility() 
             * exp( -d_model_parameters.q0_phase_mobility() / Rt );

         }
      }
   }
}

//=======================================================================

// Computes eta mobility at cell centers.

// eta_id is CellData with NGHOSTS
// mobility_id is CellData with no ghosts

void QuatModel::computeEtaMobility(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& phase_id,
   int& mobility_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeEtaMobility(
         patch_level, phase_id, mobility_id, time );
   }
}

void QuatModel::computeEtaMobility(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& phase_id,
   int& mobility_id,
   const double time )
{
   (void)time;

   if ( phase_id < 0 ) phase_id = d_phase_scratch_id;
   if ( mobility_id < 0 ) mobility_id = d_eta_mobility_id;

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;

      const hier::Box& pbox = patch->getBox();
      
      boost::shared_ptr< pdat::CellData<double> > temp_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_temperature_id) ) );
      assert( temp_data );

      boost::shared_ptr< pdat::CellData<double> > phase_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( phase_id) ) );
      assert( phase_data );

      boost::shared_ptr< pdat::CellData<double> > mobility_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( mobility_id) ) );
      assert( mobility_data );

      if ( d_model_parameters.min_eta_mobility() == d_model_parameters.eta_mobility() &&
           d_model_parameters.q0_eta_mobility() == 0.0 ) {

         mobility_data->fillAll( d_model_parameters.eta_mobility() );

      }
      else {

         computeEtaMobilityPatch(
            pbox,
            temp_data,
            mobility_data,
            phase_data );

      }
   }
}

//-----------------------------------------------------------------------

void QuatModel::computeEtaMobilityPatch(
   const hier::Box& pbox,
   boost::shared_ptr< pdat::CellData<double> > cd_temp,
   boost::shared_ptr< pdat::CellData<double> > cd_mobility,
   boost::shared_ptr< pdat::CellData<double> > cd_phi )
{
   double* ptr_temp = cd_temp->getPointer();
   double* ptr_m = cd_mobility->getPointer();
   double* ptr_phi = cd_phi->getPointer();
   
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

   const hier::Box& m_gbox = cd_mobility->getGhostBox();
   int imin_m = m_gbox.lower(0);
   int jmin_m = m_gbox.lower(1);
   int jp_m = m_gbox.numberCells(0);
   int kmin_m = 0;
   int kp_m = 0;
#if (NDIM == 3)
   kmin_m = m_gbox.lower(2);
   kp_m = jp_m * m_gbox.numberCells(1);
#endif

   const hier::Box& phi_gbox = cd_phi->getGhostBox();
   int imin_phi = phi_gbox.lower(0);
   int jmin_phi = phi_gbox.lower(1);
   int jp_phi = phi_gbox.numberCells(0);
   int kmin_phi = 0;
   int kp_phi = 0;
#if (NDIM == 3)
   kmin_phi = phi_gbox.lower(2);
   kp_phi = jp_phi * phi_gbox.numberCells(1);
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

            const int idx_m = (ii - imin_m) +
               (jj - jmin_m) * jp_m + (kk - kmin_m) * kp_m;

            const int idx_phi = (ii - imin_phi) +
               (jj - jmin_phi) * jp_phi + (kk - kmin_phi) * kp_phi;

            double t = ptr_temp[idx_temp];
            double phi = ptr_phi[idx_phi];

            double h_phi = FORT_INTERP_FUNC( phi, "p" );

            double Rt = t * gas_constant_R_JpKpmol;

            double m = d_model_parameters.eta_mobility() 
                     * exp( -d_model_parameters.q0_eta_mobility() / Rt );

            ptr_m[idx_m] =
               m - h_phi * ( m - d_model_parameters.min_eta_mobility() );

         }
      }
   }
}

//=======================================================================

// Computes Quaternion mobility at cell centers.

// phase_id is CellData with NGHOSTS
// mobility_id is CellData with no ghosts

void QuatModel::computeQuatMobility(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& phase_id,
   int& mobility_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatMobility(
         patch_level, phase_id, mobility_id, time );
   }
}

void QuatModel::computeQuatMobility(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& phase_id,
   int& mobility_id,
   const double time )
{
   (void)time;

   if ( phase_id < 0 ) phase_id = d_phase_scratch_id;
   if ( mobility_id < 0 ) mobility_id = d_quat_mobility_id;

   double alt_scale_factor=d_model_parameters.quatMobilityScaleFactor();
   
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      boost::shared_ptr<hier::Patch > patch = *p;

      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::CellData<double> > phase_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( phase_id) ) );
      assert( phase_data );
      assert( phase_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::CellData<double> > mobility_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( mobility_id) ) );
      assert( mobility_data );
      assert( mobility_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      FORT_QUATMOBILITY(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         phase_data->getPointer(),
         NGHOSTS,
         mobility_data->getPointer(),
         NGHOSTS,
         d_model_parameters.quat_mobility(),
         d_model_parameters.min_quat_mobility(),
         d_model_parameters.quat_mobility_func_type().c_str(),
         alt_scale_factor );

   }

}

//=======================================================================

// Computes Derivative of Quaternion mobility versus Phase at cell centers.

// phase_id is CellData with NGHOSTS
// mobility_deriv_id is CellData with no ghosts

void QuatModel::computeQuatMobilityDeriv(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int& phase_id,
   int& mobility_deriv_id,
   const double time,
   const CACHE_TYPE cache )
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if ( time == old_time && cache == CACHE ) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      boost::shared_ptr<hier::PatchLevel > patch_level =
         hierarchy->getPatchLevel( ln );

      computeQuatMobilityDeriv(
         patch_level, phase_id, mobility_deriv_id, time );
   }
}

void QuatModel::computeQuatMobilityDeriv(
   const boost::shared_ptr< hier::PatchLevel > level,
   int& phase_id,
   int& mobility_deriv_id,
   const double time )
{
   (void)time;

   if ( phase_id < 0 ) phase_id = d_phase_scratch_id;
   assert( mobility_deriv_id >= 0 );

   double alt_scale_factor=d_model_parameters.quatMobilityScaleFactor();

   for ( hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p ) {
      boost::shared_ptr<hier::Patch > patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast =  pbox.upper();

      boost::shared_ptr< pdat::CellData<double> > phase_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( phase_id) ) );
      assert( phase_data );
      assert( phase_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

      boost::shared_ptr< pdat::CellData<double> > mobility_deriv_data (
         BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( mobility_deriv_id) ) );
      assert( mobility_deriv_data );
      assert( mobility_deriv_data->getGhostCellWidth() ==
              hier::IntVector(tbox::Dimension(NDIM),0) );

      FORT_QUATMOBILITYDERIV(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         phase_data->getPointer(),
         NGHOSTS,
         mobility_deriv_data->getPointer(),
         0,
         d_model_parameters.quat_mobility(),
         d_model_parameters.min_quat_mobility(),
         d_model_parameters.quat_mobility_func_type().c_str(),
         alt_scale_factor );
   }
}

//=======================================================================

void QuatModel::checkQuatNorm(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const double tol)
{
   assert( tol >= 0. );
   
   if ( ! d_model_parameters.with_orientation() )return;
   if ( d_qlen == 1 ) return;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++ ) {

      boost::shared_ptr<hier::PatchLevel >
         level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
         boost::shared_ptr<hier::Patch > patch = *p;

         boost::shared_ptr< pdat::CellData<double> > y (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData(d_quat_id) ) );
         const hier::Box& pbox = patch->getBox();

         pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox)); ic != icend; ++ic) {
            pdat::CellIndex cell = *ic;
            double qnorm2 = 0.;
            for ( int q = 0; q < d_qlen; q++ ) {
               qnorm2 += (*y)(cell,q) * (*y)(cell,q);
            }
            double qnorm = sqrt(qnorm2);
            if ( fabs(qnorm-1.) > tol ) {
               cerr << setprecision(10) << scientific;
               cerr << "WARNING: q norm=" << qnorm << endl;
            }
         }
      }
   }
}

/*
 * Set weight appropriate for computing vector norms.
 *
 * If you this function to set the weights used when you
 * SAMRAIVectorReal::addComponent, you can use the
 * vector norm functions of SAMRAIVectorReal, and
 * the weights will be used to blank out coarse grid
 * regions under fine grids.
 *
 * The weights computed are specific to the cell-centered
 * discretization used by this class.  The weight is equal
 * to the cell volume if the cell has not been refined,
 * and zero if it has.
 *
 * This function is state-independent.  All inputs are in
 * the argument list.
 *
 * hierarchy:   Hierarchy configuration to compute weights for
 * coarsest_ln: Coarsest level number.  Must be included
 *              in hierarchy.  Must not be greater than finest_ln.
 *              Default to 0.
 * finest_ln:   Finest level number.  Must be included
 *              in hierarchy.  Must not be less than coarsest_ln.
 *              Default to finest level in hierarchy.
 */
void QuatModel::computeVectorWeights(
   boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   int coarsest_ln,
   int finest_ln )
{
   assert( d_weight_id!=-1 );
   
   if ( coarsest_ln == -1 ) coarsest_ln = 0;
   if ( finest_ln == -1 ) finest_ln = hierarchy->getFinestLevelNumber();
   if ( finest_ln < coarsest_ln ) {
      TBOX_ERROR(d_object_name
                 << ": Illegal level number range.  finest_ln < coarsest_ln.");
   }

   for ( int ln=finest_ln; ln >= coarsest_ln; --ln ) {

      /*
       * On every level, first assign cell volume to vector weight.
       */

      boost::shared_ptr< hier::PatchLevel > level =
         hierarchy->getPatchLevel(ln);
      for ( hier::PatchLevel::iterator p(level->begin());
            p != level->end(); ++p ) {
         boost::shared_ptr< hier::Patch > patch = *p;
         boost::shared_ptr< geom::CartesianPatchGeometry > patch_geometry (
            BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(
               patch->getPatchGeometry()) );
         TBOX_ASSERT(patch_geometry);

         const double* dx = patch_geometry->getDx();
         double cell_vol = dx[0];
	 if (NDIM > 1) {
	    cell_vol *= dx[1];
	 }
	 if (NDIM > 2) {
	    cell_vol *= dx[2];
	 }

         boost::shared_ptr< pdat::CellData<double> > w (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(d_weight_id) ) );
         if ( !w ) {
            TBOX_ERROR(d_object_name
                       << ": weight id must refer to a pdat::CellVariable");
         }
         w->fillAll(cell_vol);
      }

      /*
       * On all but the finest level, assign 0 to vector
       * weight to cells covered by finer cells.
       */

      if (ln < finest_ln) {

         /*
          * First get the boxes that describe index space of the next finer
          * level and coarsen them to describe corresponding index space
          * at this level.
          */

         boost::shared_ptr< hier::PatchLevel > next_finer_level =
            hierarchy->getPatchLevel(ln+1);
         hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
         hier::IntVector coarsen_ratio = next_finer_level->getRatioToLevelZero();
         coarsen_ratio /= level->getRatioToLevelZero();
         coarsened_boxes.coarsen( coarsen_ratio );

         /*
          * Then set vector weight to 0 wherever there is
          * a nonempty intersection with the next finer level.
          * Note that all assignments are local.
          */

         for ( hier::PatchLevel::iterator p(level->begin()); p != level->end(); ++p ) {

            boost::shared_ptr< hier::Patch > patch = *p;
            for (hier::BoxContainer::const_iterator i=coarsened_boxes.begin();
                 i != coarsened_boxes.end(); ++i) {

               hier::Box coarse_box = *i;
               hier::Box intersection = coarse_box*( patch->getBox() );
               if ( !intersection.empty() ) {
                  boost::shared_ptr< pdat::CellData<double> > w (
                     BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData(d_weight_id) ) );
                  w->fillAll(0.0, intersection);

               }  // assignment only in non-empty intersection
            }  // loop over coarsened boxes from finer level
         }  // loop over patches in level
      }  // all levels except finest
   }  // loop over levels
}

//=======================================================================

void QuatModel::evaluateEnergy(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const double time,
   double& total_energy,
   double& total_phase_e,
   double& total_eta_e,
   double& total_orient_e,
   double& total_qint_e,
   double& total_well_e,
   double& total_free_e,
   const bool gp )
{
   assert( d_weight_id != -1 );
   if( d_model_parameters.with_visit_energy_output() )
      assert( d_energy_diag_id != -1 );
   if( d_model_parameters.with_concentration() )
      assert( d_phase_conc_strategy!=nullptr );
   
   total_energy = 0.;
   total_phase_e = 0.;
   total_eta_e = 0.;
   total_orient_e = 0.;
   total_qint_e = 0.;
   total_free_e = 0.;
   total_well_e = 0.;
   
   copyCurrentToScratch( hierarchy, time, d_all_refine_patch_strategy );

   if ( d_model_parameters.with_orientation() ) {
      int diff_id = -1;
      d_quat_grad_strategy->computeDiffs(
         hierarchy,
         d_quat_scratch_id,
         diff_id,
         time,
         QuatGradStrategy::FORCE );
      assert( diff_id >= 0 );

      // Compute gradients on cell faces
      d_quat_grad_strategy->computeGradSide(
         hierarchy,
         diff_id,
         d_quat_grad_side_id,
         time,
         QuatGradStrategy::FORCE );
   }
   // compute free energies function of phase, concentration and temperature

   if( d_model_parameters.with_concentration() ){
      if( d_model_parameters.partition_phase_concentration() ){
         //computeVelocity(hierarchy,ydot_phase_id);
         
         assert( d_partition_coeff_strategy!=nullptr );
         
         d_partition_coeff_strategy->evaluate(hierarchy);
         if( d_model_parameters.needGhosts4PartitionCoeff() )
            fillPartitionCoeffGhosts();
      }
      d_phase_conc_strategy->computePhaseConcentrations(
         hierarchy,
         d_temperature_scratch_id,
         d_phase_scratch_id,
         d_eta_scratch_id,
         d_conc_scratch_id );
   }
   
   d_free_energy_strategy->computeFreeEnergyLiquid(
      hierarchy,
      d_temperature_id,
      d_f_l_id,
      gp );

   d_free_energy_strategy->computeFreeEnergySolidA(
      hierarchy, 
      d_temperature_id,
      d_f_a_id,
      gp );
   
   if ( d_model_parameters.with_third_phase() ) {
      d_free_energy_strategy->computeFreeEnergySolidB(
         hierarchy, 
         d_temperature_id,
         d_f_b_id,
         gp );
   }

   const double epsilon_anisotropy = d_model_parameters.epsilon_anisotropy();

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {

      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {

         boost::shared_ptr<hier::Patch > patch = *p;
         boost::shared_ptr< geom::CartesianPatchGeometry > patch_geom ( 
            BOOST_CAST< geom::CartesianPatchGeometry,
                        hier::PatchGeometry>(patch->getPatchGeometry()) );
         TBOX_ASSERT(patch_geom);

         const double * dx = patch_geom->getDx();

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast =  pbox.upper();

         double* pgrad_quat[NDIM];
         if ( d_model_parameters.with_orientation() ){
            boost::shared_ptr< pdat::SideData<double> > grad_quat (
               BOOST_CAST< pdat::SideData<double>,
                           hier::PatchData>(patch->getPatchData(
                              d_quat_grad_side_id) ) );
            assert( grad_quat );
            assert( grad_quat->getGhostCellWidth() == hier::IntVector(tbox::Dimension(NDIM),0) );
            for ( int d = 0; d < NDIM; d++ ) {
               pgrad_quat[d] = grad_quat->getPointer( d );
            }
#ifdef DEBUG_CHECK_ASSERTIONS
            SAMRAI::math::PatchSideDataNormOpsReal<double> sops;
            double l2gq=sops.L2Norm(grad_quat,pbox);
            assert( l2gq==l2gq );
#endif
         }
         else {
            for ( int d = 0; d < NDIM; d++ ) {
               pgrad_quat[d] = nullptr;
            }
         }

         if( d_model_parameters.with_phase() ){

         boost::shared_ptr< pdat::CellData<double> > phase (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_phase_scratch_id) ) );
         boost::shared_ptr< pdat::CellData<double> > weight (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_weight_id) ) );
         boost::shared_ptr< pdat::CellData<double> > fl (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_f_l_id) ) );
         boost::shared_ptr< pdat::CellData<double> > fa (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_f_a_id) ) );
         boost::shared_ptr< pdat::CellData<double> > temperature (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_temperature_id) ) );

         double* quat_ptr = nullptr;
         if( epsilon_anisotropy>=0. ){
            boost::shared_ptr< pdat::CellData<double> > quat (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
                  patch->getPatchData( d_quat_scratch_id) ) );
            quat_ptr = quat->getPointer();
            assert( quat_ptr != nullptr );
         }

         assert( phase );
         assert( weight );
         assert( fl );
         assert( fa );

#ifdef DEBUG_CHECK_ASSERTIONS
         SAMRAI::math::PatchCellDataNormOpsReal<double> ops;
         double l2phi=ops.L2Norm(phase,pbox);
         assert( l2phi==l2phi );

         double l2t=ops.L2Norm(temperature,pbox);
         assert( l2t==l2t );

         double l2fl=ops.L2Norm(fl,pbox);
         assert( l2fl==l2fl );

         double l2fa=ops.L2Norm(fa,pbox);
         assert( l2fa==l2fa );
#endif

         int three_phase = 0;
         double* ptr_fb = nullptr;
         double* ptr_eta = nullptr;
         if ( d_model_parameters.with_third_phase() ) {
            three_phase = 1;
            boost::shared_ptr< pdat::CellData<double> > fb (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_f_b_id) ) );
            ptr_fb = fb->getPointer();
            boost::shared_ptr< pdat::CellData<double> > eta (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_eta_scratch_id) ) );
            ptr_eta = eta->getPointer();
         }

         int per_cell = 0;
         double* ptr_energy = nullptr;
         if ( d_model_parameters.with_visit_energy_output() ) {
            per_cell = 1; 
            boost::shared_ptr< pdat::CellData<double> > energy (
               BOOST_CAST< pdat::CellData<double>,
                           hier::PatchData>(patch->getPatchData(
                              d_energy_diag_id) ) );
            ptr_energy = energy->getPointer();
         }

         assert( phase->getGhostCellWidth()==hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );
         assert( weight->getGhostCellWidth()==hier::IntVector(tbox::Dimension(NDIM),0) );
#if (NDIM == 3)
         if ( d_model_parameters.with_orientation() )assert( pgrad_quat[2]!=nullptr );
#endif

         const char interpf = energyInterpChar(d_model_parameters.energy_interp_func_type());
         const char interpe = energyInterpChar(d_model_parameters.eta_interp_func_type());

         FORT_QUATENERGY(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            d_qlen,
            dx,
            pgrad_quat[0],
            pgrad_quat[1],
#if (NDIM == 3)
            pgrad_quat[2],
#endif
            0,
            phase->getPointer(), NGHOSTS,
            ptr_eta, NGHOSTS,
            quat_ptr, NGHOSTS,
            d_model_parameters.epsilon_phase(),
            d_model_parameters.epsilon_eta(),
            d_model_parameters.epsilon_q(),
            epsilon_anisotropy, 4,
            2.*d_model_parameters.H_parameter(),
            temperature->getPointer(),
            temperature->getGhostCellWidth()[0],
            d_model_parameters.phase_well_scale(),
            d_model_parameters.eta_well_scale(),
            fl->getPointer(),
            fa->getPointer(),
            ptr_fb,
            three_phase,
            weight->getPointer(),
            total_energy,
            total_phase_e,
            total_eta_e,
            total_orient_e,
            total_qint_e,
            total_well_e,
            total_free_e,
            ptr_energy,
            per_cell,
            &interpf, &interpe,
            d_model_parameters.phase_well_func_type().c_str(),
            d_model_parameters.eta_well_func_type().c_str(),
            d_model_parameters.orient_interp_func_type().c_str(),
            d_model_parameters.avg_func_type().c_str(),
            d_model_parameters.quat_grad_floor_type().c_str(),
            d_model_parameters.quat_grad_floor());
         } // with_phase
      }
   }

   total_energy = sumReduction( total_energy );
   total_phase_e = sumReduction( total_phase_e );
   total_orient_e = sumReduction( total_orient_e );
   total_qint_e = sumReduction( total_qint_e );
   total_well_e = sumReduction( total_well_e );
   total_free_e = sumReduction( total_free_e );
   
   math::HierarchyCellDataOpsReal<double> mathops( hierarchy );
    
   if( d_model_parameters.with_visit_energy_output() )
   {
      double emin = mathops.min( d_energy_diag_id );
      double emax = mathops.max( d_energy_diag_id );
      tbox::plog<<"Min. energy density = "<<emin<<endl;
      tbox::plog<<"Max. energy density = "<<emax<<endl;
   }
   
}

//=======================================================================

void QuatModel::computeSymmetryRotations(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const double time )
{
   (void)time;

   assert( d_quat_scratch_id>=0 );
   assert( d_quat_symm_rotation_id>=0 );
   
   tbox::pout<<"compute symmetry rotations..."<<endl;
   
   // Fill ghosts of original quat data
   copyCurrentToScratch(
      d_patch_hierarchy,
      d_time,
      d_all_refine_patch_strategy );

   int maxln = hierarchy->getFinestLevelNumber();
   for ( int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );

      for ( hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p ) {
         boost::shared_ptr<hier::Patch > patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast =  pbox.upper();

         boost::shared_ptr< pdat::CellData<double> > quat (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_quat_scratch_id) ) );
         assert( quat );
         assert( quat->getGhostCellWidth() ==
                 hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

         boost::shared_ptr< pdat::SideData<int> > rotation_index (
            BOOST_CAST< pdat::SideData<int>, hier::PatchData>(patch->getPatchData( d_quat_symm_rotation_id) ) );
         assert( rotation_index );
         assert( rotation_index->getGhostCellWidth() ==
                 hier::IntVector(tbox::Dimension(NDIM),NGHOSTS) );

         FORT_QUATSYMMROTATION(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            quat->getPointer(),
            quat->getGhostCellWidth()[0],
            d_qlen,
            rotation_index->getPointer(0),
            rotation_index->getPointer(1),
#if (NDIM == 3)
            rotation_index->getPointer(2),
#endif
            rotation_index->getGhostCellWidth()[0]
            );

         if ( d_model_parameters.with_extra_visit_output() ) {
            assert( d_quat_symm_rotation_cell_id>=0 );
            boost::shared_ptr< pdat::CellData<int> > cell_rot_data (
               BOOST_CAST< pdat::CellData<int>, hier::PatchData>(patch->getPatchData( d_quat_symm_rotation_cell_id) ) );
            assert( cell_rot_data );

            pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
            for ( pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i!= iend; ++i ) {
               const pdat::CellIndex ccell = *i;
               const pdat::SideIndex xside(
                  ccell, pdat::SideIndex::X, pdat::SideIndex::Lower );
               const pdat::SideIndex yside(
                  ccell, pdat::SideIndex::Y, pdat::SideIndex::Lower );
#if (NDIM == 3)
               const pdat::SideIndex zside(
                  ccell, pdat::SideIndex::Z, pdat::SideIndex::Lower );
#endif
               (*cell_rot_data)(ccell,0) = (*rotation_index)(xside);
               (*cell_rot_data)(ccell,1) = (*rotation_index)(yside);
#if (NDIM == 3)
               (*cell_rot_data)(ccell,2) = (*rotation_index)(zside);
#endif
            }
         }
      }
   }
}

//=======================================================================

void QuatModel::makeQuatFundamental(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const double time )
{
   (void)time;

   assert( d_quat_id>=0 );
   
   if ( d_verbosity->notSilent() ) {
      tbox::pout << "Setting fundamental orientation" << endl;
   }

   int maxln = hierarchy->getFinestLevelNumber();
   for ( int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );

      for ( hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p ) {
         boost::shared_ptr<hier::Patch > patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast =  pbox.upper();

         boost::shared_ptr< pdat::CellData<double> > quat (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_quat_id) ) );
         assert( quat );
         const hier::Box & quat_gbox = quat->getGhostBox();
         const hier::Index& q_lower = quat_gbox.lower();
         const hier::Index& q_upper = quat_gbox.upper();

         FORT_QUATFUNDAMENTAL(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            quat->getPointer(),
            q_lower[0], q_upper[0],
            q_lower[1], q_upper[1],
#if (NDIM == 3)
            q_lower[2], q_upper[2],
#endif
            d_qlen
            );

      }
   }
}

//=======================================================================

double QuatModel::evaluateIntegralConcentration(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int component )
{
   assert( d_weight_id != -1 );

   if ( ! d_model_parameters.with_concentration() ) return 0.;

   math::HierarchyCellDataOpsReal<double> mathops( hierarchy );

   int conc_id = d_conc_id;
   if( d_ncompositions>1 ){
      assert( d_work_id != -1 );

      copyDepthCellData(hierarchy, d_work_id, 0,
                                   d_conc_id, component);
      conc_id = d_work_id;
   }

   double value = mathops.L1Norm( conc_id, d_weight_id );

   return value;
}

//=======================================================================

double QuatModel::evaluateIntegralPhaseConcentration(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
   const int component )
{
   assert( d_weight_id != -1 );

   if ( ! d_model_parameters.with_concentration() ) return 0.;

   math::HierarchyCellDataOpsReal<double> mathops( hierarchy );

   int conc_id = d_conc_id;
   if( d_ncompositions>1 ){
      copyDepthCellData(hierarchy, d_work_id, 0,
                                   d_conc_id, component);
      conc_id = d_work_id;
   }

   mathops.multiply(d_phase_scratch_id, conc_id, d_phase_id);
        
   double value = mathops.L1Norm( d_phase_scratch_id, d_weight_id );
    
   return value;
}

//=======================================================================

double QuatModel::evaluateVolumeSolid(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy )
{
   assert( d_weight_id != -1 );

   math::HierarchyCellDataOpsReal<double> mathops( hierarchy );
    
   applyPolynomial( hierarchy, d_phase_id, d_phase_scratch_id );
   double value = mathops.L1Norm( d_phase_scratch_id, d_weight_id );
    
   return value;
}

//=======================================================================

double QuatModel::evaluateVolumeEta(
   const boost::shared_ptr< hier::PatchHierarchy > hierarchy )
{
   assert( d_weight_id != -1 );
   assert( d_eta_id    != -1 );

   math::HierarchyCellDataOpsReal<double> mathops( hierarchy );
    
   double value = mathops.L1Norm( d_eta_id, d_weight_id );
    
   return value;
}

//=======================================================================

void QuatModel::fillPartitionCoeffGhosts( void )
{
   assert( d_partition_coeff_id>=0 );
   assert( d_partition_coeff_scratch_id>=0 );
   if ( ! d_all_periodic )
      assert( d_partition_coeff_refine_patch_strategy!=nullptr );
   
   //tbox::pout<<"QuatModel::fillPartitionCoeffGhosts"<<endl;
   
   xfer::RefineAlgorithm copy_to_scratch;

   boost::shared_ptr<hier::RefineOperator > refine_op =
      d_grid_geometry->lookupRefineOperator(
         d_partition_coeff_var,
         "LINEAR_REFINE" );

   copy_to_scratch.registerRefine(
      d_partition_coeff_scratch_id,  // destination
      d_partition_coeff_id,          // source
      d_partition_coeff_scratch_id,  // temporary work space
      refine_op );

   const int maxl = d_patch_hierarchy->getNumberOfLevels();

   for ( int ln = 0; ln < maxl; ln++ ) {
      boost::shared_ptr< hier::PatchLevel > level =
         d_patch_hierarchy->getPatchLevel( ln );

      copy_to_scratch.createSchedule(
         level,
         ln-1,
         d_patch_hierarchy,
         d_partition_coeff_refine_patch_strategy )
         ->fillData( d_time );
   }
}

//=======================================================================

void QuatModel::resetRefPhaseConcentrations()
{
   assert( d_conc_l_scratch_id>=0 );
   assert( d_conc_a_scratch_id>=0 );
   assert( d_conc_l_ref_id>=0 );
   assert( d_conc_a_ref_id>=0 );
   
   //tbox::pout<<"QuatModel::resetRefPhaseConcentrations()"<<endl;
   
   math::HierarchyCellDataOpsReal<double> cellops( d_patch_hierarchy );
   cellops.copyData( d_conc_l_ref_id, d_conc_l_scratch_id, false );
   cellops.copyData( d_conc_a_ref_id, d_conc_a_scratch_id, false );
   if( d_model_parameters.with_third_phase() )
      cellops.copyData( d_conc_b_ref_id, d_conc_b_scratch_id, false );

}

//=======================================================================

void QuatModel::setPhaseConcentrationsToEquilibrium(const double* const ceq)
{
   assert( d_conc_l_id>=0 );
   assert( d_conc_a_id>=0 );
   
   tbox::pout<<"QuatModel::setPhaseConcentrationsToEquilibrium(ceq)"<<endl;
#if 0
   math::HierarchyCellDataOpsReal<double> cellops( d_patch_hierarchy );
   cellops.setToScalar( d_conc_l_id, ceq[0], false );
   cellops.setToScalar( d_conc_a_id, ceq[1], false );
   if( d_model_parameters.with_third_phase() )
      cellops.setToScalar( d_conc_b_id, ceq[2], false );
#else
   int offset=d_ncompositions;
   for (int ln=0; ln<=d_patch_hierarchy->getFinestLevelNumber(); ln++) {
      boost::shared_ptr< hier::PatchLevel > level =
         d_patch_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip!=level->end(); ip++) {
         boost::shared_ptr<hier::Patch > patch = *ip;

         boost::shared_ptr< pdat::CellData<double> > concl (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_conc_l_id ) ) );
         assert( concl );

         for(int i=0;i<d_ncompositions;i++)
            concl->fill(ceq[i],i);

         boost::shared_ptr< pdat::CellData<double> > conca (
            BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_conc_a_id ) ) );
         assert( conca );

         for(int i=0;i<d_ncompositions;i++)
            conca->fill(ceq[offset+i],i);
         if( d_model_parameters.with_third_phase() ){
            boost::shared_ptr< pdat::CellData<double> > concb (
               BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_conc_b_id ) ) );
            assert( concb );

            for(int i=0;i<d_ncompositions;i++)
               concb->fill(ceq[2*offset+i],i);

         }
      }
   }
#endif
}

//=======================================================================

void QuatModel::setRefPhaseConcentrationsToEquilibrium(const double* const ceq)
{
   assert( d_conc_l_ref_id>=0 );
   assert( d_conc_a_ref_id>=0 );
   
   tbox::pout<<"QuatModel::setRefPhaseConcentrationsToEquilibrium(ceq)"<<endl;

   math::HierarchyCellDataOpsReal<double> cellops( d_patch_hierarchy );
   cellops.setToScalar( d_conc_l_ref_id, ceq[0], false );
   cellops.setToScalar( d_conc_a_ref_id, ceq[1], false );
   if( d_model_parameters.with_third_phase() )
      cellops.setToScalar( d_conc_b_ref_id, ceq[2], false );

}

//=======================================================================

//evaluate velocity field at every cell
void QuatModel::computeVelocity(boost::shared_ptr<hier::Patch > patch,
   int phi_dot_id)
{
   assert( d_phase_grad_cell_id>=0 );
   assert( d_velocity_id>=0 );
   assert( phi_dot_id>=0 );

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast =  pbox.upper();

   boost::shared_ptr< pdat::CellData<double> > grad_cell_data (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_phase_grad_cell_id) ) );

   boost::shared_ptr< pdat::CellData<double> > phi_dot_data (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( phi_dot_id) ) );

   boost::shared_ptr< pdat::CellData<double> > velocity_data (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_velocity_id) ) );

   assert( grad_cell_data );
   assert( phi_dot_data );
   assert( velocity_data );
   
   assert( grad_cell_data->getGhostCellWidth()[0]==0 );
   assert( phi_dot_data->getGhostCellWidth()[0]==0 );
   assert( velocity_data->getGhostCellWidth()[0]==0 );

   boost::shared_ptr< geom::CartesianPatchGeometry > patch_geom ( 
      BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(patch->getPatchGeometry()) );
   TBOX_ASSERT(patch_geom);

   const double * dx = patch_geom->getDx();
   const double threshold = 0.02/dx[0];
   //tbox::pout<<"QuatModel::computeVelocity() with threshold "<<threshold<<endl;
   
   FORT_VELOCITY(
         ifirst(0), ilast(0),
         ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         threshold,
         grad_cell_data->getPointer( 0 ),
         grad_cell_data->getPointer( 1 ),
#if (NDIM == 3)
         grad_cell_data->getPointer( 2 ),
#endif
         phi_dot_data->getPointer(),
         velocity_data->getPointer()
         );
}

//=======================================================================

void QuatModel::computeVelocity(const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
                                int phi_dot_id)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for ( int ln = 0; ln <= maxln; ln++ ) {
      boost::shared_ptr<hier::PatchLevel > level =
         hierarchy->getPatchLevel( ln );

      for ( hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); ++p ) {
         boost::shared_ptr<hier::Patch > patch = *p;
         computeVelocity(patch, phi_dot_id);
      }
   }
}

//=======================================================================

double QuatModel::computeThermalEnergy( const boost::shared_ptr<hier::PatchHierarchy > hierarchy )
{
   math::HierarchyCellDataOpsReal<double> cellops( hierarchy );

   double lenergy = 0.;
   if( d_model_parameters.with_phase() ){
      lenergy = cellops.integral(d_phase_id,d_weight_id);
      lenergy *= (-1.*d_model_parameters.latent_heat() );
   }
   //double refenergy=d_model_parameters.rescale_factorT()*cellops.integral(d_cp_id, d_weight_id );
   //if ( d_model_parameters.with_rescaled_temperature() )refenergy/d_model_parameters.rescale_factorT();

   // store product cp*T in d_fl_id
   cellops.multiply( d_f_l_id, d_cp_id, d_temperature_id ); // rescaling of cp compensates rescaling of T

   double cenergy = cellops.integral(d_f_l_id,d_weight_id);

   return lenergy+cenergy;//-refenergy; 
}

