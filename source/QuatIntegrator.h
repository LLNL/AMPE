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
#ifndef included_QuatIntegrator
#define included_QuatIntegrator

#define USE_CPODE

#include "QuatSysSolver.h"
#include "HeatCapacityStrategy.h"
#include "PartitionCoefficientStrategy.h"
#include "TemperatureStrategy.h"
#include "QuatModelParameters.h"
#include "TemperatureFACOps.h"
#include "PhaseFACOps.h"
#include "PhaseTemperatureFACOps.h"
#include "EllipticFACSolver.h"

// Headers for SAMRAI objects
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/IEEE.h"
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideVariable.h"

#ifdef USE_CPODE
#include "CPODESSolver.h"
#include "CPODESAbstractFunctions.h"
#else
#include "SAMRAI/solv/CVODESolver.h"
#include "SAMRAI/solv/CVODEAbstractFunctions.h"
#endif

#ifndef included_solv_SundialsAbstractVector
#include "SAMRAI/solv/SundialsAbstractVector.h"
#endif


#include <boost/make_shared.hpp>
#include <set>
#include <vector>

using namespace SAMRAI;

class QuatRefinePatchStrategy;
class GradStrategy;
class QuatGradStrategy;
class QuatMobilityStrategy;
class FreeEnergyStrategy;
class QuatModel;
class TemperatureStrategy;
class PhaseFACSolver;
class EtaFACSolver;
class ConcFACSolver;
class TemperatureFACSolver;
class CompositionRHSStrategy;
class PhaseFluxStrategy;
class PhaseConcentrationsStrategy;

class QuatIntegrator
   : public mesh::StandardTagAndInitStrategy
#ifdef USE_CPODE
   , public CPODESAbstractFunctions
#else
   , public solv::CVODEAbstractFunctions
#endif
{
public :
   QuatIntegrator(
      const std::string& name,
      const QuatModelParameters& model_parameters,
      QuatModel* model,
      const boost::shared_ptr< hier::VariableContext > current,
      const boost::shared_ptr< hier::VariableContext > scratch,
      const int qlen,
      const int ncompositions,
      boost::shared_ptr<tbox::Database> input_db,
      boost::shared_ptr<geom::CartesianGridGeometry > grid_geom,
      boost::shared_ptr<tbox::Database> bc_db,
      const bool with_phase,
      const bool with_orientation,
      const bool with_concentration,
      const bool with_third_phase,
      const bool with_heat_equation,
      const bool with_steady_temperature,
      const bool with_gradT,
      const bool with_antitrapping,
      const bool with_partition_coeff,
      const bool use_warm_start,
      const bool symmetry_aware,
      const bool use_gradq_for_flux );

   virtual ~QuatIntegrator();

   void initialize(
      const boost::shared_ptr<hier::PatchHierarchy >& hierarchy );

   void initializeCoarseRefineOperators(
      boost::shared_ptr<mesh::GriddingAlgorithm > gridding_alg,
      boost::shared_ptr<hier::RefineOperator > quat_refine_op,
      boost::shared_ptr<hier::CoarsenOperator > quat_coarsen_op );

   void setVerbosity( const int v );
   void setupPreconditioners(boost::shared_ptr<tbox::Database>);

   void setModelParameters(
      const double current_time,
      const double end_time,
      const double h_parameter,
      const double epsilon_phase,
      const double epsilon_eta,
      const double epsilon_q,
      const double quat_grad_floor,
      const std::string quat_smooth_floor_type,
      const double phase_well_scale,
      const double eta_well_scale,
      const std::string orient_interp_func_type,
      const std::string avg_func_type,
      const std::string phase_well_func_type,
      const std::string phase_interp_func_type,
      const std::string eta_well_func_type,
      const std::string eta_interp_func_type );

   void setAbsTol(const double atol)
   {
      d_atol = atol;
   }

   void setRelTol(const double rtol)
   {
      d_rtol = rtol;
   }

   void setConcentrationModelParameters(
      const double mobility);

   void RegisterVariables(
      const boost::shared_ptr< pdat::CellVariable<double> > phase_var,
      const boost::shared_ptr< pdat::CellVariable<double> > eta_var,
      const boost::shared_ptr< pdat::CellVariable<double> > quat_var,
      const boost::shared_ptr< pdat::CellVariable<double> > quat_grad_cell_var,
      const boost::shared_ptr< pdat::SideVariable<double> > quat_grad_side_var,
      const boost::shared_ptr< pdat::CellVariable<double> > quat_grad_modulus_var,
      const boost::shared_ptr< pdat::CellVariable<double> > phase_mobility_var,
      const boost::shared_ptr< pdat::CellVariable<double> > eta_mobility_var,
      const boost::shared_ptr< pdat::CellVariable<double> > quat_mobility_var,
      const boost::shared_ptr< pdat::SideVariable<double> > quat_diffusion_var,
      const boost::shared_ptr< pdat::SideVariable<double> > quat_diffs_var,
      const boost::shared_ptr< pdat::SideVariable<int> >    quat_symm_rotation_var,
      const boost::shared_ptr< pdat::CellVariable<double> > weight_var,
      const boost::shared_ptr< pdat::CellVariable<double> > temperature_var,
      const boost::shared_ptr< pdat::CellVariable<double> > cp_var);

   void RegisterConcentrationVariables(
      const boost::shared_ptr< pdat::CellVariable<double> > conc_var,
      const boost::shared_ptr< pdat::SideVariable<double> > conc_diffusion0_var,
      const boost::shared_ptr< pdat::SideVariable<double> > conc_phase_coupling_diffusion_var,
      const boost::shared_ptr< pdat::SideVariable<double> > conc_eta_coupling_diffusion_var,
      const boost::shared_ptr< pdat::SideVariable<double> > conc_diffusion_var
   );

   void RegisterFreeEnergyVariables(
      const boost::shared_ptr< pdat::CellVariable<double> > f_l_var,
      const boost::shared_ptr< pdat::CellVariable<double> > f_a_var,
      const boost::shared_ptr< pdat::CellVariable<double> > f_b_var );

   void RegisterWithVisit(
      boost::shared_ptr<appu::VisItDataWriter > visit_data_writer);

   double Advance(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const double beginning_time,
      const int cycle );

   void updateSolution(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int coarsest_level,
      const int finest_level );

   void printSolverTotals( void );

   void setTimestep( const double dt );

   /////////////////////////////////////////////////////////////
   //
   // Methods inherited from StandardTagAndInitStrategy
   //

   void initializeLevelData(
      const boost::shared_ptr<hier::PatchHierarchy >& hierarchy,
      const int level_number,
      const double time,
      const bool can_be_refined,
      const bool initial_time,
      const boost::shared_ptr<hier::PatchLevel >& old_level,
      const bool allocate_data );

   void resetHierarchyConfiguration(
      const boost::shared_ptr<hier::PatchHierarchy >& hierarchy,
      const int coarsest_level,
      const int finest_level);
   
   /////////////////////////////////////////////////////////////
   
   /*!
    * compute gradient of of class quat_var 
    * (which may not be the same as quat_var in class QuatModel)
    */
   void computeQuatGradients(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
      double time,
      const bool recompute_quat_sidegrad );

   /**
    * User-supplied right-hand side function evaluation inherited
    * from CPODEAbstractFunctions or CVODEAbstractFunctions.
    *
    * The function arguments are:
    *
    
    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    *
    
    *
    * IMPORTANT: This function must not modify the vector y.
    */
   int evaluateRHSFunction(
      double time,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* y_dot,
      int fd_flag);

#ifdef USE_CPODE
   
   //
   // Methods inherited from CPODEAbstractFunctions
   //
   
   int applyProjection(
      double time,
      solv::SundialsAbstractVector * y,
      solv::SundialsAbstractVector * corr,
      double epsProj,
      solv::SundialsAbstractVector * err);

   int CPSpgmrPrecondSet(
      double t,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* fy,
      int jok,
      int *jcurPtr,
      double gamma,
      solv::SundialsAbstractVector* vtemp1,
      solv::SundialsAbstractVector* vtemp2,
      solv::SundialsAbstractVector* vtemp3);

   int CPSpgmrPrecondSolve(
      double t,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* fy,
      solv::SundialsAbstractVector* r,
      solv::SundialsAbstractVector* z,
      double gamma,
      double delta,
      int lr,
      solv::SundialsAbstractVector* vtemp);
#else

   //
   // Methods inherited from CVODEAbstractFunctions
   //

   int CVSpgmrPrecondSet(
      double t,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* fy,
      int jok,
      int *jcurPtr,
      double gamma,
      solv::SundialsAbstractVector* vtemp1,
      solv::SundialsAbstractVector* vtemp2,
      solv::SundialsAbstractVector* vtemp3);

   int CVSpgmrPrecondSolve(
      double t,
      solv::SundialsAbstractVector* y,
      solv::SundialsAbstractVector* fy,
      solv::SundialsAbstractVector* r,
      solv::SundialsAbstractVector* z,
      double gamma,
      double delta,
      int lr,
      solv::SundialsAbstractVector* vtemp);
#endif

   void setQuatGradStrategy( QuatGradStrategy* quat_grad_strategy );
   void setMobilityStrategy( QuatMobilityStrategy* mobility_strategy );
   void setFreeEnergyStrategy( FreeEnergyStrategy* free_energy_strategy);
   void setPhaseConcentrationsStrategy(
      PhaseConcentrationsStrategy* phase_conc_strategy)
   {
      assert( phase_conc_strategy != NULL ); 
      d_phase_conc_strategy = phase_conc_strategy;
   }
   void setCompositionRHSStrategy( CompositionRHSStrategy* composition_rhs_strategy);
   void setPhaseFluxStrategy( PhaseFluxStrategy* );
   void setTemperatureStrategy( TemperatureStrategy* );
   void setHeatCapacityStrategy( HeatCapacityStrategy* );
   
   void setPartitionCoefficientStrategy(
      PartitionCoefficientStrategy* partition_coeff_strategy)
   {
      assert( partition_coeff_strategy != NULL ); 
      d_partition_coeff_strategy = partition_coeff_strategy;
   }

   void getCPODESIdsRequiringRegrid( 
      std::set< int >& cpode_id_set,
      std::set< int >& phase_id_set,
      std::set< int >& eta_id_set,
      std::set< int >& orient_id_set,
      std::set< int >& conc_id_set,
      std::set< int >& temp_id_set );

   void setupBC();

protected:

   void evaluatePhaseRHS(
      const double time,
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > y_dot_samvect,
      int fd_flag )
   {
      const int ydot_phase_id =
         y_dot_samvect->getComponentDescriptorIndex( d_phase_component_index );
      
      evaluatePhaseRHS(
         time,
         hierarchy,
         d_phase_scratch_id,
         d_eta_scratch_id,
         d_conc_scratch_id,
         d_quat_scratch_id,
         ydot_phase_id,
         d_temperature_scratch_id,
         fd_flag==0 );
   }

   void evaluateQuatRHS(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > y_dot_samvect,
      int fd_flag )
   {
      const int ydot_quat_id =
         y_dot_samvect->getComponentDescriptorIndex( d_quat_component_index );

      evaluateQuatRHS(
         hierarchy,
         d_quat_scratch_id,
         d_phase_scratch_id,
         ydot_quat_id,
         fd_flag==0 );         
   }

   void evaluateTemperatureRHS(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > y_dot_samvect,
      int fd_flag )
   {
      const int ydot_phase_id = 
         d_with_phase ? y_dot_samvect->getComponentDescriptorIndex( d_phase_component_index )
                    : -1;
      const int ydot_temperature_id =
         y_dot_samvect->getComponentDescriptorIndex( d_temperature_component_index );
      evaluateTemperatureRHS(
         hierarchy,
         d_temperature_scratch_id,
         ydot_phase_id,
         ydot_temperature_id,
         fd_flag==0 );
   }
   
   void setTemperatureField(boost::shared_ptr<hier::PatchHierarchy > hierarchy,
                            double time)
   {
      if ( ! d_model_parameters.with_unsteady_heat_equation() )
      {
         assert( d_temperature_strategy!=NULL );
         d_temperature_strategy->setCurrentTemperature( hierarchy, time );
      }
   }

   void evaluateConcentrationRHS(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int phase_id,
      const int conc_rhs_id,
      const int temperature_id,
      const int phase_rhs_id );

   void resetSolutionVector(const boost::shared_ptr<hier::PatchHierarchy > hierarchy);
   void fillScratch(
      double time,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > y );

   void computeMobilities(double time,boost::shared_ptr<hier::PatchHierarchy > hierarchy);

   void setDiffusionCoeffForQuat(
      const boost::shared_ptr< hier::PatchHierarchy >,
      const double time );
   void setUniformDiffusionCoeffForQuat(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy );
   void setDerivDiffusionCoeffForQuat(
      const boost::shared_ptr< hier::PatchHierarchy >,
      const double time );
   
   boost::shared_ptr<solv::SAMRAIVectorReal<double> > d_solution_vec;

   // SAMRAIVector component index
   int d_phase_component_index;
   int d_eta_component_index;
   int d_quat_component_index;
   int d_conc_component_index;
   int d_temperature_component_index;

   //
   // Variables owned by QuatModel
   //

   boost::shared_ptr< pdat::CellVariable<double> > d_phase_var;
   int d_phase_id;
   int d_phase_scratch_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_temperature_var;
   int d_temperature_id;
   int d_temperature_scratch_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_quat_var;
   int d_quat_id;
   int d_quat_scratch_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_conc_var;
   int d_conc_id;
   int d_conc_scratch_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_weight_var;
   int d_weight_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_phase_mobility_var;
   int d_phase_mobility_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_phase_temperature_mobility_var;
   int d_phase_temperature_mobility_id;
   const int d_ncompositions;

   bool d_with_phase;
   const bool d_with_concentration;
   bool d_with_orientation;
   bool d_with_unsteady_temperature;

   bool d_precond_has_dquatdphi;
   bool d_precond_has_dTdphi;
   bool d_precond_has_dPhidT;

   boost::shared_ptr<geom::CartesianGridGeometry > d_grid_geometry;
   
   boost::shared_ptr<hier::CoarsenOperator > d_quat_coarsen_op;
   boost::shared_ptr<hier::CoarsenOperator > d_phase_coarsen_op;
   boost::shared_ptr<hier::CoarsenOperator > d_eta_coarsen_op;
   boost::shared_ptr<hier::CoarsenOperator > d_temperature_coarsen_op;
   boost::shared_ptr<hier::CoarsenOperator > d_conc_coarsen_op;

   boost::shared_ptr<hier::RefineOperator > d_conc_refine_op;

   /*!
       Flag controlling the lagging of the quaternion
       side gradients in a step:

       false = no lagging
       true  = side gradients only recomputed when 
               the integrator computes the system 
	       right-hand side
   */       
   bool d_lag_quat_sidegrad;

   std::string d_phase_interp_func_type;

   FreeEnergyStrategy*   d_free_energy_strategy;

   hier::ComponentSelector d_local_data;

   double d_uniform_diffusion_time_threshold;

   double d_conc_mobility;

   bool d_show_conc_sys_stats;

   solv::LocationIndexRobinBcCoefs* d_conc_bc_coefs;

   // Timers
   boost::shared_ptr<tbox::Timer> t_rhs_timer;
   boost::shared_ptr<tbox::Timer> t_set_coeff_timer;

private :

   //
   // Local methods
   //
   virtual void setCoefficients(
      double time,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > y,
      const bool recompute_quat_sidegrad );

   void createSundialsSolver();
   void computePhaseConcentrations(const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > y);
   virtual void setupPreconditionersConcentration(boost::shared_ptr<tbox::Database> integrator_db);
   void setupPreconditionersPhase(boost::shared_ptr<tbox::Database> integrator_db);
   void setupPreconditionersEta(boost::shared_ptr<tbox::Database> integrator_db);
   void setupPreconditionersTemperature(boost::shared_ptr<tbox::Database> integrator_db);

   void setSolversBoundaries();
   virtual void setConcentrationSolverBoundaries();

   void RegisterLocalPhaseVariables();
   virtual void RegisterLocalConcentrationVariables();
   void RegisterLocalQuatVariables();
   void RegisterLocalUnsteadyTemperatureVariables();
   void RegisterLocalEtaVariables();
   void RegisterLocalVisitVariables();

   virtual int applyPhasePreconditioner(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const double t,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > r_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > ewt_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > z_samvect,
      const double delta, const double gamma);
   int PhasePrecondSolve(boost::shared_ptr<hier::PatchHierarchy > hierarchy,
                         int r_phase_id, int ewt_phase_id, int z_phase_id,
                         const double delta, const double gamma);
   int EtaPrecondSolve(boost::shared_ptr<hier::PatchHierarchy > hierarchy,
                       int r_eta_id, int ewt_eta_id, int z_eta_id, 
                       const double delta);
   virtual int applyTemperaturePreconditioner(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const double t,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > r_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > ewt_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > z_samvect,
      const double delta, const double gamma);
   int TemperaturePrecondSolve(boost::shared_ptr<hier::PatchHierarchy > hierarchy,
                               int r_temperature_id, int ewt_temperature_id, int z_temperature_id, 
                               const double delta, const double gamma);
   int ConcentrationPrecondSolve(boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > r_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > ewt_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > z_samvect,
      const double delta);
   virtual int applyConcentrationPreconditioner(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > r_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > ewt_samvect,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > z_samvect,
      const double delta);
   int QuatPrecondSolve(boost::shared_ptr<hier::PatchHierarchy > hierarchy,
                        int r_quat_id, int ewt_quat_id, int z_quat_id, 
                        const double delta, const double gamma);

   void updateDependentVariables(
      const boost::shared_ptr<hier::PatchHierarchy >& hierarchy,
      const boost::shared_ptr<hier::VariableContext> src_context,
      const boost::shared_ptr<hier::VariableContext> dst_context);

   void updateDependentVariables(
      const boost::shared_ptr<hier::PatchLevel> level,
      const boost::shared_ptr<hier::VariableContext> src_context,
      const boost::shared_ptr<hier::VariableContext> dst_context);

   
   void coarsenData(
      const int phase_id,
      const int eta_id,
      const int quat_id,
      const int conc_id,
      const int temperature_id,
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy );

   void updateSolutionTime(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const double time);

   void updateTimeForPatchID(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int id );

   void setDiffusionCoeffForQuatPatch( hier::Patch& patch );
   void setDerivDiffusionCoeffForQuatPatch( hier::Patch& patch );
   void setDiffusionCoeffForConcentration(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
      const double                                               time );

   void evaluatePhaseRHS(
      const double time,
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int phase_id,
      const int eta_id,
      const int conc_id,
      const int quat_id,
      const int phase_rhs_id,
      const int temperature_id,
      const bool visit_flag );
   void evaluateEtaRHS( 
      const double time,
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int phase_id,
      const int eta_id,
      const int conc_id,
      const int quat_id,
      const int eta_rhs_id,
      const int temperature_id );
   void evaluateTemperatureRHS(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int phase_rhs_id,
      const int temperature_rhs_id,
      const bool visit_flag );
   void evaluateQuatRHS(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int quat_id,
      const int phase_id,
      const int quat_rhs_id,
      const bool visit_flag );

   void correctRhsForSymmetry(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int,
      const int);

   void resetIntegrator(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int coarsest_level,
      const int finest_level );

   void resetAfterRegrid(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int coarsest_level,
      const int finest_level );

   void resetSolversState(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int coarsest_level,
      const int finest_level );
  
   void setSundialsOptions();

   std::vector< boost::shared_ptr< solv::SAMRAIVectorReal<double> > >*
      getCPODESVectorsRequiringRegrid( void );

   bool d_compute_velocity;
   void computeVelocity(const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
                        int phi_dot_id);
   
   void initializeNonPeriodicBC();
   void initializeSolvers(const boost::shared_ptr< hier::PatchHierarchy >& hierarchy );
   
   virtual void createSolutionvector(const boost::shared_ptr<hier::PatchHierarchy > hierarchy);
   virtual void fillScratchComposition( 
      double time,
      boost::shared_ptr< solv::SAMRAIVectorReal<double> > y,
      xfer::RefineAlgorithm& copy_to_scratch);
   virtual void resetSolversStateConcentration(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int coarsest_level,
      const int finest_level );
   virtual void initializeConcentrationNonPeriodicBC();
   virtual void initializeConcentrationSolver(const boost::shared_ptr< hier::PatchHierarchy >&);
   virtual void setCompositionOperatorCoefficients(const double gamma);
   
   std::string d_name;

   const QuatModelParameters& d_model_parameters;

   int d_qlen;
   
   boost::shared_ptr<mesh::GriddingAlgorithm > d_gridding_algorithm;

   int d_verbosity;

   bool d_all_periodic;

   boost::shared_ptr<hier::VariableContext> d_scratch;
   boost::shared_ptr<hier::VariableContext> d_current;
   
   QuatGradStrategy*     d_quat_grad_strategy;
   QuatMobilityStrategy* d_mobility_strategy;
   PhaseConcentrationsStrategy* d_phase_conc_strategy;
   PartitionCoefficientStrategy* d_partition_coeff_strategy;
   
   // use temperature strategy if temperature function of time only,
   // or if solving for steady state
   // Not used if temperature evolves with an ODE
   TemperatureStrategy*  d_temperature_strategy;
   HeatCapacityStrategy* d_heat_capacity_strategy;
   
   
   CompositionRHSStrategy* d_composition_rhs_strategy;

   PhaseFluxStrategy* d_phase_flux_strategy;

   double d_current_time;
   double d_previous_timestep;

   //
   // Variables owned by QuatModel
   //

   boost::shared_ptr< pdat::CellVariable<double> > d_eta_var;
   int d_eta_id;
   int d_eta_scratch_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_cp_var;
   int d_cp_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_quat_grad_cell_var;
   int d_quat_grad_cell_id;
   boost::shared_ptr< pdat::SideVariable<double> > d_quat_grad_side_var;
   int d_quat_grad_side_id;
   boost::shared_ptr< pdat::CellVariable<double> > d_quat_grad_modulus_var;
   int d_quat_grad_modulus_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_eta_mobility_var;
   int d_eta_mobility_id;
   boost::shared_ptr< pdat::CellVariable<double> > d_quat_mobility_var;
   int d_quat_mobility_id;

   boost::shared_ptr< pdat::SideVariable<double> > d_quat_diffusion_var;
   int d_quat_diffusion_id;
   boost::shared_ptr< pdat::SideVariable<double> > d_quat_diffusion_deriv_var;
   int d_quat_diffusion_deriv_id;

   boost::shared_ptr< pdat::SideVariable<int> > d_quat_symm_rotation_var;
   int d_quat_symm_rotation_id;

   /*!
    * diffusion coefficient in preconditioner for composition equation
    */
   boost::shared_ptr< pdat::SideVariable<double> > d_conc_diffusion0_var;
   int d_conc_diffusion0_id;

   boost::shared_ptr< pdat::SideVariable<double> > d_conc_diffusion_var;
   int d_conc_diffusion_id;
   boost::shared_ptr< pdat::SideVariable<double> > d_conc_phase_coupling_diffusion_var;
   int d_conc_phase_coupling_diffusion_id;
   boost::shared_ptr< pdat::SideVariable<double> > d_conc_eta_coupling_diffusion_var;
   int d_conc_eta_coupling_diffusion_id;

   boost::shared_ptr< pdat::SideVariable<double> > d_quat_diffs_var;
   int d_quat_diffs_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_f_l_var;
   int d_f_l_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_f_a_var;
   int d_f_a_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_f_b_var;
   int d_f_b_id;

   //
   // Variables owned locally
   //

   boost::shared_ptr< pdat::CellVariable<double> > d_phase_rhs_visit_var;
   int d_phase_rhs_visit_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_q_rhs_visit_var;
   int d_q_rhs_visit_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_modulus_q_rhs_visit_var;
   int d_modulus_q_rhs_visit_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_temperature_rhs_visit_var;
   int d_temperature_rhs_visit_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_q_rhs1_visit_var;
   int d_q_rhs1_visit_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_phase_sol_var;
   int d_phase_sol_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_phase_rhs_var;
   int d_phase_rhs_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_eta_sol_var;
   int d_eta_sol_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_eta_rhs_var;
   int d_eta_rhs_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_temperature_sol_var;
   int d_temperature_sol_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_temperature_rhs_var;
   int d_temperature_rhs_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_quat_sol_var;
   int d_quat_sol_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_quat_rhs_var;
   int d_quat_rhs_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_u_sol_var;
   int d_u_sol_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_u_rhs_var;
   int d_u_rhs_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_conc_sol_var;
   int d_conc_sol_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_conc_rhs_var;
   int d_conc_rhs_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_quat_mobility_deriv_var;
   int d_quat_mobility_deriv_id;

   boost::shared_ptr< pdat::SideVariable<double> > d_flux_var;
   int d_flux_id;

   boost::shared_ptr< pdat::SideVariable<double> > d_flux_conc_var;
   int d_flux_conc_id;
   
   boost::shared_ptr< pdat::CellVariable<double> > d_velocity_var;
   int d_velocity_id;
   
   /*!
    * Possibly lagged copy of d_quat_grad_side_var
    */
   boost::shared_ptr< pdat::SideVariable<double> > d_quat_grad_side_copy_var;
   int d_quat_grad_side_copy_id;

   boost::shared_ptr< pdat::CellVariable<double> > d_tmp1_var;
   int d_tmp1_id;
   boost::shared_ptr< pdat::CellVariable<double> > d_tmp2_var;
   int d_tmp2_id;
   
   xfer::CoarsenAlgorithm d_quat_diffusion_coarsen;
   std::vector<boost::shared_ptr<xfer::CoarsenSchedule > >
      d_quat_diffusion_coarsen_schedule;

   xfer::CoarsenAlgorithm d_quat_diffusion_deriv_coarsen;
   std::vector<boost::shared_ptr<xfer::CoarsenSchedule > >
      d_quat_diffusion_deriv_coarsen_schedule;

   xfer::CoarsenAlgorithm d_conc_diffusion_coarsen;
   std::vector<boost::shared_ptr<xfer::CoarsenSchedule > >
      d_conc_diffusion_coarsen_schedule;

   xfer::CoarsenAlgorithm d_flux_coarsen_algorithm;
   std::vector<boost::shared_ptr<xfer::CoarsenSchedule > > d_flux_coarsen_schedule;

   xfer::CoarsenAlgorithm d_flux_conc_coarsen_algorithm;
   std::vector<boost::shared_ptr<xfer::CoarsenSchedule > > d_flux_conc_coarsen_schedule;

   xfer::CoarsenAlgorithm d_coarsen_alg;

   // model parameters
   double d_H_parameter;
   double d_epsilon_phase;
   double d_epsilon_eta;
   double d_phase_well_scale;
   double d_eta_well_scale;
   
   double d_thermal_diffusivity;
   double d_latent_heat;
   std::vector<double> d_T_source;

   // norm(grad q)^2 coefficient
   double d_epsilon_q;
   
   // antitrapping coefficient
   double d_alpha_AT;

   double d_quat_grad_floor;

   std::string d_quat_smooth_floor_type;

#ifdef USE_CPODE
   CPODESSolver * d_sundials_solver;
#else
   solv::CVODESolver * d_sundials_solver;
#endif

   boost::shared_ptr< QuatSysSolver > d_quat_sys_solver;

   boost::shared_ptr< PhaseFACOps >    d_phase_fac_ops;
   boost::shared_ptr< PhaseFACSolver > d_phase_sys_solver;

   boost::shared_ptr< EtaFACSolver >   d_eta_sys_solver;
   boost::shared_ptr< ConcFACSolver >  d_conc_sys_solver;

   boost::shared_ptr< TemperatureFACOps >     d_temperature_fac_ops;
   boost::shared_ptr< TemperatureFACSolver >  d_temperature_sys_solver;

   boost::shared_ptr<PhaseTemperatureFACOps> d_phase_temperature_fac_ops;
   boost::shared_ptr< EllipticFACSolver >    d_phase_temperature_sys_solver;

   double d_end_time;

   boost::shared_ptr<tbox::Database> d_boundary_cond_db;

   QuatRefinePatchStrategy* d_all_refine_patch_strategy;
   
   boost::shared_ptr<hier::RefineOperator > d_quat_refine_op;
   boost::shared_ptr<hier::RefineOperator > d_phase_refine_op;
   boost::shared_ptr<hier::RefineOperator > d_eta_refine_op;
   boost::shared_ptr<hier::RefineOperator > d_temperature_refine_op;
   boost::shared_ptr<hier::CoarsenOperator > d_flux_coarsen_op;

   bool d_with_third_phase;
   bool d_with_heat_equation;
   bool d_with_steady_temperature;
   bool d_with_gradT;
   bool d_with_antitrapping;
   bool d_with_partition_coeff;
   bool d_use_warm_start;
   bool d_symmetry_aware;

   bool d_show_integrator_stats;
   bool d_show_solver_stats;
   bool d_show_phase_sys_stats;
   bool d_show_eta_sys_stats;
   bool d_show_quat_sys_stats;
   bool d_show_temperature_sys_stats;

   bool d_use_preconditioner;
   bool d_precondition_quat;

   double d_max_step_size;

   double d_rtol;

   double d_atol;

   int d_max_order;

   int d_max_krylov_dimension;

   int d_max_precond_steps;

   // Cumulative number of Newton iterations
   int d_cum_newton_iter;

   // Cumulative number of linear iterations
   int d_cum_lin_iter;

   // Cumulative number of Newton iteration failures
   int d_cum_newton_fail;

   // Cumulative number of linear iteration failures
   int d_cum_lin_fail;

   // Cumulative number of error test failures
   int d_cum_err_test_fail;

   // Cumulative number of function evaluations
   int d_cum_f_eval;

   // Cumulative number of preconditioner setups
   int d_cum_p_setup;

   // Cumulative number of preconditioner applications
   int d_cum_p_apply;


   std::string d_orient_interp_func_type;
   
   /*!
    * function type to take averages between two cell centered quantities
    * when computing D_q in time-evolution equation for q
    */
   std::string d_avg_func_type;
   std::string d_phase_well_func_type;
   std::string d_eta_well_func_type;
   std::string d_eta_interp_func_type;

   QuatModel* d_quat_model;
   
   bool d_use_gradq_for_flux;

   solv::LocationIndexRobinBcCoefs* d_phase_bc_coefs;
   solv::LocationIndexRobinBcCoefs* d_eta_bc_coefs;
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_coefs;
   solv::LocationIndexRobinBcCoefs* d_quat_bc_coefs;
   
   // Timers
   boost::shared_ptr<tbox::Timer> t_advance_timer;
   boost::shared_ptr<tbox::Timer> t_phase_rhs_timer;
   boost::shared_ptr<tbox::Timer> t_eta_rhs_timer;
   boost::shared_ptr<tbox::Timer> t_conc_rhs_timer;
   boost::shared_ptr<tbox::Timer> t_symm_rhs_timer;
   boost::shared_ptr<tbox::Timer> t_set_diffcoeff_conc_timer;
   boost::shared_ptr<tbox::Timer> t_psolve_setup_timer;
   boost::shared_ptr<tbox::Timer> t_psolve_solve_timer;
};

#endif
