// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_QuatIntegrator
#define included_QuatIntegrator

#include "QuatSysSolver.h"
#include "HeatCapacityStrategy.h"
#include "PartitionCoefficientStrategy.h"
#include "TemperatureStrategy.h"
#include "QuatModelParameters.h"
#include "TemperatureFACOps.h"
#include "PhaseFACOps.h"
#include "PhaseTemperatureFACOps.h"
#include "EllipticFACSolver.h"
#include "InterpolationType.h"
#include "DerivDiffusionCoeffForQuat.h"
#include "PhaseRHSStrategy.h"
#include "TemperatureRHSStrategy.h"
#include "MovingFrameRHS.h"
#include "AdaptMovingFrame.h"
#include "QuatFaceCoeff.h"
#include "RigidBodyMotionRHS.h"
#include "RigidBodyMotionConcRHS.h"
#include "RigidBodyForces.h"

// Headers for SAMRAI objects
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/IEEE.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/solv/CartesianRobinBcHelper.h"
#include "SAMRAI/solv/Sundials_SAMRAIVector.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"


#include "CVODESolver.h"
#include "CVODEAbstractFunctions.h"
#include "ARKODESolver.h"
#include "ARKODEAbstractFunctions.h"

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

class QuatIntegrator : public mesh::StandardTagAndInitStrategy,
                       public CVODEAbstractFunctions,
                       public ARKODEAbstractFunctions
{
 public:
   QuatIntegrator(const std::string& name,
                  const QuatModelParameters& model_parameters, QuatModel* model,
                  const std::shared_ptr<hier::VariableContext> current,
                  const std::shared_ptr<hier::VariableContext> scratch,
                  const int qlen, const int ncompositions,
                  std::shared_ptr<tbox::Database> input_db,
                  std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
                  std::shared_ptr<tbox::Database> bc_db, const bool with_phase,
                  const bool with_concentration, const bool with_heat_equation,
                  const bool with_steady_temperature, const bool with_gradT,
                  const bool with_antitrapping, const bool with_partition_coeff,
                  const bool use_warm_start, const bool symmetry_aware,
                  const bool use_gradq_for_flux);

   virtual ~QuatIntegrator();

   void initialize(const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   void initializeCoarseRefineOperators(
       std::shared_ptr<mesh::GriddingAlgorithm> gridding_alg,
       std::shared_ptr<hier::RefineOperator> quat_refine_op,
       std::shared_ptr<hier::CoarsenOperator> quat_coarsen_op);

   void setVerbosity(const int v);
   void setupPreconditioners();

   void setModelParameters(
       const double current_time, const double end_time,
       const double h_parameter, const double epsilon_phase,
       const double epsilon_eta, const double epsilon_q,
       const double quat_grad_floor, const std::string quat_smooth_floor_type,
       const double phase_well_scale, const double eta_well_scale,
       const std::string orient_interp_func_type1,
       const std::string orient_interp_func_type2,
       const std::string avg_func_type, const std::string phase_well_func_type,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       const std::string eta_well_func_type);

   void setAbsTol(const double atol) { d_atol = atol; }

   void setRelTol(const double rtol) { d_rtol = rtol; }

   void setConcentrationModelParameters(const double mobility);

   void RegisterVariables(
       const std::shared_ptr<pdat::CellVariable<double>> phase_var,
       const std::shared_ptr<pdat::CellVariable<double>> eta_var,
       const std::shared_ptr<pdat::CellVariable<double>> quat_var,
       const std::shared_ptr<pdat::CellVariable<double>> quat_grad_cell_var,
       const std::shared_ptr<pdat::SideVariable<double>> quat_grad_side_var,
       const std::shared_ptr<pdat::CellVariable<double>> quat_grad_modulus_var,
       const std::shared_ptr<pdat::CellVariable<double>> phase_mobility_var,
       const std::shared_ptr<pdat::CellVariable<double>> eta_mobility_var,
       const std::shared_ptr<pdat::CellVariable<double>> quat_mobility_var,
       const std::shared_ptr<pdat::SideVariable<double>> quat_diffs_var,
       const std::shared_ptr<pdat::SideVariable<int>> quat_symm_rotation_var,
       const std::shared_ptr<pdat::CellVariable<double>> weight_var,
       const std::shared_ptr<pdat::CellVariable<double>> temperature_var,
       const std::shared_ptr<pdat::CellVariable<double>> cp_var);

   void RegisterConcentrationVariables(
       const std::shared_ptr<pdat::CellVariable<double>> conc_var,
       const std::vector<std::shared_ptr<pdat::SideVariable<double>>>
           conc_pfm_diffusion_var,
       const std::shared_ptr<pdat::SideVariable<double>>
           conc_phase_coupling_diffusion_var,
       const std::shared_ptr<pdat::SideVariable<double>>
           conc_eta_coupling_diffusion_var,
       const std::shared_ptr<pdat::SideVariable<double>> conc_diffusion_var);

   void RegisterFreeEnergyVariables(
       const std::shared_ptr<pdat::CellVariable<double>> f_l_var,
       const std::shared_ptr<pdat::CellVariable<double>> f_a_var,
       const std::shared_ptr<pdat::CellVariable<double>> f_b_var);

   void RegisterWithVisit(
       std::shared_ptr<appu::VisItDataWriter> visit_data_writer);

   double Advance(const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void updateSolution(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                       const int coarsest_level, const int finest_level);

   void printSolverTotals(void);

   void setTimestep(const double dt);

   /////////////////////////////////////////////////////////////
   //
   // Methods inherited from StandardTagAndInitStrategy
   //

   void initializeLevelData(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int level_number, const double time, const bool can_be_refined,
       const bool initial_time,
       const std::shared_ptr<hier::PatchLevel>& old_level,
       const bool allocate_data);

   void resetHierarchyConfiguration(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int coarsest_level, const int finest_level);

   /////////////////////////////////////////////////////////////

   /*!
    * compute gradient of of class quat_var
    * (which may not be the same as quat_var in class QuatModel)
    */
   void computeQuatGradients(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, double time,
       const bool recompute_quat_sidegrad);

   /**
    * User-supplied right-hand side function evaluation inherited
    * from CVODEAbstractFunctions.
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
   int evaluateRHSFunction(double time, solv::SundialsAbstractVector* y,
                           solv::SundialsAbstractVector* y_dot, int fd_flag);

   //
   // Methods inherited from CVODEAbstractFunctions
   //
   int applyProjection(double time, solv::SundialsAbstractVector* y,
                       solv::SundialsAbstractVector* corr, double epsProj,
                       solv::SundialsAbstractVector* err);
   int evaluateJTimesRHSFunction(double t, solv::SundialsAbstractVector* y,
                                 solv::SundialsAbstractVector* y_dot)
   {
      return evaluateRHSFunction(t, y, y_dot, 1);
   }
   int evaluateRHSFunction(double t, solv::SundialsAbstractVector* y,
                           solv::SundialsAbstractVector* y_dot)
   {
      return evaluateRHSFunction(t, y, y_dot, 0);
   }

   int CVSpgmrPrecondSet(double t, solv::SundialsAbstractVector* y,
                         solv::SundialsAbstractVector* fy, int jok,
                         int* jcurPtr, double gamma);

   int CVSpgmrPrecondSolve(double t, solv::SundialsAbstractVector* y,
                           solv::SundialsAbstractVector* fy,
                           solv::SundialsAbstractVector* r,
                           solv::SundialsAbstractVector* z, double gamma,
                           double delta, int lr);

   //
   // Methods inherited from ARKODEAbstractFunctions
   //
   int evaluateRHSFunctionImp(double time, solv::SundialsAbstractVector* y,
                              solv::SundialsAbstractVector* y_dot)
   {
      return 1;
   }

   int evaluateRHSFunctionExp(double time, solv::SundialsAbstractVector* y,
                              solv::SundialsAbstractVector* y_dot)
   {
      return evaluateRHSFunction(time, y, y_dot, 0);
   }

   int ARKSpgmrPrecondSet(double t, solv::SundialsAbstractVector* y,
                          solv::SundialsAbstractVector* fy, int jok,
                          int* jcurPtr, double gamma)
   {
      return CVSpgmrPrecondSet(t, y, fy, jok, jcurPtr, gamma);
   }

   int ARKSpgmrPrecondSolve(double t, solv::SundialsAbstractVector* y,
                            solv::SundialsAbstractVector* fy,
                            solv::SundialsAbstractVector* r,
                            solv::SundialsAbstractVector* z, double gamma,
                            double delta, int lr)
   {
      return CVSpgmrPrecondSolve(t, y, fy, r, z, gamma, delta, lr);
   }

   /*
    * return time of last accepted step
    */
   double lastAcceptedTime() const
   {
      return (
          d_sundials_solver
              ? d_sundials_solver->getActualFinalValueOfIndependentVariable()
              : d_arkode_solver->getActualFinalValueOfIndependentVariable());
   }

   void setQuatGradStrategy(QuatGradStrategy* quat_grad_strategy);
   void setMobilityStrategy(std::shared_ptr<QuatMobilityStrategy>&);
   void setFreeEnergyStrategy(
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy);
   void setCompositionRHSStrategy(
       std::shared_ptr<CompositionRHSStrategy> composition_rhs_strategy);
   void setCompositionDiffusionStrategy(
       std::shared_ptr<CompositionDiffusionStrategy>);

   void setPhaseFluxStrategy(std::shared_ptr<PhaseFluxStrategy>);
   void setTemperatureStrategy(std::shared_ptr<TemperatureStrategy>);
   void setHeatCapacityStrategy(HeatCapacityStrategy*);

   void setRigidBodyForces(std::shared_ptr<RigidBodyForces> rigid_body_forces)
   {
      assert(rigid_body_forces);
      d_rigid_body_forces = rigid_body_forces;
   };

   void setPartitionCoefficientStrategy(
       PartitionCoefficientStrategy* partition_coeff_strategy)
   {
      assert(partition_coeff_strategy != nullptr);
      d_partition_coeff_strategy = partition_coeff_strategy;
   }
#ifdef USE_CPODE
   void getCPODESIdsRequiringRegrid(std::set<int>& cpode_id_set,
                                    std::set<int>& phase_id_set,
                                    std::set<int>& eta_id_set,
                                    std::set<int>& orient_id_set,
                                    std::set<int>& conc_id_set,
                                    std::set<int>& temp_id_set);
#endif
   void setupBC();

   bool needDphiDt() const
   {
      return (d_with_phase &&
              (d_with_antitrapping || d_model_parameters.inMovingFrame() ||
               d_with_heat_equation));
   }

   void fillDphiDt(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                   const double time, const int phase_rhs_id);

 protected:
   void evaluatePhaseRHS(
       const double time, std::shared_ptr<hier::PatchHierarchy> hierarchy,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> y_dot_samvect,
       int fd_flag)
   {
      const int ydot_phase_id =
          y_dot_samvect->getComponentDescriptorIndex(d_phase_component_index);

      evaluatePhaseRHS(time, hierarchy, d_phase_scratch_id, ydot_phase_id,
                       fd_flag == 0);
   }

   void evaluateQuatRHS(
       std::shared_ptr<hier::PatchHierarchy> hierarchy,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> y_dot_samvect,
       int fd_flag)
   {
      const int ydot_quat_id =
          y_dot_samvect->getComponentDescriptorIndex(d_quat_component_index);

      evaluateQuatRHS(hierarchy, ydot_quat_id, fd_flag == 0);
   }

   void evaluateTemperatureRHS(
       std::shared_ptr<hier::PatchHierarchy> hierarchy,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> y_dot_samvect,
       int fd_flag)
   {
      const int phase_rhs_id = d_with_phase ? d_dphidt_scratch_id : -1;
      const int ydot_temperature_id =
          y_dot_samvect->getComponentDescriptorIndex(
              d_temperature_component_index);
      evaluateTemperatureRHS(hierarchy, d_temperature_scratch_id, phase_rhs_id,
                             ydot_temperature_id, fd_flag == 0);
   }

   void setTemperatureField(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                            double time)
   {
      if (!d_model_parameters.with_unsteady_heat_equation()) {
         assert(d_temperature_strategy);
         d_temperature_strategy->setCurrentTemperature(hierarchy, time);
      }
   }

   void evaluateConcentrationRHS(
       std::shared_ptr<hier::PatchHierarchy> hierarchy, const int phase_id,
       const int conc_id, const int conc_rhs_id, const int temperature_id,
       const bool visit_flag);

   void resetSolutionVector(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);
   void fillScratch(double time,
                    std::shared_ptr<solv::SAMRAIVectorReal<double>> y);

   void computeMobilities(double time,
                          std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void setUniformDiffusionCoeffForQuat(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   std::shared_ptr<solv::SAMRAIVectorReal<double>> d_solution_vec;

   // SAMRAIVector component index
   int d_phase_component_index;
   int d_eta_component_index;
   int d_quat_component_index;
   int d_conc_component_index;
   int d_temperature_component_index;

   //
   // Variables owned by QuatModel
   //

   std::shared_ptr<pdat::CellVariable<double>> d_phase_var;
   int d_phase_id;
   int d_phase_scratch_id;

   std::shared_ptr<pdat::CellVariable<double>> d_temperature_var;
   int d_temperature_id;
   int d_temperature_scratch_id;

   std::shared_ptr<pdat::CellVariable<double>> d_quat_var;
   int d_quat_id;
   int d_quat_scratch_id;

   std::shared_ptr<pdat::CellVariable<double>> d_conc_var;
   int d_conc_id;
   int d_conc_scratch_id;

   std::shared_ptr<pdat::CellVariable<double>> d_weight_var;
   int d_weight_id;

   std::shared_ptr<pdat::CellVariable<double>> d_phase_mobility_var;
   int d_phase_mobility_id;

   std::shared_ptr<pdat::CellVariable<double>> d_phase_temperature_mobility_var;
   int d_phase_temperature_mobility_id;
   const int d_ncompositions;

   bool d_with_phase;
   const bool d_with_concentration;
   const bool d_with_orientation;
   const bool d_evolve_quat;
   bool d_with_unsteady_temperature;

   bool d_precond_has_dquatdphi;
   bool d_precond_has_dTdphi;
   bool d_precond_has_dPhidT;

   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   std::shared_ptr<hier::CoarsenOperator> d_quat_coarsen_op;
   std::shared_ptr<hier::CoarsenOperator> d_phase_coarsen_op;
   std::shared_ptr<hier::CoarsenOperator> d_eta_coarsen_op;
   std::shared_ptr<hier::CoarsenOperator> d_temperature_coarsen_op;
   std::shared_ptr<hier::CoarsenOperator> d_conc_coarsen_op;

   std::shared_ptr<hier::RefineOperator> d_conc_refine_op;

   /*!
       Flag controlling the lagging of the quaternion
       side gradients in a step:

       false = no lagging
       true  = side gradients only recomputed when
               the integrator computes the system
               right-hand side
   */
   bool d_lag_quat_sidegrad;

   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type;

   Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;

   std::shared_ptr<FreeEnergyStrategy> d_free_energy_strategy;

   hier::ComponentSelector d_local_data;

   double d_uniform_diffusion_time_threshold;

   double d_conc_mobility;

   bool d_show_conc_sys_stats;

   solv::LocationIndexRobinBcCoefs* d_conc_bc_coefs;
   solv::LocationIndexRobinBcCoefs* d_dphidt_bc_coefs;
   solv::CartesianRobinBcHelper* d_dphidt_bc_helper;

   /*!
    * diffusion coefficient in preconditioner for composition equation
    */
   std::vector<std::shared_ptr<pdat::SideVariable<double>>>
       d_conc_pfm_diffusion_var;
   std::vector<int> d_conc_pfm_diffusion_id;

   // Timers
   std::shared_ptr<tbox::Timer> t_rhs_timer;
   std::shared_ptr<tbox::Timer> t_set_coeff_timer;

 private:
   double currentTime()
   {
      return d_sundials_solver
                 ? d_sundials_solver->getActualFinalValueOfIndependentVariable()
                 : d_arkode_solver->getActualFinalValueOfIndependentVariable();
   }

   double currentDT()
   {
      return d_sundials_solver
                 ? d_sundials_solver->getStepSizeForLastInternalStep()
                 : d_arkode_solver->getStepSizeForLastInternalStep();
   }

   int numberRHSeval()
   {
      return d_sundials_solver
                 ? d_sundials_solver->getNumberOfRHSFunctionCalls()
                 : d_arkode_solver->getNumberOfRHSFunctionExCalls();
   }

   virtual void setCoefficients(
       double time, std::shared_ptr<solv::SAMRAIVectorReal<double>> y,
       const bool recompute_quat_sidegrad);

   void createSundialsSolver();
   virtual void setupPreconditionersConcentration(
       std::shared_ptr<tbox::Database> integrator_db);
   void setupPreconditionersPhase(
       std::shared_ptr<tbox::Database> integrator_db);
   void setupPreconditionersEta(std::shared_ptr<tbox::Database> integrator_db);
   void setupPreconditionersTemperature(
       std::shared_ptr<tbox::Database> integrator_db);

   void setSolversBoundaries();
   virtual void setConcentrationSolverBoundaries();

   void RegisterLocalPhaseVariables();
   virtual void RegisterLocalConcentrationVariables();
   void RegisterLocalQuatVariables();
   void RegisterLocalUnsteadyTemperatureVariables();
   void RegisterLocalEtaVariables();
   void RegisterLocalVisitVariables();

   void RegisterQuatVariables(
       const std::shared_ptr<pdat::CellVariable<double>> quat_var,
       const std::shared_ptr<pdat::CellVariable<double>> quat_grad_cell_var,
       const std::shared_ptr<pdat::SideVariable<double>> quat_grad_side_var,
       const std::shared_ptr<pdat::CellVariable<double>> quat_grad_modulus_var,
       const std::shared_ptr<pdat::CellVariable<double>> quat_mobility_var,
       const std::shared_ptr<pdat::SideVariable<double>> quat_diffs_var,
       const std::shared_ptr<pdat::SideVariable<int>> quat_symm_rotation_var);

   virtual int applyPhasePreconditioner(
       std::shared_ptr<hier::PatchHierarchy> hierarchy, const double t,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> r_samvect,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> z_samvect,
       const double delta, const double gamma);
   int PhasePrecondSolve(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                         int r_phase_id, int z_phase_id, const double delta,
                         const double gamma);
   int EtaPrecondSolve(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                       int r_eta_id, int z_eta_id, const double delta);
   virtual int applyTemperaturePreconditioner(
       std::shared_ptr<hier::PatchHierarchy> hierarchy, const double t,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> r_samvect,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> z_samvect,
       const double delta, const double gamma);
   int TemperaturePrecondSolve(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                               int r_temperature_id, int z_temperature_id,
                               const double delta, const double gamma);
   int ConcentrationPrecondSolve(
       std::shared_ptr<hier::PatchHierarchy> hierarchy,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> r_samvect,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> z_samvect,
       const double delta);
   virtual int applyConcentrationPreconditioner(
       std::shared_ptr<hier::PatchHierarchy> hierarchy,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> r_samvect,
       std::shared_ptr<solv::SAMRAIVectorReal<double>> z_samvect,
       const double delta);
   int QuatPrecondSolve(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int r_quat_id, int z_quat_id, const double delta,
                        const double gamma);

   void updateDependentVariables(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const std::shared_ptr<hier::VariableContext> src_context,
       const std::shared_ptr<hier::VariableContext> dst_context);

   void updateDependentVariables(
       const std::shared_ptr<hier::PatchLevel> level,
       const std::shared_ptr<hier::VariableContext> src_context,
       const std::shared_ptr<hier::VariableContext> dst_context);


   void coarsenData(const int phase_id, const int eta_id, const int quat_id,
                    const int conc_id, const int temperature_id,
                    const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void updateSolutionTime(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

   void updateTimeForPatchID(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int id);

   void setDerivDiffusionCoeffForQuatPatch(hier::Patch& patch);
   void setDiffusionCoeffForConcentration(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

   void evaluatePhaseRHS(const double time,
                         std::shared_ptr<hier::PatchHierarchy> hierarchy,
                         const int phase_id, const int phase_rhs_id,
                         const bool eval_flag);
   void evaluateTemperatureRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                               const int temperature_id, const int phase_rhs_id,
                               const int temperature_rhs_id,
                               const bool visit_flag);
   void evaluateQuatRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        const int quat_rhs_id, const bool visit_flag);

   void correctRhsForSymmetry(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int,
       const int);

   void resetIntegrator(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        const int coarsest_level, const int finest_level);

   void resetAfterRegrid(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                         const int coarsest_level, const int finest_level);

   void resetSolversState(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          const int coarsest_level, const int finest_level);

   void setSundialsOptions();
   void setARKODEOptions();

#ifdef USE_CPODE
   std::vector<std::shared_ptr<solv::SAMRAIVectorReal<double>>>*
   getCPODESVectorsRequiringRegrid(void);
#endif

   void computeVelocity(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int phi_dot_id);

   void initializeNonPeriodicBC();
   void initializeSolvers(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy);

   virtual void createSolutionvector(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);
   virtual void fillScratchComposition(
       double time, std::shared_ptr<solv::SAMRAIVectorReal<double>> y,
       xfer::RefineAlgorithm& copy_to_scratch);
   virtual void resetSolversStateConcentration(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);
   virtual void initializeConcentrationNonPeriodicBC();
   virtual void initializeConcentrationSolver(
       const std::shared_ptr<hier::PatchHierarchy>&);
   virtual void setCompositionOperatorCoefficients(const double gamma);

   std::string d_name;

   const QuatModelParameters& d_model_parameters;

   int d_qlen;

   std::shared_ptr<mesh::GriddingAlgorithm> d_gridding_algorithm;

   int d_verbosity;

   bool d_all_periodic;

   bool d_compute_velocity;

   std::shared_ptr<hier::VariableContext> d_current;
   std::shared_ptr<hier::VariableContext> d_scratch;

   QuatGradStrategy* d_quat_grad_strategy;
   std::shared_ptr<QuatMobilityStrategy> d_mobility_strategy;
   PartitionCoefficientStrategy* d_partition_coeff_strategy;

   // use temperature strategy if temperature function of time only,
   // or if solving for steady state
   // Not used if temperature evolves with an ODE
   std::shared_ptr<TemperatureStrategy> d_temperature_strategy;
   HeatCapacityStrategy* d_heat_capacity_strategy;


   std::shared_ptr<CompositionRHSStrategy> d_composition_rhs_strategy;

   std::shared_ptr<CompositionDiffusionStrategy>
       d_composition_diffusion_strategy;

   std::shared_ptr<PhaseFluxStrategy> d_phase_flux_strategy;

   std::shared_ptr<PhaseRHSStrategy> d_phase_rhs_strategy;

   std::shared_ptr<TemperatureRHSStrategy> d_temperature_rhs_strategy;

   std::shared_ptr<MovingFrameRHS> d_movingframe_phi;
   std::shared_ptr<MovingFrameRHS> d_movingframe_conc;
   std::shared_ptr<MovingFrameRHS> d_movingframe_temp;

   std::shared_ptr<RigidBodyMotionRHS> d_rb_motion;
   std::shared_ptr<RigidBodyMotionConcRHS> d_rb_motion_conc;

   /*!
    * NDIM-dimensional forces acting on each grain
    */
   std::vector<std::array<double, NDIM>> d_rb_forces;

   double d_current_time;
   double d_previous_timestep;

   //
   // Variables owned by QuatModel
   //

   std::shared_ptr<pdat::CellVariable<double>> d_eta_var;
   int d_eta_id;
   int d_eta_scratch_id;

   std::shared_ptr<pdat::CellVariable<double>> d_cp_var;
   int d_cp_id;

   std::shared_ptr<pdat::CellVariable<double>> d_quat_grad_cell_var;
   int d_quat_grad_cell_id;
   std::shared_ptr<pdat::SideVariable<double>> d_quat_grad_side_var;
   int d_quat_grad_side_id;
   std::shared_ptr<pdat::CellVariable<double>> d_quat_grad_modulus_var;
   int d_quat_grad_modulus_id;

   std::shared_ptr<pdat::CellVariable<double>> d_eta_mobility_var;
   int d_eta_mobility_id;
   std::shared_ptr<pdat::CellVariable<double>> d_quat_mobility_var;
   int d_quat_mobility_id;

   std::shared_ptr<pdat::SideVariable<double>> d_quat_diffusion_deriv_var;
   int d_quat_diffusion_deriv_id;

   std::shared_ptr<pdat::SideVariable<int>> d_quat_symm_rotation_var;
   int d_quat_symm_rotation_id;

   std::shared_ptr<pdat::SideVariable<double>> d_conc_diffusion_var;
   int d_conc_diffusion_id;
   std::shared_ptr<pdat::SideVariable<double>>
       d_conc_phase_coupling_diffusion_var;
   int d_conc_phase_coupling_diffusion_id;
   std::shared_ptr<pdat::SideVariable<double>>
       d_conc_eta_coupling_diffusion_var;
   int d_conc_eta_coupling_diffusion_id;

   std::shared_ptr<pdat::SideVariable<double>> d_quat_diffs_var;
   int d_quat_diffs_id;

   std::shared_ptr<pdat::CellVariable<double>> d_f_l_var;
   int d_f_l_id;

   std::shared_ptr<pdat::CellVariable<double>> d_f_a_var;
   int d_f_a_id;

   std::shared_ptr<pdat::CellVariable<double>> d_f_b_var;
   int d_f_b_id;

   //
   // Variables owned locally
   //

   std::shared_ptr<pdat::CellVariable<double>> d_phase_rhs_visit_var;
   int d_phase_rhs_visit_id;

   std::shared_ptr<pdat::CellVariable<double>> d_driving_force_visit_var;
   int d_driving_force_visit_id;

   std::shared_ptr<pdat::CellVariable<double>> d_q_rhs_visit_var;
   int d_q_rhs_visit_id;

   std::shared_ptr<pdat::CellVariable<double>> d_modulus_q_rhs_visit_var;
   int d_modulus_q_rhs_visit_id;

   std::shared_ptr<pdat::CellVariable<double>> d_temperature_rhs_visit_var;
   int d_temperature_rhs_visit_id;

   std::shared_ptr<pdat::CellVariable<double>> d_conc_rhs_visit_var;
   int d_conc_rhs_visit_id;

   std::shared_ptr<pdat::CellVariable<double>> d_q_rhs1_visit_var;
   int d_q_rhs1_visit_id;

   std::shared_ptr<pdat::CellVariable<double>> d_phase_sol_var;
   int d_phase_sol_id;

   std::shared_ptr<pdat::CellVariable<double>> d_phase_rhs_var;
   int d_phase_rhs_id;

   std::shared_ptr<pdat::CellVariable<double>> d_eta_sol_var;
   int d_eta_sol_id;

   std::shared_ptr<pdat::CellVariable<double>> d_eta_rhs_var;
   int d_eta_rhs_id;

   std::shared_ptr<pdat::CellVariable<double>> d_temperature_sol_var;
   int d_temperature_sol_id;

   std::shared_ptr<pdat::CellVariable<double>> d_temperature_rhs_var;
   int d_temperature_rhs_id;

   std::shared_ptr<pdat::CellVariable<double>> d_quat_sol_var;
   int d_quat_sol_id;

   std::shared_ptr<pdat::CellVariable<double>> d_quat_rhs_var;
   int d_quat_rhs_id;

   std::shared_ptr<pdat::CellVariable<double>> d_u_sol_var;
   int d_u_sol_id;

   std::shared_ptr<pdat::CellVariable<double>> d_u_rhs_var;
   int d_u_rhs_id;

   std::shared_ptr<pdat::CellVariable<double>> d_conc_sol_var;
   int d_conc_sol_id;

   std::shared_ptr<pdat::CellVariable<double>> d_conc_rhs_var;
   int d_conc_rhs_id;

   std::shared_ptr<pdat::CellVariable<double>> d_dphidt_scratch_var;
   int d_dphidt_scratch_id;

   std::shared_ptr<pdat::CellVariable<double>> d_quat_mobility_deriv_var;
   int d_quat_mobility_deriv_id;

   std::shared_ptr<pdat::SideVariable<double>> d_flux_var;
   int d_flux_id;

   std::shared_ptr<pdat::SideVariable<double>> d_flux_conc_var;
   int d_flux_conc_id;

   std::shared_ptr<pdat::CellVariable<double>> d_velocity_var;
   int d_velocity_id;

   /*!
    * Possibly lagged copy of d_quat_grad_side_var
    */
   std::shared_ptr<pdat::SideVariable<double>> d_quat_grad_side_copy_var;
   int d_quat_grad_side_copy_id;

   std::shared_ptr<pdat::CellVariable<double>> d_noise_var;
   int d_noise_id;

   std::shared_ptr<pdat::CellVariable<double>> d_tmp1_var;
   int d_tmp1_id;
   std::shared_ptr<pdat::CellVariable<double>> d_tmp2_var;
   int d_tmp2_id;

   std::shared_ptr<DerivDiffusionCoeffForQuat> d_diffusion4quatderiv;

   xfer::CoarsenAlgorithm d_conc_diffusion_coarsen;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule>>
       d_conc_diffusion_coarsen_schedule;

   xfer::CoarsenAlgorithm d_flux_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule>> d_flux_coarsen_schedule;

   xfer::CoarsenAlgorithm d_flux_conc_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule>>
       d_flux_conc_coarsen_schedule;

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

   std::shared_ptr<CVODESolver> d_sundials_solver;
   std::shared_ptr<ARKODESolver> d_arkode_solver;

   std::shared_ptr<QuatSysSolver> d_quat_sys_solver;
   std::shared_ptr<QuatFaceCoeff> d_quat_face_coeff_strategy;

   std::shared_ptr<PhaseFACOps> d_phase_fac_ops;
   std::shared_ptr<PhaseFACSolver> d_phase_sys_solver;

   std::shared_ptr<EtaFACSolver> d_eta_sys_solver;
   std::shared_ptr<ConcFACSolver> d_conc_sys_solver;

   std::shared_ptr<TemperatureFACOps> d_temperature_fac_ops;
   std::shared_ptr<TemperatureFACSolver> d_temperature_sys_solver;

   std::shared_ptr<PhaseTemperatureFACOps> d_phase_temperature_fac_ops;
   std::shared_ptr<EllipticFACSolver> d_phase_temperature_sys_solver;

   double d_end_time;

   std::shared_ptr<tbox::Database> d_boundary_cond_db;

   QuatRefinePatchStrategy* d_all_refine_patch_strategy;

   std::shared_ptr<hier::RefineOperator> d_quat_refine_op;
   std::shared_ptr<hier::RefineOperator> d_phase_refine_op;
   std::shared_ptr<hier::RefineOperator> d_eta_refine_op;
   std::shared_ptr<hier::RefineOperator> d_temperature_refine_op;
   std::shared_ptr<hier::CoarsenOperator> d_flux_coarsen_op;

   bool d_with_third_phase;
   bool d_with_heat_equation;
   bool d_with_steady_temperature;
   bool d_with_gradT;
   bool d_with_antitrapping;
   bool d_with_partition_coeff;
   bool d_use_warm_start;
   bool d_use_cvode;

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

   // interpolation polynomial types for 1st and 2nd order grad q terms
   std::string d_orient_interp_func_type1;
   std::string d_orient_interp_func_type2;

   /*!
    * function type to take averages between two cell centered quantities
    * when computing D_q in time-evolution equation for q
    */
   std::string d_avg_func_type;
   std::string d_phase_well_func_type;
   std::string d_eta_well_func_type;

   QuatModel* d_quat_model;

   std::shared_ptr<AdaptMovingFrame> d_adapt_moving_frame;

   bool d_use_gradq_for_flux;

   /*
    * Velocity of moving frame
    */
   double d_frame_velocity;

   std::shared_ptr<RigidBodyForces> d_rigid_body_forces;

   solv::LocationIndexRobinBcCoefs* d_phase_bc_coefs;
   solv::LocationIndexRobinBcCoefs* d_eta_bc_coefs;
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_coefs;
   solv::LocationIndexRobinBcCoefs* d_quat_bc_coefs;

   std::shared_ptr<tbox::Database> d_integrator_db;

   // Timers
   std::shared_ptr<tbox::Timer> t_advance_timer;
   std::shared_ptr<tbox::Timer> t_eta_rhs_timer;
   std::shared_ptr<tbox::Timer> t_conc_rhs_timer;
   std::shared_ptr<tbox::Timer> t_symm_rhs_timer;
   std::shared_ptr<tbox::Timer> t_set_diffcoeff_conc_timer;
   std::shared_ptr<tbox::Timer> t_psolve_setup_timer;
   std::shared_ptr<tbox::Timer> t_psolve_solve_timer;
   std::shared_ptr<tbox::Timer> t_phase_precond_timer;
   std::shared_ptr<tbox::Timer> t_conc_precond_timer;
   std::shared_ptr<tbox::Timer> t_quat_grad_timer;
   std::shared_ptr<tbox::Timer> t_phase_rhs_timer;
};

#endif
