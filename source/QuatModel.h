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
#ifndef included_QuatModel
#define included_QuatModel

#include "PhaseFreeEnergyStrategy.h"
#include "QuatModelParameters.h"

// Headers for basic SAMRAI objects
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/IEEE.h"

// Headers for major algorithm/data structure objects
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"

#include "FieldsInitializer.h"

#include "PFModel.h"
#include "QuatIntegrator.h"
#include "QuatRefinePatchStrategy.h"
#include "PartitionCoeffRefinePatchStrategy.h"
#include "CALPHADMobility.h"
#include "QuatParams.h"
#include "Grains.h"
#include "CALPHADFreeEnergyFunctionsWithPenaltyBinary.h"

#include "EventInterval.h"

#include <map>

using namespace SAMRAI;

class GradStrategy;
class QuatGradStrategy;
class FreeEnergyStrategy;
class CompositionRHSStrategy;
class QuatVerbosity;
class MolarVolumeStrategy;
class HeatCapacityStrategy;
class PhaseFluxStrategy;
class PhaseConcentrationsStrategy;
class MeltingTemperatureStrategy;
class CompositionStrategyMobilities;
class CompositionDiffusionStrategy;

class QuatModel : public PFModel
{
 public:
   enum CACHE_TYPE { CACHE = 0, FORCE = 1 };

   QuatModel(int qlen);
   virtual ~QuatModel();

   virtual void Initialize(std::shared_ptr<tbox::MemoryDatabase>& input_db,
                           const std::string& run_name,
                           const bool is_from_restart,
                           const std::string& restart_read_dirname,
                           const int restore_num);

   virtual void Run(void);

   virtual void CreateIntegrator(std::shared_ptr<tbox::Database> input_db);

   virtual void RegisterVariables(void);

   virtual void initializeCoarseRefineOperators(void);

   virtual void InitializeIntegrator(void);

   virtual void RegisterWithVisit(void);

   virtual double Advance(void);

   virtual void postAdvanceDiagnostics(void);

   virtual void preRunDiagnostics(void);

   virtual void postRunDiagnostics(void);

   virtual void writeRestartFile(void);

   virtual void Regrid(const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void printScalarDiagnostics(void);

   void listLocalToGlobal(std::map<int, double>& local_map,
                          std::map<int, double>& global_map);

   void readInitialDatabase(std::shared_ptr<tbox::Database> main_input_db);

   void WriteInitialConditionsFile(void);

   std::shared_ptr<hier::PatchLevel> FlattenHierarchy(
       const std::shared_ptr<hier::PatchHierarchy> src_hierarchy,
       const int level_number, const double time);

   void setupInitialDataLevel(void);

   void setupHierarchy(void);

   void computeMinMaxQModulus(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void tagGradientDetectorCells(hier::Patch& patch, const double regrid_time,
                                 const bool initial_error, const int tag_index,
                                 const bool uses_richardson_extrapolation_too);

   void computePhaseDiffs(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& phase_diffs_id, const double time, const CACHE_TYPE cache = CACHE);
   void computePhaseDiffs(const std::shared_ptr<hier::PatchLevel> patch_level,
                          int& phase_id, int& phase_diffs_id,
                          const double time);

   void computeEtaDiffs(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int& eta_id, int& eta_diffs_id, const double time,
                        const CACHE_TYPE cache = CACHE);
   void computeEtaDiffs(const std::shared_ptr<hier::PatchLevel> patch_level,
                        int& eta_id, int& eta_diffs_id, const double time);

   void computeVarDiffs(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int& var_id, int& diffs_id, const double time,
                        const CACHE_TYPE cache = CACHE);
   void computeVarDiffs(const std::shared_ptr<hier::PatchLevel> patch_level,
                        int& var_id, int& diffs_id, const double time);

   void computePhaseGradCell(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       int& phase_diffs_id, int& phase_grad_cell_id, const double time,
       const CACHE_TYPE cache = CACHE);
   void computePhaseGradCell(
       const std::shared_ptr<hier::PatchLevel> patch_level,
       int& phase_diffs_id, int& phase_grad_cell_id, const double time);

   void computeEtaGradCell(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       int& eta_diffs_id, int& eta_grad_cell_id, const double time,
       const CACHE_TYPE cache = CACHE);
   void computeEtaGradCell(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& eta_diffs_id,
       int& eta_grad_cell_id, const double time);

   void computeVarGradCell(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_cell_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeVarGradCell(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_cell_id, const double time);

   void computeVarGradSide(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_side_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeVarGradSide(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_side_id, const double time);

   void computeSymmetryRotations(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

   void makeQuatFundamental(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

   void computeQuatDiffs(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& quat_id,
       int& diffs_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeQuatDiffs(const std::shared_ptr<hier::PatchLevel> patch_level,
                         int& quat_id, int& diffs_id, const double time);

   void computeQuatGradCell(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_cell_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeQuatGradCell(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_cell_id, const double time);

   void computeQuatGradSide(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_side_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeQuatGradSide(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_side_id, const double time);

   void computeQuatGradModulus(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       int& grad_cell_id, int& grad_modulus_id, const double time,
       const CACHE_TYPE cache = CACHE);
   void computeQuatGradModulus(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& grad_cell_id,
       int& grad_modulus_id, const double time);

   void computeQuatGradModulusFromSides(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       int& grad_side_id, int& grad_modulus_id, const double time,
       const CACHE_TYPE cache = CACHE);
   void computeQuatGradModulusFromSides(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& grad_cell_id,
       int& grad_modulus_id, const double time);

   void normalizeQuat(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                      const int quat_id);
   void normalizeQuat(const std::shared_ptr<hier::PatchLevel> patch_level,
                      const int quat_id);

   void computeUniformPhaseMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeUniformPhaseMobility(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& phase_id,
       int& mobility_id, const double time);

   void computeEtaMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeEtaMobility(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& phase_id,
       int& mobility_id, const double time);

   void computeQuatMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);
   void computeQuatMobility(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& phase_id,
       int& mobility_id, const double time);

   void computeQuatMobilityDeriv(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_deriv_id, const double time,
       const CACHE_TYPE cache = CACHE);
   void computeQuatMobilityDeriv(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& phase_id,
       int& mobility_deriv_id, const double time);

   void computeVectorWeights(std::shared_ptr<hier::PatchHierarchy>, int, int);

   void evaluateEnergy(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                       const double time, double& total_energy,
                       double& total_phase_e, double& total_eta_e,
                       double& total_orient_e, double& total_qint_e,
                       double& total_well_e, double& total_free_e,
                       const bool gp = false);
   double evaluateIntegralConcentration(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int component);
   double evaluateVolumeSolid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);
   double evaluateVolumeEta(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);
   double evaluateIntegralPhaseConcentration(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int component);

   void computeVelocity(std::shared_ptr<hier::Patch> patch, int phi_dot_id);
   void computeVelocity(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int phi_dot_id);
   int getConcDiffusionLid() const { return d_conc_pfm_diffusion_l_id; }
   int getConcDiffusionAid() const { return d_conc_pfm_diffusion_a_id; }

   //-----------------------------------------------------------------------
   //
   // Methods inherited from Serializable (through PFModel)
   //
   virtual void putToRestart(const std::shared_ptr<tbox::Database>& db) const;

   //-----------------------------------------------------------------------
   //
   // Methods inherited from StandardTagAndInitStrategy (through PFModel)
   //
   virtual void initializeLevelData(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int level_number, const double time, const bool can_be_refined,
       const bool initial_time,
       const std::shared_ptr<hier::PatchLevel>& old_level,
       const bool allocate_data);

   virtual void resetHierarchyConfiguration(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int coarsest_level, const int finest_level);

   virtual void applyGradientDetector(
       std::shared_ptr<hier::PatchHierarchy>& hierarchy, int level_number,
       double time, int tag_index, const bool initial_time,
       const bool uses_richardson_extrapolation_too);

   //-----------------------------------------------------------------------

   bool resetGrains(void);

   void fillPhaseConcentrationGhosts(void);
   void fillPartitionCoeffGhosts(void);

   bool isSymmetryAware(void);

   // deallocate some temporary data to free some memory,
   // for example before some high memory footprint postprocessing
   void DeallocateIntermediateLocalPatchData(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void registerPhaseConcentrationVariables(
       std::shared_ptr<pdat::CellVariable<double> > conc_l_var,
       std::shared_ptr<pdat::CellVariable<double> > conc_a_var,
       std::shared_ptr<pdat::CellVariable<double> > conc_b_var);

   double computeThermalEnergy(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

 private:
   void copyCurrentToScratch(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time, QuatRefinePatchStrategy* patch_strategy);

   void copyCurrentToScratch(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int ln,
       const double time, QuatRefinePatchStrategy* patch_strategy);

   void AllocateLocalPatchData(const std::shared_ptr<hier::PatchLevel> level,
                               const double time, const bool initial_time);

   void DeallocateIntermediateLocalPatchData(
       const std::shared_ptr<hier::PatchLevel> level);

   void AllocateQuatLocalPatchData(
       const std::shared_ptr<hier::PatchLevel> level, const double time,
       const bool initial_time);

   template <typename T>
   void AllocateAndZeroData(const int data_id,
                            const std::shared_ptr<hier::PatchLevel> level,
                            const double time, const bool initial_time);

   void checkQuatNorm(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                      const double tol);

   void computePhaseMobilityPatch(
       const hier::Box& pbox,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_mobility);

   void computeEtaMobilityPatch(
       const hier::Box& pbox,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_mobility,
       std::shared_ptr<pdat::CellData<double> > cd_phi);

   void preRunDiagnosticsMobilityInPhases(const double temperature);
   bool computeCeq(const double temperature, const PhaseIndex pi0,
                   const PhaseIndex pi1, double* ceq) const;

   void applyPolynomial(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        const int src_cell_data_id, const int dst_cell_data_id);
   void applyPolynomial(const std::shared_ptr<hier::PatchLevel> patch_level,
                        const int src_cell_data_id, const int dst_cell_data_id);
   void readAMRdatabase(std::shared_ptr<tbox::Database> amr_db);
   void registerPatchDataForRestart(hier::VariableDatabase*);
   void registerConcentrationVariables(void);
   void registerPhaseConcentrationVariables();
   void registerPhaseVariables(void);
   void registerEtaVariables(void);
   void registerOrientationVariables(void);
   void registerGrainVariables(void);
   void registerPatchDataForRestart(void);
   void initializeRefineCoarsenAlgorithms();
   void initializeTemperature(std::shared_ptr<tbox::Database>,
                              std::shared_ptr<tbox::Database>);
   void initializeAmr(std::shared_ptr<tbox::Database> amr_db);

   void resetRefPhaseConcentrations();
   void setPhaseConcentrationsToEquilibrium(const double* const ceq);
   void setRefPhaseConcentrationsToEquilibrium(const double* const ceq);

   void findAndNumberGrains(void);
   void computeGrainDiagnostics(void);
   void extendGrainOrientation(void);

   void smoothQuat(const std::shared_ptr<hier::PatchLevel> patch_level);
   void smoothQuat(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                   const double time);
   void initializeRHSandEnergyStrategies(
       std::shared_ptr<tbox::MemoryDatabase>& quat_db);
   void initializeCompositionRHSStrategy();

   QuatModelParameters d_model_parameters;

   int d_qlen;
   int d_ncompositions;

   QuatVerbosity* d_verbosity;

   bool d_symmetry_aware;
   bool d_extra_energy_detail;

   std::shared_ptr<EventInterval> d_test_interval;
   std::shared_ptr<EventInterval> d_fundamental_interval;
   std::shared_ptr<EventInterval> d_scalar_diag_interval;
   std::shared_ptr<EventInterval> d_grain_extend_interval;

   std::shared_ptr<pdat::CellVariable<double> > d_phase_var;
   int d_phase_id;
   int d_phase_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_eta_var;
   int d_eta_id;
   int d_eta_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_quat_var;
   int d_quat_id;
   int d_quat_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_quat_relax_var;
   int d_quat_relax_id;
   int d_quat_relax_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_conc_var;
   int d_conc_id;
   int d_conc_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_conc_l_var;
   int d_conc_l_id;
   int d_conc_l_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_conc_a_var;
   int d_conc_a_id;
   int d_conc_a_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_conc_b_var;
   int d_conc_b_id;
   int d_conc_b_scratch_id;
   std::shared_ptr<pdat::CellVariable<double> > d_conc_l_ref_var;
   int d_conc_l_ref_id;
   std::shared_ptr<pdat::CellVariable<double> > d_conc_a_ref_var;
   int d_conc_a_ref_id;
   std::shared_ptr<pdat::CellVariable<double> > d_conc_b_ref_var;
   int d_conc_b_ref_id;

   std::shared_ptr<pdat::SideVariable<double> > d_phase_diffs_var;
   int d_phase_diffs_id;
   std::shared_ptr<pdat::SideVariable<double> > d_eta_diffs_var;
   int d_eta_diffs_id;
   std::shared_ptr<pdat::SideVariable<double> > d_quat_diffs_var;
   int d_quat_diffs_id;

   std::shared_ptr<pdat::CellVariable<double> > d_phase_diffs_cell_var;
   int d_phase_diffs_cell_id;

   std::shared_ptr<pdat::SideVariable<int> > d_quat_symm_rotation_var;
   int d_quat_symm_rotation_id;

   std::shared_ptr<pdat::CellVariable<int> > d_quat_symm_rotation_cell_var;
   int d_quat_symm_rotation_cell_id;

   std::shared_ptr<pdat::CellVariable<double> > d_quat_diffs_cell_var;
   int d_quat_diffs_cell_id;
   std::shared_ptr<pdat::CellVariable<double> > d_quat_nonsymm_diffs_cell_var;
   int d_quat_nonsymm_diffs_cell_id;

   std::shared_ptr<pdat::CellVariable<double> > d_quat_norm_error_var;
   int d_quat_norm_error_id;

   std::shared_ptr<pdat::CellVariable<double> > d_phase_grad_cell_var;
   int d_phase_grad_cell_id;
   std::shared_ptr<pdat::SideVariable<double> > d_phase_grad_side_var;
   int d_phase_grad_side_id;

   std::shared_ptr<pdat::CellVariable<double> > d_eta_grad_cell_var;
   int d_eta_grad_cell_id;
   std::shared_ptr<pdat::SideVariable<double> > d_eta_grad_side_var;
   int d_eta_grad_side_id;

   std::shared_ptr<pdat::CellVariable<double> > d_quat_grad_cell_var;
   int d_quat_grad_cell_id;
   std::shared_ptr<pdat::SideVariable<double> > d_quat_grad_side_var;
   int d_quat_grad_side_id;
   std::shared_ptr<pdat::CellVariable<double> > d_quat_grad_modulus_var;
   int d_quat_grad_modulus_id;

   std::shared_ptr<pdat::CellVariable<double> > d_phase_mobility_var;
   int d_phase_mobility_id;
   std::shared_ptr<pdat::CellVariable<double> > d_eta_mobility_var;
   int d_eta_mobility_id;
   std::shared_ptr<pdat::CellVariable<double> > d_quat_mobility_var;
   int d_quat_mobility_id;

   std::shared_ptr<pdat::SideVariable<double> > d_quat_diffusion_var;
   int d_quat_diffusion_id;

   std::shared_ptr<pdat::SideVariable<double> > d_conc_diffusion_var;
   int d_conc_diffusion_id;

   std::vector<std::shared_ptr<pdat::SideVariable<double> > >
       d_conc_pfm_diffusion_var;
   std::vector<int> d_conc_pfm_diffusion_id;

   /*!
    * holds data for diffusion coefficients in composition equation
    * according to EBS scheme, including weight due to phase fraction
    */
   std::shared_ptr<pdat::SideVariable<double> > d_conc_pfm_diffusion_l_var;
   int d_conc_pfm_diffusion_l_id;
   std::shared_ptr<pdat::SideVariable<double> > d_conc_pfm_diffusion_a_var;
   int d_conc_pfm_diffusion_a_id;
   std::shared_ptr<pdat::SideVariable<double> > d_conc_pfm_diffusion_b_var;
   int d_conc_pfm_diffusion_b_id;


   /*!
    * holds data for diffusion coefficients in each individual phase
    */
   std::shared_ptr<pdat::SideVariable<double> > d_conc_diffusion_coeff_l_var;
   int d_conc_diffusion_coeff_l_id;
   std::shared_ptr<pdat::SideVariable<double> > d_conc_diffusion_coeff_a_var;
   int d_conc_diffusion_coeff_a_id;
   std::shared_ptr<pdat::SideVariable<double> > d_conc_diffusion_coeff_b_var;
   int d_conc_diffusion_coeff_b_id;

   std::shared_ptr<pdat::SideVariable<double> > d_conc_Mq_var;
   int d_conc_Mq_id;

   std::shared_ptr<pdat::SideVariable<double> >
       d_conc_phase_coupling_diffusion_var;
   int d_conc_phase_coupling_diffusion_id;
   std::shared_ptr<pdat::SideVariable<double> >
       d_conc_eta_coupling_diffusion_var;
   int d_conc_eta_coupling_diffusion_id;

   std::shared_ptr<pdat::CellVariable<double> > d_f_l_var;
   int d_f_l_id;
   std::shared_ptr<pdat::CellVariable<double> > d_f_a_var;
   int d_f_a_id;
   std::shared_ptr<pdat::CellVariable<double> > d_f_b_var;
   int d_f_b_id;

   std::shared_ptr<pdat::CellVariable<double> > d_cp_var;
   int d_cp_id;

   std::shared_ptr<pdat::CellVariable<double> > d_velocity_var;
   int d_velocity_id;

   std::shared_ptr<pdat::CellVariable<double> > d_partition_coeff_var;
   int d_partition_coeff_id;
   int d_partition_coeff_scratch_id;

   std::shared_ptr<pdat::CellVariable<double> > d_equilibrium_temperature_var;
   int d_equilibrium_temperature_id;

   std::shared_ptr<pdat::CellVariable<double> > d_energy_diag_var;
   int d_energy_diag_id;


   /*
    * Variable containing volume weights used in composite grid norm
    * calculations.
    */
   std::shared_ptr<pdat::CellVariable<double> > d_weight_var;
   int d_weight_id;

   /*
    * Variable containing temporary work data of depth 1
    */
   std::shared_ptr<pdat::CellVariable<double> > d_work_var;
   int d_work_id;

   /*!
    * Temperature field.
    * Typically needs ghosts values, even if not solution of PDE.
    * Ghosts values are also use to get side values through averaging
    */
   std::shared_ptr<pdat::CellVariable<double> > d_temperature_var;
   int d_temperature_id;
   int d_temperature_scratch_id;

   /*!
    * RHS needed for temperature when result of steady state equation
    */
   std::shared_ptr<pdat::CellVariable<double> > d_temperature_rhs_steady_var;
   int d_temperature_rhs_steady_id;

   bool d_is_from_restart;

   std::shared_ptr<QuatIntegrator> d_integrator;
   std::shared_ptr<QuatIntegrator> d_integrator_quat_only;

   QuatRefinePatchStrategy* d_all_refine_patch_strategy;
   PartitionCoeffRefinePatchStrategy* d_partition_coeff_refine_patch_strategy;

   std::shared_ptr<hier::RefineOperator> d_phase_refine_op;
   std::shared_ptr<hier::RefineOperator> d_eta_refine_op;
   std::shared_ptr<hier::RefineOperator> d_quat_refine_op;
   std::shared_ptr<hier::RefineOperator> d_conc_refine_op;
   std::shared_ptr<hier::CoarsenOperator> d_quat_coarsen_op;

   std::shared_ptr<xfer::RefineAlgorithm> d_regrid_refine_alg;
   std::shared_ptr<xfer::RefineAlgorithm> d_curr_to_curr_refine_alg;
   std::shared_ptr<xfer::RefineAlgorithm> d_curr_to_scr_refine_alg;

   std::vector<std::shared_ptr<xfer::RefineSchedule> >
       d_curr_to_scr_refine_sched;

   QuatGradStrategy* d_quat_grad_strategy;
   std::shared_ptr<QuatMobilityStrategy> d_mobility_strategy;
   FreeEnergyStrategy* d_free_energy_strategy;
   FreeEnergyStrategy* d_free_energy_strategy_for_diffusion;
   PhaseConcentrationsStrategy* d_phase_conc_strategy;
   PartitionCoefficientStrategy* d_partition_coeff_strategy;
   TemperatureStrategy* d_temperature_strategy;
   TemperatureStrategy* d_temperature_strategy_quat_only;
   CompositionRHSStrategy* d_composition_rhs_strategy;

   HeatCapacityStrategy* d_heat_capacity_strategy;
   MolarVolumeStrategy* d_mvstrategy;
   FreeEnergyFunctions* d_cafe;
   PhaseFluxStrategy* d_phase_flux_strategy;
   MeltingTemperatureStrategy* d_meltingT_strategy;
   CompositionStrategyMobilities* d_composition_strategy_mobilities;
   std::shared_ptr<CompositionDiffusionStrategy>
       d_diffusion_for_conc_in_phase;

   bool d_tag_phase;
   bool d_tag_eta;
   bool d_tag_quat;
   double d_phase_threshold_tagged;
   double d_phase_threshold_untagged;
   double d_eta_threshold_tagged;
   double d_eta_threshold_untagged;
   double d_quat_threshold_tagged;
   double d_quat_threshold_untagged;

   std::string d_fenergy_diag_filename;

   bool d_use_warm_start;

   std::shared_ptr<tbox::Database> d_conc_db;
   std::shared_ptr<tbox::Database> d_calphad_db;
   std::shared_ptr<tbox::Database> d_newton_db;

   std::shared_ptr<Grains> d_grains;
   int d_number_of_grains;
   double d_phase_threshold;

   // Timers
   std::shared_ptr<tbox::Timer> t_resetGrains_timer;
   std::shared_ptr<tbox::Timer> t_phase_diffs_timer;
};

#endif
