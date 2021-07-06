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
#ifndef included_PhaseRHSStrategyWithQ
#define included_PhaseRHSStrategyWithQ

#include "PhaseRHSStrategy.h"
#include "PhaseFluxStrategy.h"
#include "QuatIntegrator.h"
#include "FreeEnergyStrategy.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include <string>

class PhaseRHSStrategyWithQ : public PhaseRHSStrategy
{
 public:
   PhaseRHSStrategyWithQ(
       const QuatModelParameters& _model_parameters, const int phase_scratch_id,
       const int conc_scratch_id, const int quat_scratch_id,
       const int temperature_scratch_id, const int eta_scratch_id,
       const int f_l_id, const int f_a_id, const int f_b_id,
       const int phase_mobility_id, const int flux_id,
       const int quat_grad_modulus_id, const int noise_id,
       const int phase_rhs_visit_id, const int driving_force_visit_id,
       QuatIntegrator* integrator,
#ifdef USE_CPODE
       CPODESSolver* sundials_solver,
#else
       CVODESolver* sundials_solver,
#endif
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
       std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
       std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy);

   void setup(std::shared_ptr<hier::PatchHierarchy> hierarchy);

   virtual void evaluateRHS(const double time,
                            std::shared_ptr<hier::PatchHierarchy> hierarchy,
                            const int ydot_phase_id, const bool eval_flag,
                            const double frame_velocity);

 private:
   const QuatModelParameters& d_model_parameters;

   const double d_phase_well_scale;
   const double d_eta_well_scale;
   const double d_H_parameter;

   const EnergyInterpolationType d_energy_interp_func_type;
   std::string d_orient_interp_func_type;

   const int d_phase_scratch_id;
   const int d_conc_scratch_id;
   const int d_quat_scratch_id;
   const int d_temperature_scratch_id;
   const int d_eta_scratch_id;

   const int d_f_l_id;
   const int d_f_a_id;
   const int d_f_b_id;

   const int d_phase_mobility_id;
   const int d_flux_id;
   const int d_quat_grad_modulus_id;

   const int d_noise_id;

   const int d_phase_rhs_visit_id;
   const int d_driving_force_visit_id;

   QuatIntegrator* d_integrator;
#ifdef USE_CPODE
   CPODESSolver* d_sundials_solver;
#else
   CVODESolver* d_sundials_solver;
#endif
   std::shared_ptr<FreeEnergyStrategy> d_free_energy_strategy;

   double d_deltat;
   bool d_newtime;

   // Timers
   std::shared_ptr<tbox::Timer> t_phase_rhs_timer;

   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   std::shared_ptr<hier::CoarsenOperator> d_flux_coarsen_op;
   xfer::CoarsenAlgorithm d_flux_coarsen_algorithm;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> > d_flux_coarsen_schedule;

   std::shared_ptr<PhaseFluxStrategy> d_phase_flux_strategy;

   // internal functions
   void evaluateRHS(const double time, std::shared_ptr<hier::Patch> patch,
                    const int ydot_phase_id, const bool eval_flag);
};

#endif
