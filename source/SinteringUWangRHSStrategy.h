// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_SinteringUWangRHSStrategy
#define included_SinteringUWangRHSStrategy

#include "PhaseRHSStrategy.h"
#include "PhaseFluxStrategy.h"
#include "FreeEnergyStrategy.h"
#include "QuatModelParameters.h"
#include "QuatIntegrator.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include <string>

class SinteringUWangRHSStrategy : public PhaseRHSStrategy
{
 public:
   SinteringUWangRHSStrategy(
       const QuatModelParameters& _model_parameters, const int phase_id,
       const int conc_id, const int temperature_id, const int phase_mobility_id,
       const int flux_id, QuatIntegrator* time_integrator,
       std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
       std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy);

   void setup(std::shared_ptr<hier::PatchHierarchy> hierarchy);

   virtual void evaluateRHS(const double time,
                            std::shared_ptr<hier::PatchHierarchy> hierarchy,
                            const int ydot_phase_id, const bool eval_flag);

 private:
   const QuatModelParameters& d_model_parameters;

   /*
    *  model parameters
    */
   const double d_B;

   const int d_phase_id;
   const int d_conc_id;
   const int d_temperature_id;

   const int d_phase_mobility_id;
   const int d_flux_id;

   std::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

   QuatIntegrator* d_time_integrator;

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
