// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_EnergyEvaluationStrategy
#define included_EnergyEvaluationStrategy

#include "QuatModelParameters.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Patch.h"

#include <memory>

using namespace SAMRAI;

class EnergyEvaluationStrategy
{
 public:
   EnergyEvaluationStrategy(const QuatModelParameters& model_parameters,
                            const int energy_diag_id)
       : d_model_parameters(model_parameters),
         d_energy_diag_id(energy_diag_id){};

   void evaluateEnergy(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                       const double time, double& total_energy,
                       double& total_phase_e, double& total_orient_e,
                       double& total_qint_e, double& total_well_e,
                       double& total_free_e, const bool gp);

   virtual void evaluateEnergy(std::shared_ptr<hier::Patch> patch,
                               const double time, double& total_energy,
                               double& total_phase_e, double& total_orient_e,
                               double& total_qint_e, double& total_well_e,
                               double& total_free_e, const bool gp) = 0;

 private:
   const QuatModelParameters& d_model_parameters;

   const int d_energy_diag_id;
};

#endif
