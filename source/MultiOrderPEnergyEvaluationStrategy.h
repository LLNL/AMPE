// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_MultiOrderPEnergyEvaluationStrategy
#define included_MultiOrderPEnergyEvaluationStrategy

#include "QuatModelParameters.h"

class MultiOrderPEnergyEvaluationStrategy
{
 public:
   MultiOrderPEnergyEvaluationStrategy(
       const QuatModelParameters& model_parameters, const int phase_scratch_id,
       const int weight_id, const int energy_diag_id);

   void evaluatePairEnergy(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void printPairEnergy(std::ostream& os);

 private:
   const QuatModelParameters& d_model_parameters;

   const int d_phase_id;
   const int d_weight_id;

   const int d_energy_diag_id;

   std::vector<double> d_total_energy;

   void evaluatePairEnergy(std::shared_ptr<hier::Patch> patch,
                           std::vector<double>& total_energy);
};

#endif
