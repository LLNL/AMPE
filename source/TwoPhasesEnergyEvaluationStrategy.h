// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_TwoPhasesEnergyEvaluationStrategy
#define included_TwoPhasesEnergyEvaluationStrategy

#include "EnergyEvaluationStrategy.h"
#include "QuatModelParameters.h"

class TwoPhasesEnergyEvaluationStrategy : public EnergyEvaluationStrategy
{
 public:
   TwoPhasesEnergyEvaluationStrategy(
       const QuatModelParameters& model_parameters, const int qlen,
       const int phase_scratch_id, const int quat_scratch_id,
       const int quat_grad_side_id, const int weight_id, const int f_l_id,
       const int f_a_id, const int temperature_id, const int energy_diag_id);

   void evaluateEnergy(std::shared_ptr<hier::Patch> patch, const double time,
                       double& total_energy, double& total_phase_e,
                       double& total_orient_e, double& total_qint_e,
                       double& total_well_e, double& total_free_e,
                       const bool gp);

 private:
   const QuatModelParameters& d_model_parameters;
   const int d_qlen;

   const int d_phase_scratch_id;
   const int d_quat_scratch_id;
   const int d_quat_grad_side_id;
   const int d_weight_id;
   const int d_f_l_id;
   const int d_f_a_id;
   const int d_temperature_id;

   const int d_energy_diag_id;
};

#endif
