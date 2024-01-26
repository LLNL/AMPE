// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "PhaseFACSolver.h"
#include "PhaseFACOps.h"


PhaseFACSolver::PhaseFACSolver(const std::string& object_name,
                               std::shared_ptr<PhaseFACOps>& fac_ops,
                               std::shared_ptr<tbox::Database> database)
    : EllipticFACSolver(object_name, fac_ops, database)
{
   t_set_op_coef = tbox::TimerManager::getManager()->getTimer(
       "AMPE::PhaseFACSolver::setOperatorCoefficients");
}

void PhaseFACSolver::setOperatorCoefficients(
    const int phase_id, const int phase_mobility_id, const double epsilon_phase,
    const double gamma, const double phase_well_scale,
    const std::string phase_well_func_type)
{
   assert(phase_id >= 0);
   assert(phase_mobility_id >= 0);

   t_set_op_coef->start();

   std::shared_ptr<PhaseFACOps> phase_fac_ops(
       std::dynamic_pointer_cast<PhaseFACOps, EllipticFACOps>(d_fac_ops));

   phase_fac_ops->setOperatorCoefficients(phase_id, phase_mobility_id,
                                          epsilon_phase, gamma,
                                          phase_well_scale,
                                          phase_well_func_type);

   finalizeCoefficients();

   t_set_op_coef->stop();
}
