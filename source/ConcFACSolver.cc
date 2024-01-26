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
//
#include "ConcFACSolver.h"


ConcFACSolver::ConcFACSolver(const std::string& object_name,
                             std::shared_ptr<ConcFACOps> fac_ops,
                             const std::shared_ptr<tbox::Database>& database)
    : EllipticFACSolver(object_name, fac_ops, database)
{
   t_set_op_coef = tbox::TimerManager::getManager()->getTimer(
       "AMPE::ConcFACSolver::setOperatorCoefficients");
}

void ConcFACSolver::setOperatorCoefficients(
    const double gamma, const std::vector<int>& diffusion_id,
    const double mobility)
{
   t_set_op_coef->start();

   std::shared_ptr<ConcFACOps> conc_fac_ops(
       std::dynamic_pointer_cast<ConcFACOps, EllipticFACOps>(d_fac_ops));

   conc_fac_ops->setOperatorCoefficients(gamma, diffusion_id, mobility);

   finalizeCoefficients();

   t_set_op_coef->stop();
}
