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
#include "TemperatureFACSolver.h"
#include "TemperatureFACOps.h"


TemperatureFACSolver::TemperatureFACSolver(
    const std::string &object_name, std::shared_ptr<TemperatureFACOps> fac_ops,
    const std::shared_ptr<tbox::Database> database)
    : EllipticFACSolver(object_name, fac_ops, database)
{
   t_set_op_coef = tbox::TimerManager::getManager()->getTimer(
       "AMPE::TemperatureFACSolver::setOperatorCoefficients");
}

// Set coefficiients for Eq. M div (D grad u) + C u = f
void TemperatureFACSolver::setOperatorCoefficients(const double m,
                                                   const double c,
                                                   const double d)
{
   assert(d < 0.);

   t_set_op_coef->start();

   std::shared_ptr<TemperatureFACOps> Temperature_fac_ops(
       std::dynamic_pointer_cast<TemperatureFACOps, EllipticFACOps>(d_fac_ops));

   Temperature_fac_ops->setOperatorCoefficients(m, c, d);

   finalizeCoefficients();

   t_set_op_coef->stop();
}
