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
#include "EtaFACSolver.h"
#include "EtaFACOps.h"

EtaFACSolver::EtaFACSolver(const std::string &object_name,
                           std::shared_ptr<EtaFACOps> fac_ops,
                           std::shared_ptr<tbox::Database> database)
    : EllipticFACSolver(object_name, fac_ops, database)
{
   return;
}

void EtaFACSolver::setOperatorCoefficients(
    const int phase_id, const int eta_id, const int eta_mobility_id,
    const double epsilon_eta, const double gamma,
    const Thermo4PFM::EnergyInterpolationType phase_interp_func_type,
    const double eta_well_scale, const std::string eta_well_func_type)
{
   std::shared_ptr<EtaFACOps> eta_fac_ops(
       std::dynamic_pointer_cast<EtaFACOps, EllipticFACOps>(d_fac_ops));

   eta_fac_ops->setOperatorCoefficients(phase_id, eta_id, eta_mobility_id,
                                        epsilon_eta, gamma,
                                        phase_interp_func_type, eta_well_scale,
                                        eta_well_func_type);

   finalizeCoefficients();
}
