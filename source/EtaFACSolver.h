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
#ifndef included_EtaFACSolver
#define included_EtaFACSolver

#include "EllipticFACSolver.h"
#include "InterpolationType.h"
#include "EtaFACOps.h"

using namespace SAMRAI;

class EtaFACSolver : public EllipticFACSolver
{

 public:
   EtaFACSolver(const std::string &object_name,
                std::shared_ptr<EtaFACOps> fac_ops,
                const std::shared_ptr<tbox::Database> database =
                    std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(
       const int phase_id, const int eta_id, const int eta_mobility_id,
       const double epsilon_eta, const double gamma,
       const Thermo4PFM::EnergyInterpolationType phase_interp_func_type,
       const double eta_well_scale, const std::string eta_well_func_type);
};

#endif  // included_EtaFACSolver
