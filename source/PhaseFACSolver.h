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
#ifndef included_PhaseFACSolver
#define included_PhaseFACSolver

#include "EllipticFACSolver.h"
#include "InterpolationType.h"
#include "PhaseFACOps.h"

using namespace SAMRAI;

class PhaseFACSolver : public EllipticFACSolver
{

 public:
   PhaseFACSolver(const std::string& object_name,
                  std::shared_ptr<PhaseFACOps>& fac_ops,
                  const std::shared_ptr<tbox::Database> database =
                      std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(const int phase_id, const int phase_mobility_id,
                                const double epsilon_phase, const double gamma,
                                const double phase_well_scale,
                                const std::string phase_well_func_type);

 private:
   std::shared_ptr<tbox::Timer> t_set_op_coef;
};

#endif  // included_PhaseFACSolver
