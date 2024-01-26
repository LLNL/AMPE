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
#ifndef included_ConcFACSolver
#define included_ConcFACSolver

#include "EllipticFACSolver.h"
#include "ConcFACOps.h"

using namespace SAMRAI;

class ConcFACSolver : public EllipticFACSolver
{

 public:
   ConcFACSolver(const std::string& object_name,
                 std::shared_ptr<ConcFACOps> fac_ops,
                 const std::shared_ptr<tbox::Database>& database =
                     std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(const double gamma,
                                const std::vector<int>& diffusion_id,
                                const double mobility);

 private:
   std::shared_ptr<tbox::Timer> t_set_op_coef;
};

#endif  // included_ConcFACSolver
