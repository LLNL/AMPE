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
#ifndef included_TemperatureFACSolver
#define included_TemperatureFACSolver

#include "EllipticFACSolver.h"

class TemperatureFACOps;

using namespace SAMRAI;

class TemperatureFACSolver : public EllipticFACSolver
{

 public:
   TemperatureFACSolver(const std::string &object_name,
                        std::shared_ptr<TemperatureFACOps> fac_ops,
                        const std::shared_ptr<tbox::Database> database =
                            std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(const double, const double, const double);

 private:
   std::shared_ptr<tbox::Timer> t_set_op_coef;
};

#endif  // included_TemperatureFACSolver
