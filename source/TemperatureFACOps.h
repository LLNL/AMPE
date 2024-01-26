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
#ifndef included_TemperatureFACOps
#define included_TemperatureFACOps

#include "EllipticFACOps.h"

#include <string>

using namespace SAMRAI;

class TemperatureFACOps : public EllipticFACOps
{

 public:
   TemperatureFACOps(const std::string &object_name = std::string(),
                     const std::shared_ptr<tbox::Database> database =
                         std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(const double, const double, const double);

 private:
};

#endif  // included_TemperatureFACOps
