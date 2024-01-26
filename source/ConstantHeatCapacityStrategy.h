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
#ifndef included_ConstantHeatCapacityStrategy
#define included_ConstantHeatCapacityStrategy

#include "HeatCapacityStrategy.h"

// Simple HeatCapacityStrategy:
// set same value at every cell once for all
class ConstantHeatCapacityStrategy : public HeatCapacityStrategy
{
 public:
   ConstantHeatCapacityStrategy(const double cp, const int cp_id);
   void setCurrentValue(std::shared_ptr<hier::PatchHierarchy> patch_hierarchy);

 private:
   double d_cp;

   int d_cp_id;
};

#endif
