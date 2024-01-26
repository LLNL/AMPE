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
#ifndef included_HeatCapacityStrategy
#define included_HeatCapacityStrategy

#include "SAMRAI/hier/PatchHierarchy.h"
using namespace SAMRAI;

class HeatCapacityStrategy
{
 public:
   HeatCapacityStrategy(){};

   virtual ~HeatCapacityStrategy(){};

   virtual void setCurrentValue(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy) = 0;
};

#endif
