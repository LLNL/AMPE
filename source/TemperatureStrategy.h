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
#ifndef included_TemperatureStrategy
#define included_TemperatureStrategy

#include "SAMRAI/hier/PatchHierarchy.h"
using namespace SAMRAI;

class TemperatureStrategy
{
 public:
   TemperatureStrategy(){};

   virtual ~TemperatureStrategy(){};

   virtual void initialize(
       const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy)
   {
      (void)patch_hierarchy;
   };

   virtual double getCurrentMaxTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time) = 0;
   virtual double getCurrentMinTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time) = 0;
   virtual double getCurrentAverageTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time) = 0;

   // set temperature field according to specific strategy
   virtual void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time) = 0;

   virtual void resetSolversState(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int coarsest_level, const int finest_level)
   {
      (void)hierarchy;
      (void)coarsest_level;
      (void)finest_level;
   };
};

#endif
