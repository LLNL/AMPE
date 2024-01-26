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
#ifndef included_QuatMobilityStrategy
#define included_QuatMobilityStrategy

#include "SAMRAI/hier/PatchHierarchy.h"
using namespace SAMRAI;

class QuatMobilityStrategy
{
 public:
   enum CACHE_TYPE { CACHE = 0, FORCE = 1 };

   QuatMobilityStrategy(){};

   virtual ~QuatMobilityStrategy(){};

   virtual void computePhaseMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE) = 0;

   virtual void computeEtaMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& eta_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE) = 0;

   virtual void computeQuatMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE) = 0;

   virtual void computeQuatMobilityDeriv(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_deriv_id, const double time,
       const CACHE_TYPE cache = CACHE) = 0;

   virtual void computePhaseTemperatureMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int&, int&, int&);
};

#endif
