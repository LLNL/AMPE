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
#ifndef included_SimpleQuatMobilityStrategy
#define included_SimpleQuatMobilityStrategy

#include "QuatMobilityStrategy.h"
#include "QuatModel.h"

class SimpleQuatMobilityStrategy : public QuatMobilityStrategy
{
 public:
   SimpleQuatMobilityStrategy(QuatModel* quat_model);

   ~SimpleQuatMobilityStrategy();

   virtual void computePhaseMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeEtaMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& eta_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeQuatMobility(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeQuatMobilityDeriv(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
       int& mobility_deriv_id, const double time,
       const CACHE_TYPE cache = CACHE);

 private:
   QuatModel* d_quat_model;

   enum QuatModel::CACHE_TYPE translateCacheType(const CACHE_TYPE cache);
};

#endif
