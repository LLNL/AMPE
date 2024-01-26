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
#include "SimpleQuatMobilityStrategy.h"

#include <cassert>

SimpleQuatMobilityStrategy::SimpleQuatMobilityStrategy(QuatModel* quat_model)
{
   d_quat_model = quat_model;
}

//-----------------------------------------------------------------------

SimpleQuatMobilityStrategy::~SimpleQuatMobilityStrategy()
{
   d_quat_model = nullptr;
}

//-----------------------------------------------------------------------

enum QuatModel::CACHE_TYPE SimpleQuatMobilityStrategy::translateCacheType(
    const CACHE_TYPE cache)
{
   enum QuatModel::CACHE_TYPE qm_cache;

   if (cache == FORCE) {
      qm_cache = QuatModel::FORCE;
   } else {
      qm_cache = QuatModel::CACHE;
   }

   return qm_cache;
}

//-----------------------------------------------------------------------

void SimpleQuatMobilityStrategy::computePhaseMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeUniformPhaseMobility(hierarchy, phase_id, mobility_id,
                                             time, qm_cache);
}

//-----------------------------------------------------------------------

void SimpleQuatMobilityStrategy::computeEtaMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& eta_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeEtaMobility(hierarchy, eta_id, mobility_id, time,
                                    qm_cache);
}

//-----------------------------------------------------------------------

void SimpleQuatMobilityStrategy::computeQuatMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeQuatMobility(hierarchy, phase_id, mobility_id, time,
                                     qm_cache);
}

//-----------------------------------------------------------------------

void SimpleQuatMobilityStrategy::computeQuatMobilityDeriv(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_deriv_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeQuatMobilityDeriv(hierarchy, phase_id,
                                          mobility_deriv_id, time, qm_cache);
}

//-----------------------------------------------------------------------
