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
#include "SimpleQuatGradStrategy.h"

SimpleQuatGradStrategy::SimpleQuatGradStrategy(QuatModel* quat_model)
{
   d_quat_model = quat_model;
}

//-----------------------------------------------------------------------

SimpleQuatGradStrategy::~SimpleQuatGradStrategy() { d_quat_model = nullptr; }

//-----------------------------------------------------------------------

bool SimpleQuatGradStrategy::isSymmetryAware(void)
{
   return d_quat_model->isSymmetryAware();
}

//-----------------------------------------------------------------------

enum QuatModel::CACHE_TYPE SimpleQuatGradStrategy::translateCacheType(
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

void SimpleQuatGradStrategy::computeDiffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& quat_id,
    int& diffs_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeQuatDiffs(hierarchy, quat_id, diffs_id, time, qm_cache);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeDiffs(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& quat_id,
    int& diffs_id, const double time)
{
   assert(d_quat_model != nullptr);

   d_quat_model->computeQuatDiffs(patch_level, quat_id, diffs_id, time);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeGradCell(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
    int& grad_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeQuatGradCell(hierarchy, diffs_id, grad_id, time,
                                     qm_cache);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeGradCell(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
    int& grad_id, const double time)
{
   assert(d_quat_model != nullptr);

   d_quat_model->computeQuatGradCell(patch_level, diffs_id, grad_id, time);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeGradSide(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
    int& grad_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeQuatGradSide(hierarchy, diffs_id, grad_id, time,
                                     qm_cache);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeGradSide(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
    int& grad_id, const double time)
{
   assert(d_quat_model != nullptr);

   d_quat_model->computeQuatGradSide(patch_level, diffs_id, grad_id, time);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeGradModulus(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_id,
    int& mod_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeQuatGradModulus(hierarchy, grad_id, mod_id, time,
                                        qm_cache);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeGradModulus(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& grad_id,
    int& mod_id, const double time)
{
   assert(d_quat_model != nullptr);

   d_quat_model->computeQuatGradModulus(patch_level, grad_id, mod_id, time);
}

//-----------------------------------------------------------------------

void SimpleQuatGradStrategy::computeGradModulusFromSides(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_side_id,
    int& mod_id, const double time, const CACHE_TYPE cache)
{
   assert(d_quat_model != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache = translateCacheType(cache);

   d_quat_model->computeQuatGradModulusFromSides(hierarchy, grad_side_id,
                                                 mod_id, time, qm_cache);
}
