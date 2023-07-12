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
#include "SimpleGradStrategy.h"

#include <cassert>

SimpleGradStrategy::SimpleGradStrategy(QuatModel* model) { d_pfmodel = model; }

//-----------------------------------------------------------------------

SimpleGradStrategy::~SimpleGradStrategy() { d_pfmodel = nullptr; }

//-----------------------------------------------------------------------

void SimpleGradStrategy::computeDiffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& var_id,
    int& diffs_id, const double time, const CACHE_TYPE cache)
{
   assert(d_pfmodel != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache;
   if (cache == FORCE) {
      qm_cache = QuatModel::FORCE;
   } else {
      qm_cache = QuatModel::CACHE;
   }

   d_pfmodel->computeVarDiffs(hierarchy, var_id, diffs_id, time, qm_cache);
}

//-----------------------------------------------------------------------

void SimpleGradStrategy::computeGradCell(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
    int& grad_id, const double time, const CACHE_TYPE cache)
{
   assert(d_pfmodel != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache;
   if (cache == FORCE) {
      qm_cache = QuatModel::FORCE;
   } else {
      qm_cache = QuatModel::CACHE;
   }

   d_pfmodel->computeVarGradCell(hierarchy, diffs_id, grad_id, time, qm_cache);
}

//-----------------------------------------------------------------------

void SimpleGradStrategy::computeGradSide(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
    int& grad_id, const double time, const CACHE_TYPE cache)
{
   assert(d_pfmodel != nullptr);

   enum QuatModel::CACHE_TYPE qm_cache;
   if (cache == FORCE) {
      qm_cache = QuatModel::FORCE;
   } else {
      qm_cache = QuatModel::CACHE;
   }

   d_pfmodel->computeVarGradSide(hierarchy, diffs_id, grad_id, time, qm_cache);
}
