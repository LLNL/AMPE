// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
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

void SimpleGradStrategy::computeDiffs(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& var_id,
    int& diffs_id, const double time)
{
   assert(d_pfmodel != nullptr);

   d_pfmodel->computeVarDiffs(patch_level, var_id, diffs_id, time);
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

void SimpleGradStrategy::computeGradCell(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
    int& grad_id, const double time)
{
   assert(d_pfmodel != nullptr);

   d_pfmodel->computeVarGradCell(patch_level, diffs_id, grad_id, time);
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

//-----------------------------------------------------------------------

void SimpleGradStrategy::computeGradSide(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
    int& grad_id, const double time)
{
   assert(d_pfmodel != nullptr);

   d_pfmodel->computeVarGradSide(patch_level, diffs_id, grad_id, time);
}

//-----------------------------------------------------------------------
