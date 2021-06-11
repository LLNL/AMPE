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
