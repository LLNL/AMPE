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
#ifndef included_QuatGradStrategy
#define included_QuatGradStrategy

#include "SAMRAI/hier/PatchHierarchy.h"

#include <cassert>

using namespace SAMRAI;

class QuatGradStrategy
{
 public:
   enum CACHE_TYPE { CACHE = 0, FORCE = 1 };

   QuatGradStrategy(){};

   virtual ~QuatGradStrategy(){};

   virtual bool isSymmetryAware(void) = 0;

   virtual void computeDiffs(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& quat_id,
       int& diffs_id, const double time) = 0;

   virtual void computeDiffs(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& quat_id,
       int& diffs_id, const double time, const CACHE_TYPE cache = CACHE) = 0;

   virtual void computeGradCell(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_id, const double time) = 0;

   virtual void computeGradCell(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_id, const double time, const CACHE_TYPE cache = CACHE) = 0;

   virtual void computeGradSide(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_id, const double time) = 0;

   virtual void computeGradSide(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_id, const double time, const CACHE_TYPE cache = CACHE) = 0;

   virtual void computeGradModulus(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& grad_id,
       int& mod_id, const double time)
   {
      (void)patch_level;
      (void)grad_id;
      (void)mod_id;
      (void)time;
      std::cerr << "computeGradModulus unimplemented" << std::endl;
   }

   virtual void computeGradModulus(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_id,
       int& mod_id, const double time, const CACHE_TYPE cache = CACHE) = 0;

   virtual void computeGradModulusFromSides(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_id,
       int& mod_id, const double time, const CACHE_TYPE cache = CACHE) = 0;
};

#endif
