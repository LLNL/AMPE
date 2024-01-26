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
