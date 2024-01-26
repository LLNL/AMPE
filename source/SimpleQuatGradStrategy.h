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
#ifndef included_SimpleQuatGradStrategy
#define included_SimpleQuatGradStrategy

#include "QuatGradStrategy.h"
#include "QuatModel.h"

class SimpleQuatGradStrategy : public QuatGradStrategy
{
 public:
   SimpleQuatGradStrategy(QuatModel* quat_model);

   ~SimpleQuatGradStrategy();

   virtual bool isSymmetryAware(void);

   virtual void computeDiffs(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& quat_id,
       int& diffs_id, const double time);

   virtual void computeDiffs(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& quat_id,
       int& diffs_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeGradCell(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_id, const double time);

   virtual void computeGradCell(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeGradSide(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& diffs_id,
       int& grad_id, const double time);

   virtual void computeGradSide(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
       int& grad_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeGradModulus(
       const std::shared_ptr<hier::PatchLevel> patch_level, int& grad_id,
       int& mod_id, const double time);

   virtual void computeGradModulus(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_id,
       int& mod_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeGradModulusFromSides(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_id,
       int& mod_id, const double time, const CACHE_TYPE cache = CACHE);

 private:
   QuatModel* d_quat_model;

   enum QuatModel::CACHE_TYPE translateCacheType(const CACHE_TYPE cache);
};

#endif
