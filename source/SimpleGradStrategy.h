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
#ifndef included_SimpleGradStrategy
#define included_SimpleGradStrategy

#include "GradStrategy.h"
#include "QuatModel.h"

class SimpleGradStrategy : public GradStrategy
{
 public:
   SimpleGradStrategy(QuatModel* model);

   ~SimpleGradStrategy();

   virtual void computeDiffs(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& var_id,
       int& diffs_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeGradCell(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& var_id,
       int& grad_id, const double time, const CACHE_TYPE cache = CACHE);

   virtual void computeGradSide(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& var_id,
       int& grad_id, const double time, const CACHE_TYPE cache = CACHE);

 private:
   QuatModel* d_pfmodel;
};

#endif
