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
#ifndef included_HierarchyStencilOps
#define included_HierarchyStencilOps

#include "SAMRAI/hier/PatchHierarchy.h"

#include <memory>

using namespace SAMRAI;

class HierarchyStencilOps
{
 public:
   void computeGradSide(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int diffs_id, int grad_side_id);

   void computeDiffs(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                     int var_id, int diffs_id);
   void computeDiffs(std::shared_ptr<hier::PatchLevel> level, int var_id,
                     int diffs_id);
   void computeDiffs(std::shared_ptr<hier::Patch> patch, int var_id,
                     int diffs_id);

   void computeGradCell(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                        int diffs_id, int grad_id);
   void computeGradCell(std::shared_ptr<hier::PatchLevel> level, int diffs_id,
                        int grad_id);

 private:
   void computeGradSide(std::shared_ptr<hier::Patch> patch, int diffs_id,
                        int grad_side_id);

   void computeGradCell(std::shared_ptr<hier::Patch> patch, int diffs_id,
                        int grad_id);
};

#endif
