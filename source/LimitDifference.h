// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_LimitDifference
#define included_LimitDifference

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"

class LimitDifference
{
 public:
   LimitDifference(
       const std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy);

   /*!
    * Change values in cell_data so that they differ from ref_cell_data
    * by less than "diff"
    */
   void apply(const int cell_data_id, const int ref_cell_data_id,
              const int wright_id, const int depth, const double diff);

 private:
   const std::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;

   void apply(const std::shared_ptr<SAMRAI::hier::PatchLevel> level,
              const int cell_data_id, const int ref_cell_data_id,
              const int wright_id, const int depth, const double diff);
};

#endif
