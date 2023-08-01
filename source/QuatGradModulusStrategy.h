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
#ifndef included_QuatGradModulusStrategy
#define included_QuatGradModulusStrategy

#include <memory>

#include "SAMRAI/hier/PatchLevel.h"

using namespace SAMRAI;

class QuatGradModulusStrategy
{
 public:
   QuatGradModulusStrategy(const int qlen) : d_qlen(qlen) {}

   void computeQuatGradModulusFromSides(
       const std::shared_ptr<hier::PatchLevel> level, int& grad_side_id,
       int& grad_modulus_id);

   void computeQuatGradModulus(const std::shared_ptr<hier::PatchLevel> level,
                               int& grad_cell_id, int& grad_modulus_id);

 private:
   int d_qlen;
};

#endif
