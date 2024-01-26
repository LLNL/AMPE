// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_ApplyPolynomial
#define included_ApplyPolynomial

#include "InterpolationType.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"

using namespace SAMRAI;
#ifdef HAVE_THERMO4PFM
using namespace Thermo4PFM;
#else
using namespace ampe_thermo;
#endif

class ApplyPolynomial
{
 public:
   ApplyPolynomial(EnergyInterpolationType interp_func_type);

   void apply(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
              const int src_cell_data_id, const int dst_cell_data_id);

 private:
   EnergyInterpolationType d_interp_func_type;

   void apply(const std::shared_ptr<hier::PatchLevel> level,
              const int src_cell_data_id, const int dst_cell_data_id);
};

#endif
