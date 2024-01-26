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
#include "ConcFACOps.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>


ConcFACOps::ConcFACOps(const std::string& object_name, const int depth,
                       const std::shared_ptr<tbox::Database>& database)
    : EllipticFACOps(tbox::Dimension(NDIM), object_name, database, depth)
{
   t_set_op_coef = tbox::TimerManager::getManager()->getTimer(
       "AMPE::ConcFACOps::setOperatorCoefficients()");
}

void ConcFACOps::setOperatorCoefficients(const double gamma,
                                         const std::vector<int>& diffusion_id,
                                         const double mobility)
{
   assert(gamma >= 0.);
   assert(diffusion_id.size() == d_d_id.size());
   assert(diffusion_id.size() == 1 || diffusion_id.size() == 2);

   t_set_op_coef->start();

   for (unsigned ic = 0; ic < diffusion_id.size(); ic++)
      assert(diffusion_id[ic] >= 0);
   for (unsigned ic = 0; ic < d_d_id.size(); ic++)
      assert(d_d_id[ic] >= 0);
   assert(mobility > 0.);

   for (unsigned ic = 0; ic < diffusion_id.size(); ic++) {
      d_hopsside->scale(d_d_id[ic], -gamma, diffusion_id[ic], true);
#ifdef DEBUG_CHECK_ASSERTIONS
      double vmax = d_hopsside->max(diffusion_id[ic]);
      double vmin = d_hopsside->min(diffusion_id[ic]);
      if (vmax <= 0.)
         tbox::pout << "Component " << ic << ", Max. for D = " << vmax
                    << ", Min. for D = " << vmin << std::endl;
      assert(vmax > 0.);
#endif
      setDPatchDataId(d_d_id[ic], ic);
      setCConstant(1., ic);
   }

   setMConstant(mobility);

   t_set_op_coef->stop();
}
