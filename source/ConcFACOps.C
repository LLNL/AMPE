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
#include "ConcFACOps.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>
using namespace std;

ConcFACOps::ConcFACOps(const std::string& object_name, const int depth,
                       const boost::shared_ptr<tbox::Database>& database)
    : EllipticFACOps(tbox::Dimension(NDIM), object_name, database, depth)
{
   t_set_op_coef = tbox::TimerManager::getManager()->getTimer(
       "AMPE::ConcFACOps::setOperatorCoefficients()");
}

void ConcFACOps::setOperatorCoefficients(const double gamma,
                                         const vector<int>& diffusion_id,
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
      d_hopsside->scale(d_d_id[ic], -gamma, diffusion_id[ic]);
#ifdef DEBUG_CHECK_ASSERTIONS
      double vmax = d_hopsside->max(diffusion_id[ic]);
      double vmin = d_hopsside->min(diffusion_id[ic]);
      if (vmax <= 0.)
         tbox::pout << "Component " << ic << ", Max. for D = " << vmax
                    << ", Min. for D = " << vmin << endl;
      assert(vmax > 0.);
#endif
      setDPatchDataId(d_d_id[ic], ic);
      setCConstant(1., ic);
   }

   setMConstant(mobility);

   t_set_op_coef->stop();
}
