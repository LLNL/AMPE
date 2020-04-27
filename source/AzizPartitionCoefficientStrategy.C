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
#include "AzizPartitionCoefficientStrategy.h"

void AzizPartitionCoefficientStrategy::evaluate(
    hier::Patch& patch, std::shared_ptr<pdat::CellData<double> > cd_velocity,
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_partition_coeff)
{
   const hier::Box& patch_box = patch.getBox();

   pdat::CellIterator iend(pdat::CellGeometry::end(patch_box));
   for (pdat::CellIterator i(pdat::CellGeometry::begin(patch_box)); i != iend;
        ++i) {
      const pdat::CellIndex ccell = *i;

      double temperature = (*cd_temperature)(ccell);
      double vel = fabs((*cd_velocity)(ccell));
      assert(vel == vel);

      const double ke = computeKeq(temperature);
      const double scaledv = vel * d_inv_vd;
      if (ke < 1.)
         (*cd_partition_coeff)(ccell) = (ke + scaledv) / (1. + scaledv);
      else
         (*cd_partition_coeff)(ccell) = (1. + scaledv) / (1. / ke + scaledv);

      assert((*cd_partition_coeff)(ccell) > 0.);
   }
}

double AzizPartitionCoefficientStrategy::computeKeq(const double temperature)
{
   if (d_keq > 0.) return d_keq;

   assert(d_free_energy != NULL);

   double ceq[2] = {0.5, 0.5};
   bool flag = d_free_energy->computeCeqT(temperature, PhaseIndex::phaseL,
                                          PhaseIndex::phaseA, &ceq[0], 20);

   assert(flag);

   return ceq[1] / ceq[0];
}
