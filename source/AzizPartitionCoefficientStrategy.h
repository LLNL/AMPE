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
#ifndef included_AzizPartitionCoefficientStrategy
#define included_AzizPartitionCoefficientStrategy

#include "PartitionCoefficientStrategy.h"
#include "FreeEnergyFunctions.h"

class AzizPartitionCoefficientStrategy : public PartitionCoefficientStrategy
{
 public:
   AzizPartitionCoefficientStrategy(const int velocity_id,
                                    const int temperature_id,
                                    const int partition_coeff_id,
                                    FreeEnergyFunctions* free_energy,
                                    const double vd, const double keq)
       : PartitionCoefficientStrategy(velocity_id, temperature_id,
                                      partition_coeff_id),
         d_free_energy(free_energy),
         d_inv_vd(1. / vd),
         d_keq(keq)
   {
      assert(d_inv_vd == d_inv_vd);
   }

 protected:
   void evaluate(hier::Patch& patch,
                 std::shared_ptr<pdat::CellData<double> > velocity,
                 std::shared_ptr<pdat::CellData<double> > temperature,
                 std::shared_ptr<pdat::CellData<double> > partition_coeff);

 private:
   FreeEnergyFunctions* d_free_energy;

   // inverse of diffusion speed corresponding to interface
   double d_inv_vd;

   // pre-defined equilibrium partition coefficient
   // (if available). Negative value if needs to be computed
   double d_keq;

   double computeKeq(const double temperature);
};

#endif
