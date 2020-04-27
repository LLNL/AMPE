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
#include "PartitionCoefficientStrategy.h"

void PartitionCoefficientStrategy::evaluate(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   assert(d_partition_coeff_id >= 0);

   const int maxl = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {

         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double> > partition_coeff(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_partition_coeff_id)));

         std::shared_ptr<pdat::CellData<double> > velocity;
         if (d_velocity_id >= 0) {
            velocity = std::dynamic_pointer_cast<pdat::CellData<double>,
                                                   hier::PatchData>(
                patch->getPatchData(d_velocity_id));
            assert(velocity);
         }

         std::shared_ptr<pdat::CellData<double> > temperature;
         if (d_temperature_id >= 0) {
            temperature = std::dynamic_pointer_cast<pdat::CellData<double>,
                                                      hier::PatchData>(
                patch->getPatchData(d_temperature_id));
            assert(temperature);
         }

         evaluate(*patch, velocity, temperature, partition_coeff);
      }
   }
}
