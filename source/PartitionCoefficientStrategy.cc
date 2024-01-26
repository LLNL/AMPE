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
