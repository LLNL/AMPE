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
#ifndef included_PartitionCoefficientStrategy
#define included_PartitionCoefficientStrategy

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"

using namespace SAMRAI;


class PartitionCoefficientStrategy
{
 public:
   PartitionCoefficientStrategy(const int velocity_id, const int temperature_id,
                                const int partition_coeff_id)
       : d_velocity_id(velocity_id),
         d_temperature_id(temperature_id),
         d_partition_coeff_id(partition_coeff_id)
   {
      assert(d_partition_coeff_id >= 0);
   };

   virtual ~PartitionCoefficientStrategy(){};

   virtual void evaluate(
       hier::Patch& patch, std::shared_ptr<pdat::CellData<double> > cd_velocity,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_partition_coeff) = 0;

   void evaluate(const std::shared_ptr<hier::PatchHierarchy> hierarchy);

 private:
   int d_velocity_id;
   int d_temperature_id;
   int d_partition_coeff_id;
};

#endif
