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
#include "UniformPartitionCoefficientStrategy.h"

void UniformPartitionCoefficientStrategy::evaluate(
    hier::Patch& patch, std::shared_ptr<pdat::CellData<double> > cd_velocity,
    std::shared_ptr<pdat::CellData<double> > cd_temperature,
    std::shared_ptr<pdat::CellData<double> > cd_partition_coeff)
{
   (void)patch;
   (void)cd_velocity;
   (void)cd_temperature;

   cd_partition_coeff->fill(d_keq);
}
