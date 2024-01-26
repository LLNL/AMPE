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
#ifndef included_UniformPartitionCoefficientStrategy
#define included_UniformPartitionCoefficientStrategy

#include "PartitionCoefficientStrategy.h"

class UniformPartitionCoefficientStrategy : public PartitionCoefficientStrategy
{
 public:
   UniformPartitionCoefficientStrategy(const int velocity_id,
                                       const int temperature_id,
                                       const int partition_coeff_id,
                                       const double keq)
       : PartitionCoefficientStrategy(velocity_id, temperature_id,
                                      partition_coeff_id),
         d_keq(keq)
   {
      tbox::plog << "UniformPartitionCoefficientStrategy with keq=" << keq
                 << std::endl;
   }

 protected:
   void evaluate(hier::Patch& patch,
                 std::shared_ptr<pdat::CellData<double> > velocity,
                 std::shared_ptr<pdat::CellData<double> > temperature,
                 std::shared_ptr<pdat::CellData<double> > partition_coeff);

 private:
   // pre-defined equilibrium partition coefficient
   double d_keq;
};

#endif
