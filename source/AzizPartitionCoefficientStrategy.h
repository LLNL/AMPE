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
#ifndef included_AzizPartitionCoefficientStrategy
#define included_AzizPartitionCoefficientStrategy

#include "PartitionCoefficientStrategy.h"
#include "InterpolationType.h"

#include <memory>

template <class FreeEnergyType>
class AzizPartitionCoefficientStrategy : public PartitionCoefficientStrategy
{
 public:
   AzizPartitionCoefficientStrategy(
       const int velocity_id, const int temperature_id,
       const int partition_coeff_id, const double vd, const double keq,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       std::shared_ptr<tbox::Database> conc_db);

   ~AzizPartitionCoefficientStrategy() {}

 protected:
   void evaluate(hier::Patch& patch,
                 std::shared_ptr<pdat::CellData<double> > velocity,
                 std::shared_ptr<pdat::CellData<double> > temperature,
                 std::shared_ptr<pdat::CellData<double> > partition_coeff);

 private:
   std::unique_ptr<FreeEnergyType> d_fenergy;

   // inverse of diffusion speed corresponding to interface
   double d_inv_vd;

   // pre-defined equilibrium partition coefficient
   // (if available). Negative value if needs to be computed
   double d_keq;

   double computeKeq(const double temperature);
};

#endif
