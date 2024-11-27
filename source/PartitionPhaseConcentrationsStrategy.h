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
#ifndef included_PartitionPhaseConcentrationsStrategy
#define included_PartitionPhaseConcentrationsStrategy

#include "PhaseConcentrationsStrategy.h"
#include "InterpolationType.h"

// compute c_l, c_s using partition coefficient
class PartitionPhaseConcentrationsStrategy : public PhaseConcentrationsStrategy
{
 public:
   PartitionPhaseConcentrationsStrategy(
       const int conc_l_id, const int conc_a_id, const int conc_b_id,
       const Thermo4PFM::ConcInterpolationType phase_interp_func_type,
       const int partition_coeff_id)
       : PhaseConcentrationsStrategy(conc_l_id, conc_a_id, conc_b_id),
         d_phase_interp_func_type(phase_interp_func_type),
         d_partition_coeff_id(partition_coeff_id)
   {
      assert(d_partition_coeff_id >= 0);
      assert(phase_interp_func_type !=
             Thermo4PFM::ConcInterpolationType::UNDEFINED);
   };

   ~PartitionPhaseConcentrationsStrategy(){};

   virtual int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch);

 private:
   const Thermo4PFM::ConcInterpolationType d_phase_interp_func_type;
   int d_partition_coeff_id;
};

#endif
