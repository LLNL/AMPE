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
//
#ifndef included_ConcFreeEnergyStrategy
#define included_ConcFreeEnergyStrategy

#include "FreeEnergyStrategy.h"

class ConcFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
                            const PhaseIndex pi1, double* ceq) = 0;

   virtual void energyVsPhiAndC(const double temperature,
                                const double* const ceq, const bool found_ceq,
                                const double phi_well_scale,
                                const std::string& phi_well_type,
                                const int npts_phi, const int npts_c) = 0;

   virtual void computeDerivFreeEnergyLiquid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_l_id);

   virtual void computeDerivFreeEnergySolidA(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_a_id);

   virtual void computeDerivFreeEnergySolidB(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_b_id);

   virtual void computeDerivFreeEnergyLiquid(hier::Patch& patch,
                                             const int temperature_id,
                                             const int f_l_id) = 0;

   virtual void computeDerivFreeEnergySolidA(hier::Patch& patch,
                                             const int temperature_id,
                                             const int f_a_id) = 0;

   virtual void computeDerivFreeEnergySolidB(hier::Patch& patch,
                                             const int temperature_id,
                                             const int f_b_id) = 0;
};

#endif
