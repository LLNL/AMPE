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
#ifndef included_ConcFreeEnergyStrategy
#define included_ConcFreeEnergyStrategy

#include "FreeEnergyStrategy.h"

class ConcFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   ConcFreeEnergyStrategy(const int conc_l_id, const int conc_a_id,
                          const int conc_b_id);

   virtual ~ConcFreeEnergyStrategy(){};

   virtual bool computeCeqT(const double temperature,
                            const Thermo4PFM::PhaseIndex pi0,
                            const Thermo4PFM::PhaseIndex pi1, double* ceq)
   {
      return false;
   };

   // generic loop over levels and patches
   virtual void computeDerivFreeEnergyLiquid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_l_id);

   virtual void computeDerivFreeEnergySolidA(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_a_id);

   virtual void computeDerivFreeEnergySolidB(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_b_id);

   // implementation of mesh functions required by base class FreeEnergyStrategy
   void computeDerivFreeEnergyLiquid(hier::Patch& patch,
                                     const int temperature_id,
                                     const int f_l_id);

   void computeDerivFreeEnergySolidA(hier::Patch& patch,
                                     const int temperature_id,
                                     const int f_a_id);

   void computeDerivFreeEnergySolidB(hier::Patch& patch,
                                     const int temperature_id,
                                     const int f_b_id);

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp);
   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp);
   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp);

   // pure virtual functions to be implemented in derived classes
   virtual void computeFreeEnergy(hier::Patch& patch, const int temperature_id,
                                  const int f_id, const int c_i_id,
                                  const Thermo4PFM::PhaseIndex pi,
                                  const bool gp) = 0;

   virtual void computeDerivFreeEnergy(hier::Patch& patch,
                                       const int temperature_id,
                                       const int dfl_id, const int conc_id,
                                       Thermo4PFM::PhaseIndex) = 0;

 protected:
   const int d_conc_l_id;
   const int d_conc_a_id;
   const int d_conc_b_id;
};

#endif
