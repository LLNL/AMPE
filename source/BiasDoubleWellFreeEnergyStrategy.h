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
#ifndef included_BiasDoubleWellFreeEnergyStrategy
#define included_BiasDoubleWellFreeEnergyStrategy

#include "FreeEnergyStrategy.h"
#include "MeltingTemperatureStrategy.h"

/*!
 * @brief Class BiasDoubleWellFreeEnergyStrategy implements the free energy
 * used by Kobayashi in Physica D 63 (1993), p. 410-423
 */

class BiasDoubleWellFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   BiasDoubleWellFreeEnergyStrategy();

   virtual ~BiasDoubleWellFreeEnergyStrategy(){};

   // mesh functions
   virtual void computeFreeEnergyLiquid(hier::Patch& patch,
                                        const int temperature_id,
                                        const int f_l_id, const bool gp)
   {
      (void)patch;
      (void)temperature_id;
      (void)f_l_id;
      (void)gp;
   };

   virtual void computeFreeEnergySolidA(hier::Patch& patch,
                                        const int temperature_id,
                                        const int f_a_id, const bool gp)
   {
      (void)patch;
      (void)temperature_id;
      (void)f_a_id;
      (void)gp;
   };

   virtual void computeFreeEnergySolidB(hier::Patch& patch,
                                        const int temperature_id,
                                        const int f_b_id, const bool gp)
   {
      (void)patch;
      (void)temperature_id;
      (void)f_b_id;
      (void)gp;
   };

   virtual void addDrivingForce(const double time, hier::Patch& patch,
                                const int temperature_id, const int phase_id,
                                const int conc_id, const int fl_id,
                                const int fa_id, const int fb_id,
                                const int rhs_id) = 0;

   void computePhaseConcentrations(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id,
       const int concentration_id);

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

 private:
};

#endif
