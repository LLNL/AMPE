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
#ifndef included_DeltaTemperatureFreeEnergyStrategy
#define included_DeltaTemperatureFreeEnergyStrategy

#include "FreeEnergyStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cstring>

class DeltaTemperatureFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   DeltaTemperatureFreeEnergyStrategy(
       const double Tm, const double latentHeat,
       const Thermo4PFM::EnergyInterpolationType phase_interp_func_type);

   virtual ~DeltaTemperatureFreeEnergyStrategy(){};

   virtual void addDrivingForce(const double time, hier::Patch& patch,
                                const int temperature_id, const int phase_id,
                                const int conc_id, const int fl_id,
                                const int fa_id, const int fb_id,
                                const int rhs_id);

   void applydPhidTBlock(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                         const int temperature_id, const int phase_id,
                         const int rhs_id, const double phase_mobility);

   virtual void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   virtual void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   virtual void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp);

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fa_id, const bool gp);

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fb_id, const bool gp);

   void preRunDiagnostics(const double temperature){};

 private:
   // melting temperature
   const double d_Tm;

   // latent heat
   const double d_L;

   const Thermo4PFM::EnergyInterpolationType d_phase_interp_func_type;
};

#endif
