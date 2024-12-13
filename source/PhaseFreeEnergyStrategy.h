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
#ifndef included_PhaseFreeEnergyStrategy
#define included_PhaseFreeEnergyStrategy

#include "FreeEnergyStrategy.h"
#include "InterpolationType.h"

class PhaseFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   PhaseFreeEnergyStrategy(
       const Thermo4PFM::EnergyInterpolationType phase_interp_func_type,
       const double fl, const double fa, const double vml, const double vma);

   ~PhaseFreeEnergyStrategy(){};

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp) override;

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp) override;

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp) override;

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int fl_id, const int fa_id,
                        const int fb_id, const int rhs_id) override;

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;

   void preRunDiagnostics(const double temperature) override{};

 private:
   double d_f_l;
   double d_f_a;

   Thermo4PFM::EnergyInterpolationType d_phase_interp_func_type;
};

#endif
