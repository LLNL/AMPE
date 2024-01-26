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

#ifdef HAVE_THERMO4PFM
using namespace Thermo4PFM;
#else
using namespace ampe_thermo;
#endif

class PhaseFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   PhaseFreeEnergyStrategy(const EnergyInterpolationType phase_interp_func_type,
                           const EnergyInterpolationType eta_interp_func_type,
                           const double fl, const double fa, const double fb,
                           const double vml, const double vma, const double vmb,
                           const bool with_third_phase);

   ~PhaseFreeEnergyStrategy(){};

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp);

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp);

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp);

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int eta_id, const int conc_id, const int fl_id,
                        const int fa_id, const int fb_id, const int rhs_id);

   void addDrivingForceEta(const double time, hier::Patch& patch,
                           const int temperature_id, const int phase_id,
                           const int eta_id, const int conc_id, const int fl_id,
                           const int fa_id, const int fb_id, const int rhs_id);

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

   void preRunDiagnostics(const double temperature){};

 private:
   double d_f_l;
   double d_f_a;
   double d_f_b;
   EnergyInterpolationType d_phase_interp_func_type;
   EnergyInterpolationType d_eta_interp_func_type;

   bool d_with_third_phase;
};

#endif
