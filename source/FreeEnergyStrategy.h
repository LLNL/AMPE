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
#ifndef included_FreeEnergyStrategy
#define included_FreeEnergyStrategy

#include "PhysicalConstants.h"
#include "Phases.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"

#include <vector>

using namespace SAMRAI;
#ifdef HAVE_THERMO4PFM
using namespace Thermo4PFM;
#else
using namespace ampe_thermo;
#endif

class FreeEnergyStrategy
{
 public:
   FreeEnergyStrategy(){};
   virtual ~FreeEnergyStrategy(){};

   // mesh functions
   virtual void computeFreeEnergyLiquid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_l_id, const bool gp);

   virtual void computeFreeEnergySolidA(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_a_id, const bool gp);

   virtual void computeFreeEnergySolidB(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int f_b_id, const bool gp);

   // mesh functions
   virtual void computeFreeEnergyLiquid(hier::Patch& patch,
                                        const int temperature_id,
                                        const int f_l_id, const bool gp) = 0;

   virtual void computeFreeEnergySolidA(hier::Patch& patch,
                                        const int temperature_id,
                                        const int f_a_id, const bool gp) = 0;

   virtual void computeFreeEnergySolidB(hier::Patch& patch,
                                        const int temperature_id,
                                        const int f_b_id, const bool gp) = 0;

   virtual void addDrivingForce(const double time, hier::Patch& patch,
                                const int temperature_id, const int phase_id,
                                const int eta_id, const int conc_id,
                                const int f_l_id, const int f_a_id,
                                const int f_b_id, const int rhs_id) = 0;

   virtual void computeDrivingForce(
       const double time, const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id, const int eta_id,
       const int conc_id, const int f_l_id, const int f_a_id, const int f_b_id,
       const int rhs_id);

   virtual void computeDrivingForce(const double time, hier::Patch& patch,
                                    const int temperature_id,
                                    const int phase_id, const int eta_id,
                                    const int conc_id, const int f_l_id,
                                    const int f_a_id, const int f_b_id,
                                    const int rhs_id)
   {
      std::shared_ptr<pdat::CellData<double> > rhs(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(rhs_id)));
      math::PatchCellDataOpsReal<double> ops;
      ops.setToScalar(rhs, 0., patch.getBox());

      addDrivingForce(time, patch, temperature_id, phase_id, eta_id, conc_id,
                      f_l_id, f_a_id, f_b_id, rhs_id);
   };

   // pointwise functions
   virtual void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true) = 0;
   virtual void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true) = 0;
#ifndef HAVE_THERMO4PFM
   virtual void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true) = 0;
#endif
   virtual void preRunDiagnostics(const double temperature) = 0;

 private:
   void computeDerivFreeEnergy(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int df_id);
};
#endif
