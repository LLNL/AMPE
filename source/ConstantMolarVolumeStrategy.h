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
#ifndef included_ConstantMolarVolumeStrategy
#define included_ConstantMolarVolumeStrategy

#include "MolarVolumeStrategy.h"
#include "Phases.h"

class ConstantMolarVolumeStrategy : public MolarVolumeStrategy
{
 public:
   ConstantMolarVolumeStrategy(const double vml, const double vma,
                               const double vmb)
   {
      assert(vml == vml);
      assert(vml > 0.);
      assert(vma > 0.);

      d_inv_vm[static_cast<int>(Thermo4PFM::PhaseIndex::phaseL)] = 1.e-6 / vml;
      d_inv_vm[static_cast<int>(Thermo4PFM::PhaseIndex::phaseA)] = 1.e-6 / vma;
      if (vmb > 0.)
         d_inv_vm[static_cast<int>(Thermo4PFM::PhaseIndex::phaseB)] =
             1.e-6 / vmb;
   }


   virtual double computeInvMolarVolume(const double temperature,
                                        const double* const conc,
                                        const Thermo4PFM::PhaseIndex pi)
   {
      (void)temperature;
      (void)conc;
      return d_inv_vm[static_cast<int>(pi)];
   }

 private:
   double d_inv_vm[3];
};

#endif
