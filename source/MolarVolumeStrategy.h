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
#ifndef included_MolarVolumeStrategy
#define included_MolarVolumeStrategy

#ifdef HAVE_THERMO4PFM
using namespace Thermo4PFM;
#else
using namespace ampe_thermo;
#endif

class MolarVolumeStrategy
{
 public:
   virtual double computeInvMolarVolume(const double temperature,
                                        const double* const conc,
                                        const PhaseIndex pi) = 0;
};

#endif
