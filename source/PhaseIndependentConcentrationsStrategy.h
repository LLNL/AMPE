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
#ifndef included_PhaseIndependentConcentrationsStrategy
#define included_PhaseIndependentConcentrationsStrategy

#include "PhaseConcentrationsStrategy.h"

/*
 * Simply set c_l and c_s equal to c
 */
class PhaseIndependentConcentrationsStrategy
    : public PhaseConcentrationsStrategy
{
 public:
   PhaseIndependentConcentrationsStrategy(const int conc_l_id,
                                          const int conc_a_id,
                                          const int conc_b_id)
       : PhaseConcentrationsStrategy(conc_l_id, conc_a_id, conc_b_id){};

   ~PhaseIndependentConcentrationsStrategy(){};

   virtual int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch);
};

#endif
