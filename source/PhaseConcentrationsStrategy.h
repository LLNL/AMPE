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
#ifndef included_PhaseConcentrationsStrategy
#define included_PhaseConcentrationsStrategy

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchHierarchy.h"

using namespace SAMRAI;

class PhaseConcentrationsStrategy
{
 public:
   PhaseConcentrationsStrategy(const int conc_l_id, const int conc_a_id,
                               const int conc_b_id);

   virtual ~PhaseConcentrationsStrategy() {}

   void computePhaseConcentrations(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id,
       const int concentration_id);

   virtual int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch) = 0;

 protected:
   int d_conc_l_id;
   int d_conc_a_id;
   int d_conc_b_id;
};

#endif
