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
#ifndef included_QuadraticEquilibriumThreePhasesTernaryMultiOrder
#define included_QuadraticEquilibriumThreePhasesTernaryMultiOrder

#include "PhaseConcentrationsStrategy.h"
#include "QuatModelParameters.h"
#include "QuadraticFreeEnergyFunctionsTernaryThreePhase.h"

#include <string>

class QuadraticEquilibriumThreePhasesTernaryMultiOrder
    : public PhaseConcentrationsStrategy
{
 public:
   QuadraticEquilibriumThreePhasesTernaryMultiOrder(
       const short norderp_A, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const QuatModelParameters& model_parameters,
       std::shared_ptr<tbox::Database> conc_db);

   ~QuadraticEquilibriumThreePhasesTernaryMultiOrder() {}

   virtual int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch);

 private:
   // number of order parameters associated with phase A
   const short d_norderp_A;

   Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;

   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsTernaryThreePhase>
       d_fenergy;
};

#endif
