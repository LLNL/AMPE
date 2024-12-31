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
#ifndef included_EquilibriumPhaseConcentrationsBinary
#define included_EquilibriumPhaseConcentrationsBinary

#include "PhaseConcentrationsStrategy.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "InterpolationType.h"

#include "SAMRAI/tbox/InputManager.h"

/*!
 * The main role for this class is to implement the function
 * computePhaseConcentrationsOnPatch()
 * that solves the KKS problem for a binary two phases alloy
 * at every cell of a patch
 */
class EquilibriumPhaseConcentrationsBinary : public PhaseConcentrationsStrategy
{
 public:
   EquilibriumPhaseConcentrationsBinary(
       const int conc_l_id, const int conc_a_id,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type);

   virtual ~EquilibriumPhaseConcentrationsBinary() {}

   int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch);

 private:
   const Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;

   virtual int computePhaseConcentrations(const double temperature,
                                          const double conc, const double hphi,
                                          double* sol) = 0;
};

#endif
