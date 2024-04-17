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
#ifndef included_EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases
#define included_EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases

#include "PhaseConcentrationsStrategy.h"
#include "QuatModelParameters.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"

#include <string>

class EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases
    : public PhaseConcentrationsStrategy
{
 public:
   EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases(
       const short norderp_A, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const QuatModelParameters& model_parameters,
       std::shared_ptr<tbox::Database> conc_db);

   virtual ~EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases() {}

   virtual int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch);

   virtual int computePhaseConcentrations(const double t, double* c,
                                          double* hphi, double* x) = 0;

 private:
   // number of order parameters associated with phase A
   const short d_norderp_A;

   Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;
};

#endif
