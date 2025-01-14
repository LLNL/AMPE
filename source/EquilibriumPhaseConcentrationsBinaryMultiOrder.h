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
#ifndef included_EquilibriumPhaseConcentrationsBinaryMultiOrder
#define included_EquilibriumPhaseConcentrationsBinaryMultiOrder

#include "PhaseConcentrationsStrategy.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"

#include <string>

class EquilibriumPhaseConcentrationsBinaryMultiOrder
    : public PhaseConcentrationsStrategy
{
 public:
   EquilibriumPhaseConcentrationsBinaryMultiOrder(
       const int conc_l_id, const int conc_a_id,
       std::shared_ptr<tbox::Database> conc_db);

   virtual ~EquilibriumPhaseConcentrationsBinaryMultiOrder() {}

   virtual int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch);

   virtual int computePhaseConcentrations(const double t, double* c,
                                          double* hphi, double* x) = 0;

 private:
};

#endif
