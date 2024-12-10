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
#ifndef included_EquilibriumPhaseConcentrationsThreePhases
#define included_EquilibriumPhaseConcentrationsThreePhases

#include "PhaseConcentrationsStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/tbox/InputManager.h"

#include <boost/property_tree/ptree.hpp>

template <class FreeEnergyType>
class EquilibriumPhaseConcentrationsThreePhases
    : public PhaseConcentrationsStrategy
{
 public:
   EquilibriumPhaseConcentrationsThreePhases(
       const int conc_l_id, const int conc_a_id, const int conc_b_id,
       const int conc_l_ref_id, const int conc_a_ref_id,
       const int conc_b_ref_id,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       boost::property_tree::ptree calphad_pt,
       std::shared_ptr<tbox::Database> newton_db, const unsigned ncompositions);

   ~EquilibriumPhaseConcentrationsThreePhases() {}

   virtual int computePhaseConcentrationsOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<hier::Patch> patch);

 protected:
   std::unique_ptr<FreeEnergyType> d_calphad_fenergy;

   virtual int computeAuxilliaryConcentrations(const double temp, double* c,
                                               double* hphi, double* x);

 private:
   int d_conc_l_ref_id;
   int d_conc_a_ref_id;
   int d_conc_b_ref_id;

   const Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;
};

#endif
