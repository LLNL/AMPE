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
#include "KimMobilityStrategyInfMob.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "KKStools.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"

#include <iomanip>

template <class FreeEnergyType>
KimMobilityStrategyInfMob<FreeEnergyType>::KimMobilityStrategyInfMob(
    const QuatModelParameters& model_parameters, QuatModel* quat_model,
    const int conc_l_id, const int conc_s_id, const int temp_id,
    const double epsilon, const double phase_well_scale,
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions,
    const double DL, const double Q0, const double mv)
    : KimMobilityStrategy<FreeEnergyType>(quat_model, conc_l_id, conc_s_id, -1,
                                          temp_id, energy_interp_func_type,
                                          conc_interp_func_type, conc_db,
                                          ncompositions),
      d_model_parameters(model_parameters),
      d_DL(DL),
      d_Q0(Q0),
      d_mv(mv)
{
   assert(epsilon > 0.);
   assert(phase_well_scale >= 0.);
   assert(mv > 0.);

   d_factor =
       kks_mobility_factor(energy_interp_func_type, epsilon, phase_well_scale);

   d_d2fdc2.resize(this->d_ncompositions * this->d_ncompositions);
}

template <class FreeEnergyType>
double KimMobilityStrategyInfMob<FreeEnergyType>::evaluateMobility(
    const double temp, const std::vector<double>& phaseconc,
    const std::vector<double>& phi)
{
   assert(d_mv > 0.);

   const Thermo4PFM::PhaseIndex pi0 = Thermo4PFM::PhaseIndex::phaseL;

   this->d_fenergy->computeSecondDerivativeFreeEnergy(temp, &phaseconc[0], pi0,
                                                      d_d2fdc2.data());

   // std::cout<<std::setprecision(15);
   // std::cout<<"c="<<phaseconc[0]<<", d2fdc2="<<d_d2fdc2[0]<<std::endl;
   const double* const cl = &phaseconc[0];
   const double* const cs = &phaseconc[this->d_ncompositions];

   double zeta_factor = d_model_parameters.zetaFactor(temp);
   if (zeta_factor < 0.) {
      assert(temp > 0.);
      assert(d_DL > 0.);
      const double DL = d_DL * exp(-d_Q0 / (gas_constant * temp));
      double zeta = 0.;
      for (unsigned i = 0; i < this->d_ncompositions; i++)
         for (unsigned j = 0; j < this->d_ncompositions; j++)
            zeta += (cl[i] - cs[i]) * d_d2fdc2[2 * i + j] * (cl[j] - cs[j]);
      zeta_factor = zeta / DL;
   }
   // convert from J/mol to pJ/um^3
   zeta_factor *= (1.e-6 / d_mv);

   const double mob = d_factor / zeta_factor;
   // std::cout<<"mob="<<mob<<", inv_zeta="<<inv_zeta_factor<<std::endl;
   assert(!std::isnan(mob));
   return mob;
}

template class KimMobilityStrategyInfMob<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>;
template class KimMobilityStrategyInfMob<
    Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>;
template class KimMobilityStrategyInfMob<
    Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary>;
template class KimMobilityStrategyInfMob<
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary2Ph1Sl>;
