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
#ifdef HAVE_THERMO4PFM

#include "KimMobilityStrategyInfMob3Phases.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "AMPE_internal.h"

#include <iomanip>

template <class FreeEnergyType>
KimMobilityStrategyInfMob3Phases<FreeEnergyType>::
    KimMobilityStrategyInfMob3Phases(
        QuatModel* quat_model, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const int temp_id, const double epsilon,
        const double phase_well_scale,
        const EnergyThreeArgsInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions,
        const double DL, const double Q0, const double mv)
    : KimMobilityStrategy<FreeEnergyType>(
          quat_model, conc_l_id, conc_a_id, conc_b_id, temp_id,
          getTwoPhasesInterpolationType(energy_interp_func_type),
          conc_interp_func_type, conc_db, ncompositions),
      d_DL(DL),
      d_Q0(Q0)
{
   assert(epsilon > 0.);
   assert(phase_well_scale >= 0.);
   assert(mv > 0.);
   assert(conc_b_id > -1);
   assert(ncompositions > 0);

   double a2 = 0.;
   switch (energy_interp_func_type) {
      // FolchPlapp2005 reduces to PBG for two phases
      case EnergyThreeArgsInterpolationType::FOLCHPLAPP2005:
         a2 = 47. / 60.;
         break;
      case EnergyThreeArgsInterpolationType::MOELANS2011:
         TBOX_ERROR(
             "MOELANS2011 is not a valid interpolation option in "
             "KimMobilityStrategyInfMob3Phases");
      default:
         TBOX_ERROR(
             "Invalid interpolation function in "
             "KimMobilityStrategyInfMob3Phases");
   }

   const double xi = epsilon / sqrt(32. * phase_well_scale);

   d_factor = 3. * (2. * xi * xi) * a2;
   d_factor *= (1.e-6 / mv);  // convert from J/mol to pJ/um^3

   d_d2fdc2.resize(ncompositions * ncompositions);
}

template <class FreeEnergyType>
double KimMobilityStrategyInfMob3Phases<FreeEnergyType>::compute_zeta(
    const double* const cl, const double* const cs, const double temp)
{
   double zeta = 0.;
   for (unsigned i = 0; i < this->d_ncompositions; i++)
      for (unsigned j = 0; j < this->d_ncompositions; j++)
         zeta += (cl[i] - cs[i]) * d_d2fdc2[2 * i + j] * (cl[j] - cs[j]);
   return zeta;
}

template <class FreeEnergyType>
double KimMobilityStrategyInfMob3Phases<FreeEnergyType>::evaluateMobility(
    const double temp, const std::vector<double>& phaseconc,
    const std::vector<double>& phi)
{
   // return 1.5e5;
   assert(phi.size() == 3);
   assert(phaseconc.size() > 0);

   const double phiL = phi[0];
   const double phiA = phi[1];
   const double phiB = phi[2];

   const PhaseIndex pil = PhaseIndex::phaseL;

   // std::cout<<std::setprecision(15);
   // std::cout<<"c="<<phaseconc[0]<<", d2fdc2="<<d_d2fdc2[0]<<std::endl;
   const double* const cl = &phaseconc[0];
   const double* const ca = &phaseconc[this->d_ncompositions];
   const double* const cb = &phaseconc[2 * this->d_ncompositions];

   this->d_fenergy->computeSecondDerivativeFreeEnergy(temp, &phaseconc[0], pil,
                                                      d_d2fdc2.data());

   const double DL = d_DL * exp(-d_Q0 / (gas_constant * temp));

   // LA
   double zeta = compute_zeta(cl, ca, temp);
   const double mobLA = DL / (d_factor * zeta);
   // std::cout<<"DL="<<DL<<", zeta="<<zeta<<std::endl;
   assert(mobLA == mobLA);

   // LB
   zeta = compute_zeta(cl, cb, temp);
   const double mobLB = DL / (d_factor * zeta);
   assert(mobLB == mobLB);

   // AB
   // since Kim's formula was derived for DS << DL, it cannot be used for AB
   // so we use mobAB = 0.5*(mobLA+mobLB)
   // see Kim, Kim, Suzuki, Ode, J. Crystal Growth 2004
   const double mobAB = 0.5 * (mobLA + mobLB);

   double wLA, wLB, wAB;
   computeWeights3Pairs(phiL, phiA, phiB, wLA, wLB, wAB);
   assert(wLA <= 1.);
   assert(wLB <= 1.);
   assert(wAB <= 1.);

   return wLA * mobLA + wLB * mobLB + wAB * mobAB;
}

template class KimMobilityStrategyInfMob3Phases<
    CALPHADFreeEnergyFunctionsBinary3Ph2Sl>;

#endif
