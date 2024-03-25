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
#include "KimMobilityStrategyInfMob3Phases.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "AMPE_internal.h"
#include "Database2JSON.h"
#include "KKStools.h"

#include <boost/property_tree/json_parser.hpp>
#include <iomanip>

namespace pt = boost::property_tree;

template <class FreeEnergyType>
KimMobilityStrategyInfMob3Phases<FreeEnergyType>::
    KimMobilityStrategyInfMob3Phases(
        const QuatModelParameters& model_parameters, QuatModel* quat_model,
        const int conc_l_id, const int conc_a_id, const int conc_b_id,
        const int temp_id, const double epsilon, const double phase_well_scale,
        const EnergyThreeArgsInterpolationType energy_three_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions,
        const double mv)
    : KimMobilityStrategy<FreeEnergyType>(
          quat_model, conc_l_id, conc_a_id, conc_b_id, temp_id,
          getTwoPhasesInterpolationType(energy_three_interp_func_type),
          conc_interp_func_type, conc_db, ncompositions),
      d_model_parameters(model_parameters),
      d_mv(mv)
{
   assert(epsilon > 0.);
   assert(phase_well_scale >= 0.);
   assert(mv > 0.);
   assert(conc_b_id > -1);
   assert(ncompositions > 0);

   Thermo4PFM::EnergyInterpolationType energy_interp_func_type;
   switch (energy_three_interp_func_type) {
      // FolchPlapp2005 reduces to PBG for two phases
      case EnergyThreeArgsInterpolationType::FOLCHPLAPP2005:
         energy_interp_func_type = Thermo4PFM::EnergyInterpolationType::PBG;
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

   d_factor =
       kks_mobility_factor(energy_interp_func_type, epsilon, phase_well_scale);
}

template <class FreeEnergyType>
double KimMobilityStrategyInfMob3Phases<FreeEnergyType>::evaluateMobility(
    const double temp, const std::vector<double>& phaseconc,
    const std::vector<double>& phi)
{
   // return 1.5e5;
   assert(phi.size() == 3);
   assert(phaseconc.size() > 0);
   assert(this->d_ncompositions == 1);

   const double phiL = phi[0];
   const double phiA = phi[1];
   const double phiB = phi[2];

   double zeta_factor = d_model_parameters.zetaFactorLA(temp);
   assert(zeta_factor > 0.);

   // convert from J/mol to pJ/um^3
   zeta_factor *= (1.e-6 / d_mv);

   const double mobLA = d_factor / zeta_factor;
   // std::cout<<"DL="<<DL<<", zeta="<<zeta<<std::endl;
   assert(mobLA == mobLA);

   zeta_factor = d_model_parameters.zetaFactorLB(temp);
   assert(zeta_factor > 0.);

   // convert from J/mol to pJ/um^3
   zeta_factor *= (1.e-6 / d_mv);

   const double mobLB = d_factor / zeta_factor;
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
    Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>;
