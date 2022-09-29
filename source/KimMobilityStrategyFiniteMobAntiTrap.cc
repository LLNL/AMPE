// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "KimMobilityStrategyFiniteMobAntiTrap.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#ifdef HAVE_THERMO4PFM
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#endif

template <class FreeEnergyType>
KimMobilityStrategyFiniteMobAntiTrap<FreeEnergyType>::
    KimMobilityStrategyFiniteMobAntiTrap(
        QuatModel* quat_model, const int conc_l_id, const int conc_s_id,
        const int temp_id, const double interface_mobility,
        const double epsilon, const double phase_well_scale,
        const EnergyInterpolationType energy_interp_func_type,
        const ConcInterpolationType conc_interp_func_type,
        std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions,
        const double DL, const double Q0, const double mv)
    : KimMobilityStrategy<FreeEnergyType>(quat_model, conc_l_id, conc_s_id,
                                          temp_id, energy_interp_func_type,
                                          conc_interp_func_type, conc_db,
                                          ncompositions),
      d_DL(DL),
      d_Q0(Q0)
{
   assert(epsilon > 0.);
   assert(phase_well_scale >= 0.);
   assert(mv > 0.);
   assert(interface_mobility > 0.);

   double a2 = 0.;
   switch (energy_interp_func_type) {
      case EnergyInterpolationType::PBG: a2 = 47. / 60.; break;
      case EnergyInterpolationType::HARMONIC: a2 = 5. / 6.; break;
      default:
         TBOX_ERROR(
             "Invalid interpolation function in "
             "KimMobilityStrategyFiniteMobAntiTrap");
   }
   // a2=19./30.; //for Kim2007
   const double xi = epsilon / sqrt(16. * phase_well_scale);
   tbox::pout << "interface_mobility=" << interface_mobility << std::endl;
   d_alpha = 3. * sqrt(2.) * xi / interface_mobility;
   d_beta = 3. * a2 * xi * xi;
   d_beta *= (1.e-6 / mv);  // convert zeta from J/mol to pJ/um^3

   d_d2fdc2.resize(this->d_ncompositions * this->d_ncompositions);
}

template <class FreeEnergyType>
double KimMobilityStrategyFiniteMobAntiTrap<FreeEnergyType>::evaluateMobility(
    const double temp, const std::vector<double>& phaseconc)
{
   const PhaseIndex pi0 = PhaseIndex::phaseL;

   this->d_fenergy->computeSecondDerivativeFreeEnergy(temp, &phaseconc[0], pi0,
#ifdef HAVE_THERMO4PFM
                                                      d_d2fdc2.data()
#else
                                                      d_d2fdc2
#endif
   );

   const double* const cl = &phaseconc[0];
   const double* const cs = &phaseconc[this->d_ncompositions];

   double zeta = 0.;
   for (unsigned i = 0; i < this->d_ncompositions; i++)
      for (unsigned j = 0; j < this->d_ncompositions; j++)
         zeta += (cl[i] - cs[i]) * d_d2fdc2[2 * i + j] * (cl[j] - cs[j]);
   // tbox::pout<<"zeta="<<zeta<<std::endl;
   const double DL = d_DL * exp(-d_Q0 / (gas_constant * temp));
   zeta /= DL;

   return 1. / (d_alpha + d_beta * zeta);
}

template class KimMobilityStrategyFiniteMobAntiTrap<
    CALPHADFreeEnergyFunctionsBinary>;
template class KimMobilityStrategyFiniteMobAntiTrap<
    CALPHADFreeEnergyFunctionsTernary>;
template class KimMobilityStrategyFiniteMobAntiTrap<
    KKSFreeEnergyFunctionDiluteBinary>;
#ifdef HAVE_THERMO4PFM
template class KimMobilityStrategyFiniteMobAntiTrap<
    CALPHADFreeEnergyFunctionsBinary2Ph1Sl>;
#endif
