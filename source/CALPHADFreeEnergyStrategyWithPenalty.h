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
#ifndef included_CALPHADFreeEnergyStrategyWithPenalty
#define included_CALPHADFreeEnergyStrategyWithPenalty

#if 0

#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "MolarVolumeStrategy.h"

using namespace ampe_thermo;


class CALPHADFreeEnergyStrategyWithPenalty
    : public CALPHADFreeEnergyStrategyBinary<CALPHADFreeEnergyFunctionsBinary>
{
 public:
   CALPHADFreeEnergyStrategyWithPenalty(
       std::shared_ptr<tbox::Database> input_db,
       std::shared_ptr<tbox::Database> newton_db,
       const EnergyInterpolationType energy_interp_func_type,
       const ConcInterpolationType conc_interp_func_type,
       MolarVolumeStrategy* mvstrategy, const int conc_l_id,
       const int conc_a_id, const int conc_b_id, const int ncompositions,
       const bool with_third_phase);

   ~CALPHADFreeEnergyStrategyWithPenalty(){};

   virtual void setup(std::shared_ptr<tbox::Database> input_db,
                      std::shared_ptr<tbox::Database> newton_db);

   virtual bool computeCeqT(const double temperature, const PhaseIndex pi0,
                            const PhaseIndex pi1, double* ceq);

   double computeValFreeEnergyLiquid(const double temperature,
                                     const double conc, const bool gp = false);

   double computeValFreeEnergySolidA(const double temperature,
                                     const double conc, const bool gp = false);

   double computeValFreeEnergySolidB(const double temperature,
                                     const double conc, const bool gp = false);

   void computeSecondDerivativeEnergyPhaseL(const double temperature,
                                            const std::vector<double>& c_l,
                                            std::vector<double>& d2fdc2,
                                            const bool use_internal_units)
   {
      std::shared_ptr<CALPHADFreeEnergyFunctionsBinary> calphad_fenergy =
          std::dynamic_pointer_cast<CALPHADFreeEnergyFunctionsBinary>(
              d_calphad_fenergy);
      assert(calphad_fenergy);

      calphad_fenergy->computeSecondDerivativeFreeEnergy(temperature, &c_l[0],
                                                         PhaseIndex::phaseL,
                                                         d2fdc2);

      double extra_energy =
          calphad_fenergy->compute2ndDerivPenalty(PhaseIndex::phaseL, c_l[0]);

      d2fdc2[0] += extra_energy;

      if (use_internal_units)
         d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temperature, &c_l[0],
                                                           PhaseIndex::phaseL);

      if (d2fdc2[0] < 0.)
         tbox::pout << "CALPHADFreeEnergyStrategyWithPenalty --- WARNING: fcc="
                    << d2fdc2[0] << " in liquid for cl=" << c_l[0] << "!!!"
                    << std::endl;
      // tbox::pout<<"fcc="<<d2fdc2[0]<<" in liquid!!!"<<endl;
   }

   void computeSecondDerivativeEnergyPhaseA(const double temperature,
                                            const std::vector<double>& c_a,
                                            std::vector<double>& d2fdc2,
                                            const bool use_internal_units)
   {
      std::shared_ptr<CALPHADFreeEnergyFunctionsBinary> calphad_fenergy =
          std::dynamic_pointer_cast<CALPHADFreeEnergyFunctionsBinary>(
              d_calphad_fenergy);

      calphad_fenergy->computeSecondDerivativeFreeEnergy(temperature, &c_a[0],
                                                         PhaseIndex::phaseA,
                                                         d2fdc2);

      double extra_energy =
          calphad_fenergy->compute2ndDerivPenalty(PhaseIndex::phaseA, c_a[0]);

      d2fdc2[0] += extra_energy;

      if (use_internal_units)
         d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temperature, &c_a[0],
                                                           PhaseIndex::phaseA);

      if (d2fdc2[0] < 0.)
         tbox::pout << "CALPHADFreeEnergyStrategyWithPenalty --- WARNING: fcc="
                    << d2fdc2[0] << " in phase A for ca=" << c_a[0] << "!!!"
                    << std::endl;
      // tbox::pout<<"fcc="<<d2fdc2[0]<<" in liquid!!!"<<endl;
   }

   void computeSecondDerivativeEnergyPhaseB(const double temperature,
                                            const std::vector<double>& c_b,
                                            std::vector<double>& d2fdc2,
                                            const bool use_internal_units)
   {
      std::shared_ptr<CALPHADFreeEnergyFunctionsBinary> calphad_fenergy =
          std::dynamic_pointer_cast<CALPHADFreeEnergyFunctionsBinary>(
              d_calphad_fenergy);

      calphad_fenergy->computeSecondDerivativeFreeEnergy(temperature, &c_b[0],
                                                         PhaseIndex::phaseB,
                                                         d2fdc2);

      double extra_energy =
          calphad_fenergy->compute2ndDerivPenalty(PhaseIndex::phaseB, c_b[0]);

      d2fdc2[0] += extra_energy;

      if (use_internal_units)
         d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temperature, &c_b[0],
                                                           PhaseIndex::phaseB);

      if (d2fdc2[0] < 0.)
         tbox::pout << "CALPHADFreeEnergyStrategyWithPenalty --- WARNING: fcc="
                    << d2fdc2[0] << " in phase B for cb=" << c_b[0] << "!!!"
                    << std::endl;
      // tbox::pout<<"fcc="<<d2fdc2[0]<<" in liquid!!!"<<endl;
   }

 private:
   std::vector<std::vector<double> > d_penalty_parameters;
};
#endif

#endif
