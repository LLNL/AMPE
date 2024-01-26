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
//
#ifndef included_CompositionStrategyMobilities
#define included_CompositionStrategyMobilities

#include "CALPHADMobility.h"
#include "FreeEnergyStrategy.h"

#include "SAMRAI/tbox/Database.h"

#include <string>

class CompositionStrategyMobilities
{
 public:
   CompositionStrategyMobilities(
       std::shared_ptr<tbox::Database> input_db, const bool,
       const unsigned short ncompositions,
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy);

   virtual ~CompositionStrategyMobilities(){};

   void printDiagnostics(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                         const int temperature_scratch_id);
   void printDiagnostics(const double Tmin, const double Tmax);

   void computeDiffusionMobilityPhaseL(const std::vector<double>& c,
                                       const double temp,
                                       std::vector<double>& mobility);
   void computeDiffusionMobilityPhaseA(const std::vector<double>& c,
                                       const double temp,
                                       std::vector<double>& mobility);
   void computeDiffusionMobilityPhaseB(const std::vector<double>& c,
                                       const double temp,
                                       std::vector<double>& mobility);

   double computeDiffusionMobilityBinaryPhaseL(const double, const double);
   double computeDiffusionMobilityBinaryPhaseA(const double, const double);
   double computeDiffusionMobilityBinaryPhaseB(const double, const double);

 protected:
   // free energy needed to compute diffusion in each phase
   std::shared_ptr<FreeEnergyStrategy> d_free_energy_strategy;

   void printMobilitiesVsComposition(const double temperature,
                                     std::ostream& os);
   void printDiffusionVsComposition(const double temperature, std::ostream& os);

 private:
   unsigned short d_ncompositions;

   bool d_with_third_phase;

   std::vector<CALPHADMobility> d_calphad_mobilities_phaseL;
   std::vector<CALPHADMobility> d_calphad_mobilities_phaseA;
   std::vector<CALPHADMobility> d_calphad_mobilities_phaseB;

   void computeDiffusionMobilityTernaryPhaseL(const double c0, const double c1,
                                              const double temp,
                                              std::vector<double>& mobility);

   void computeDiffusionMobilityTernaryPhaseA(const double c0, const double c1,
                                              const double temp,
                                              std::vector<double>& mobility);

   void computeDiffusionMobilityTernaryPhaseB(const double c0, const double c1,
                                              const double temp,
                                              std::vector<double>& mobility);

   void computeDiffusionMobilityPhase(const char phase,
                                      const std::vector<double>& c,
                                      const double temp,
                                      std::vector<double>& mobility)
   {
      switch (phase) {
         case 'l': computeDiffusionMobilityPhaseL(c, temp, mobility); break;

         case 'a': computeDiffusionMobilityPhaseA(c, temp, mobility); break;

         case 'b': computeDiffusionMobilityPhaseB(c, temp, mobility); break;

         default:
            tbox::pout << "CompositionStrategyMobilities::"
                          "computeDiffusionMobilityPhase(), Error: phase="
                       << phase << "!!!" << std::endl;
            tbox::SAMRAI_MPI::abort();
      }
   }
};

#endif
