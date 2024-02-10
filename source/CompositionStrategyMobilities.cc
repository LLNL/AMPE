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
#include "CompositionStrategyMobilities.h"
#include "CALPHADFunctions.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/tbox/InputManager.h"

#include <cassert>

const double m2toum2 = 1.e12;

CompositionStrategyMobilities::CompositionStrategyMobilities(
    std::shared_ptr<tbox::Database> input_db, const bool with_third_phase,
    const unsigned short ncompositions,
    std::shared_ptr<FreeEnergyStrategy> free_energy_strategy)
    : d_free_energy_strategy(free_energy_strategy)
{
   assert(input_db);
   assert(d_free_energy_strategy);

   d_ncompositions = ncompositions;

   d_with_third_phase = with_third_phase;

   std::string calphad_filename = input_db->getString("filename");

   std::shared_ptr<tbox::MemoryDatabase> calphad_db(
       new tbox::MemoryDatabase("calphad_db"));
   tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                    calphad_db);

   std::shared_ptr<tbox::Database> mobility_db =
       calphad_db->getDatabase("MobilityParameters");
   std::shared_ptr<tbox::Database> species0_db =
       mobility_db->getDatabase("Species0");
   std::shared_ptr<tbox::Database> species1_db =
       mobility_db->getDatabase("Species1");
   std::shared_ptr<tbox::Database> species2_db;
   if (d_ncompositions > 1) {
      species2_db = mobility_db->getDatabase("Species2");
   }

   CALPHADMobility calphad_mobility0_phaseL("MobilitySpecies0");
   calphad_mobility0_phaseL.initialize(species0_db->getDatabase("PhaseL"));

   CALPHADMobility calphad_mobility1_phaseL("MobilitySpecies1");
   calphad_mobility1_phaseL.initialize(species1_db->getDatabase("PhaseL"));

   d_calphad_mobilities_phaseL.push_back(calphad_mobility0_phaseL);
   d_calphad_mobilities_phaseL.push_back(calphad_mobility1_phaseL);

   CALPHADMobility calphad_mobility0_phaseA("MobilitySpecies0");
   calphad_mobility0_phaseA.initialize(species0_db->getDatabase("PhaseA"));

   CALPHADMobility calphad_mobility1_phaseA("MobilitySpecies1");
   calphad_mobility1_phaseA.initialize(species1_db->getDatabase("PhaseA"));

   d_calphad_mobilities_phaseA.push_back(calphad_mobility0_phaseA);
   d_calphad_mobilities_phaseA.push_back(calphad_mobility1_phaseA);

   if (d_ncompositions > 1) {
      CALPHADMobility calphad_mobility2_phaseL("MobilitySpecies2");
      calphad_mobility2_phaseL.initialize(species2_db->getDatabase("PhaseL"));

      CALPHADMobility calphad_mobility2_phaseA("MobilitySpecies2");
      calphad_mobility2_phaseA.initialize(species2_db->getDatabase("PhaseA"));

      d_calphad_mobilities_phaseL.push_back(calphad_mobility2_phaseL);
      d_calphad_mobilities_phaseA.push_back(calphad_mobility2_phaseA);
   }

   CALPHADMobility calphad_mobility0_phaseB("MobilitySpecies0");
   CALPHADMobility calphad_mobility1_phaseB("MobilitySpecies1");

   if (d_with_third_phase) {
      calphad_mobility0_phaseB.initialize(species0_db->getDatabase("PhaseB"));
      calphad_mobility1_phaseB.initialize(species1_db->getDatabase("PhaseB"));

      d_calphad_mobilities_phaseB.push_back(calphad_mobility0_phaseB);
      d_calphad_mobilities_phaseB.push_back(calphad_mobility1_phaseB);
   }
}

//-----------------------------------------------------------------------

void CompositionStrategyMobilities::printDiagnostics(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_scratch_id)
{
   tbox::plog << "CompositionRHSStrategy::printDiagnostics()" << std::endl;

   math::HierarchyCellDataOpsReal<double> cell_ops(hierarchy, 0, 0);
   const double Tmax = cell_ops.max(temperature_scratch_id);
   const double Tmin = cell_ops.min(temperature_scratch_id);

   printDiagnostics(Tmin, Tmax);
}

//-----------------------------------------------------------------------

void CompositionStrategyMobilities::printDiagnostics(const double Tmin,
                                                     const double Tmax)
{
   tbox::plog << "CompositionRHSStrategy::printDiagnostics()" << std::endl;

   assert(Tmin > 0.);
   assert(Tmax > 0.);
   assert(Tmin < 100000.);
   assert(Tmax < 100000.);

   int nT = 0;
   double dT = 0.;
   if (Tmax - Tmin > 0.) {
      nT = 10;
      dT = (Tmax - Tmin) / nT;
   }

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getRank() == 0)
      for (int iT = 0; iT <= nT; iT++) {
         std::string filename("DvsC_T");
         double temperature = Tmin + iT * dT;
         std::stringstream ss(std::stringstream::in | std::stringstream::out);
         ss << temperature;
         std::string temps = ss.str();
         filename += temps;
         filename += "K.dat";

         std::ofstream tfile(&filename[0], std::ios::out);
         tfile << "#Diffusion vs. composition" << std::endl;
         printDiffusionVsComposition(temperature, tfile);
         // printDiffusionVsComposition( temperature, tbox::plog );
      }
}

//=======================================================================

double CompositionStrategyMobilities::computeDiffusionMobilityBinaryPhaseL(
    const double c0, const double temp)
{
   assert(d_ncompositions == 1);

   return computeDiffusionMobilityBinaryPhase(c0, temp,
                                              d_calphad_mobilities_phaseL);
}

double CompositionStrategyMobilities::computeDiffusionMobilityBinaryPhaseA(
    const double c0, const double temp)
{
   return computeDiffusionMobilityBinaryPhase(c0, temp,
                                              d_calphad_mobilities_phaseA);
}

double CompositionStrategyMobilities::computeDiffusionMobilityBinaryPhaseB(
    const double c0, const double temp)
{
   return computeDiffusionMobilityBinaryPhase(c0, temp,
                                              d_calphad_mobilities_phaseB);
}

//=======================================================================

void CompositionStrategyMobilities::computeDiffusionMobilityTernaryPhaseL(
    const double c0, const double c1, const double temp,
    std::vector<double>& mobility)
{
   assert(d_ncompositions == 2);

   computeDiffusionMobilityTernaryPhase(c0, c1, temp,
                                        d_calphad_mobilities_phaseL, mobility);
}

void CompositionStrategyMobilities::computeDiffusionMobilityTernaryPhaseA(
    const double c0, const double c1, const double temp,
    std::vector<double>& mobility)
{
   computeDiffusionMobilityTernaryPhase(c0, c1, temp,
                                        d_calphad_mobilities_phaseA, mobility);
}

void CompositionStrategyMobilities::computeDiffusionMobilityTernaryPhaseB(
    const double c0, const double c1, const double temp,
    std::vector<double>& mobility)
{
   computeDiffusionMobilityTernaryPhase(c0, c1, temp,
                                        d_calphad_mobilities_phaseB, mobility);
}

//=======================================================================

void CompositionStrategyMobilities::computeDiffusionMobilityPhaseL(
    const std::vector<double>& c, const double temp,
    std::vector<double>& mobility)
{
   switch (c.size()) {
      case 1:
         mobility[0] = computeDiffusionMobilityBinaryPhaseL(c[0], temp);
         break;

      case 2:
         computeDiffusionMobilityTernaryPhaseL(c[0], c[1], temp, mobility);
         break;

      default:
         tbox::pout << "size of c std::vector=" << c.size() << std::endl;
         tbox::pout << "Error: diffusion mobility implemented for up to 3 "
                       "species only!!!"
                    << std::endl;
         tbox::SAMRAI_MPI::abort();
   }
}

void CompositionStrategyMobilities::computeDiffusionMobilityPhaseA(
    const std::vector<double>& c, const double temp,
    std::vector<double>& mobility)
{
   switch (c.size()) {
      case 1:
         mobility[0] = computeDiffusionMobilityBinaryPhaseA(c[0], temp);
         break;

      case 2:
         computeDiffusionMobilityTernaryPhaseA(c[0], c[1], temp, mobility);
         break;

      default:
         tbox::pout << "Error: diffusion mobility implemented for up to 3 "
                       "species only!!!"
                    << std::endl;
         tbox::SAMRAI_MPI::abort();
   }

   assert(mobility[0] >= 0.);
}

void CompositionStrategyMobilities::computeDiffusionMobilityPhaseB(
    const std::vector<double>& c, const double temp,
    std::vector<double>& mobility)
{
   switch (c.size()) {
      case 1:
         mobility[0] = computeDiffusionMobilityBinaryPhaseB(c[0], temp);
         break;

      case 2:
         computeDiffusionMobilityTernaryPhaseB(c[0], c[1], temp, mobility);
         break;

      default:
         tbox::pout << "Error: diffusion mobility implemented for up to 3 "
                       "species only!!!"
                    << std::endl;
         tbox::SAMRAI_MPI::abort();
   }
}


//=======================================================================

void CompositionStrategyMobilities::printMobilitiesVsComposition(
    const double temperature, std::ostream& os)
{
   assert(temperature > 0.);
   assert(temperature < 100000.);
   assert(d_free_energy_strategy);

   int nc = 20;
   double eps = 1.e-7;
   double dc = (1. - 2 * eps) / nc;
   std::vector<double> amob(1);
   std::vector<double> conc(d_ncompositions);
   const double rt = gas_constant * temperature;

   os << std::fixed;
   os << "#Species 0: log10(atomic mobility*RT[m2/s]) vs. composition in phase "
         "L at T="
      << temperature << std::endl;
   for (int i = 0; i <= nc; i++) {
      conc[0] = eps + i * dc;
      const double m0 =
          d_calphad_mobilities_phaseL[0].getAtomicMobilityBinary(conc[0],
                                                                 temperature);
      assert(m0 > 0.);

      os << conc[0] << "  " << log10(m0 * rt) << std::endl;
   }

   os << "#Species 0: log10(atomic mobility*RT[m2/s]) vs. composition in phase "
         "A at T="
      << temperature << std::endl;
   for (int i = 0; i <= nc; i++) {
      conc[0] = eps + i * dc;

      const double m0 =
          d_calphad_mobilities_phaseA[0].getAtomicMobilityBinary(conc[0],
                                                                 temperature);
      assert(m0 > 0.);

      os << conc[0] << "  " << log10(m0 * rt) << std::endl;
   }

   if (d_with_third_phase) {
      os << "#Species 0: log10(atomic mobility*RT[m2/s]) vs. composition in "
            "phase B at T="
         << temperature << std::endl;
      for (int i = 0; i <= nc; i++) {
         conc[0] = eps + i * dc;

         const double m0 = d_calphad_mobilities_phaseB[0]
                               .getAtomicMobilityBinary(conc[0], temperature);
         assert(m0 > 0.);

         os << conc[0] << "  " << log10(m0 * rt) << std::endl;
      }
   }
}

//=======================================================================

void CompositionStrategyMobilities::printDiffusionVsComposition(
    const double temperature, std::ostream& os)
{
   assert(temperature > 0.);
   assert(temperature < 100000.);
   assert(d_free_energy_strategy);

   int nc = 50;
   double eps = 1.e-7;
   double dc = (1. - 2 * eps) / nc;
   std::vector<double> amob(1);
   std::vector<double> conc(d_ncompositions);
   std::vector<double> d2f(d_ncompositions * d_ncompositions);

   os << "#Interdiffusion[m2/s] vs. composition in phase L at T=" << temperature
      << std::endl;
   for (int i = 0; i <= nc; i++) {
      conc[0] = eps + i * dc;
      d_free_energy_strategy->computeSecondDerivativeEnergyPhaseL(temperature,
                                                                  conc, d2f,
                                                                  false);

      computeDiffusionMobilityPhaseL(conc, temperature, amob);

      assert(amob[0] > 0.);
      os << std::fixed << conc[0] << "  " << std::scientific
         << amob[0] * d2f[0] / m2toum2 << std::endl;
   }

   os << std::endl;
   os << "#Interdiffusion[m2/s]  vs. composition in phase A at T="
      << temperature << std::endl;
   for (int i = 0; i <= nc; i++) {
      conc[0] = eps + i * dc;
      d_free_energy_strategy->computeSecondDerivativeEnergyPhaseA(temperature,
                                                                  conc, d2f,
                                                                  false);

      computeDiffusionMobilityPhaseA(conc, temperature, amob);

      assert(amob[0] >= 0.);
      os << std::fixed << conc[0] << "  " << std::scientific
         << amob[0] * d2f[0] / m2toum2 << std::endl;
   }
}
