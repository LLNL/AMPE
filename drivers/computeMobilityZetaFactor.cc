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
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "CALPHADMobility.h"
#include "KKStools.h"

#include "SAMRAI/SAMRAI_config.h"
#include "Database2JSON.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <fstream>

using namespace SAMRAI;
namespace pt = boost::property_tree;
using namespace Thermo4PFM;

// compute mobility parameter according to
// S.G. Kim, Acta mat. 2007
// Run (example):
// computeMobilityCALPHADbinary mobilityAuNi.input 1400.
int main(int argc, char *argv[])
{
   // Initialize MPI, SAMRAI, and enable logging.

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      std::string input_filename = argv[1];
      double temperature = atof(argv[2]);

      //-----------------------------------------------------------------------
      // Create input database and parse all data in input file.

      std::shared_ptr<tbox::MemoryDatabase> input_db(
          new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                       input_db);

      std::cout << "Compute PFM Kim's mobility for binary alloy" << std::endl;
      std::cout << "input_filename = " << input_filename << std::endl;

      std::shared_ptr<tbox::Database> model_db =
          input_db->getDatabase("ModelParameters");

      EnergyInterpolationType energy_interp_func_type =
          EnergyInterpolationType::PBG;
      ConcInterpolationType conc_interp_func_type =
          ConcInterpolationType::LINEAR;

      // read two phases to pick from
      std::string phase0 = model_db->getString("phaseL");
      std::string phase1 = model_db->getString("phaseS");

      double DL = model_db->getDouble("D_liquid");
      double Q0 = model_db->getDouble("Q0_liquid");

      const PhaseIndex pi0 = PhaseIndex::phaseL;
      double lceq[2];
      std::vector<double> d2fdc2(1);

      if (model_db->keyExists("calphad_filename")) {
         std::string calphad_filename = model_db->getString("calphad_filename");

         pt::ptree calphad_pt;
         if (calphad_filename.compare(calphad_filename.size() - 4, 4, "json") ==
             0) {
            boost::property_tree::read_json(calphad_filename, calphad_pt);
         } else {
            std::shared_ptr<tbox::MemoryDatabase> calphad_db(
                new tbox::MemoryDatabase("calphad_db"));
            tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                             calphad_db);
            copyDatabase(calphad_db, calphad_pt);
         }

         pt::ptree newton_pt;
         CALPHADFreeEnergyFunctionsBinary2Ph1Sl cafe(calphad_pt, newton_pt,
                                                     energy_interp_func_type,
                                                     conc_interp_func_type,
                                                     phase0, phase1);

         // initial guesses
         double init_guess[2];
         model_db->getDoubleArray("initial_guess", &init_guess[0], 2);

         lceq[0] = init_guess[0];
         lceq[1] = init_guess[1];

         bool found_ceq = cafe.computeCeqT(temperature, &lceq[0], 25);
         if (lceq[0] > 1.) found_ceq = false;
         if (lceq[0] < 0.) found_ceq = false;
         if (lceq[1] > 1.) found_ceq = false;
         if (lceq[1] < 0.) found_ceq = false;

         if (found_ceq) {
            std::cout << "At temperature " << temperature
                      << ", found equilibrium concentrations: " << std::endl;
            std::cout << "Liquid: " << lceq[0] << std::endl;
            std::cout << "Solid:  " << lceq[1] << std::endl;

            cafe.computeSecondDerivativeFreeEnergy(temperature, &lceq[0], pi0,
                                                   &d2fdc2[0]);
         } else {
            std::cerr << "ERROR: Equilibrium concentrations not found... "
                      << std::endl;
            std::cerr << "Cannot compute mobility" << std::endl;
         }
      } else {
         pt::ptree conc_pt;
         copyDatabase(model_db, conc_pt);

         KKSFreeEnergyFunctionDiluteBinary fenergy(conc_pt,
                                                   energy_interp_func_type,
                                                   conc_interp_func_type);
         fenergy.computeCeqT(temperature, &lceq[0], 0, true);
         fenergy.computeSecondDerivativeFreeEnergy(temperature, &lceq[0], pi0,
                                                   &d2fdc2[0]);
      }
      std::cout << "d2fdc2 = " << d2fdc2[0] << std::endl;

      const double diffCoeff = DL * exp(-Q0 / (gas_constant * temperature));
      std::cout << "Diffusion coeff = " << diffCoeff << std::endl;
      double zeta = (lceq[0] - lceq[1]) * d2fdc2[0] * (lceq[0] - lceq[1]);

      zeta /= diffCoeff;

      std::cout << "zeta = " << zeta << std::endl;
   }
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return (0);
}
