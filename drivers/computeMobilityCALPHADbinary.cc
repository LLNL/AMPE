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

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"

#include "SAMRAI/SAMRAI_config.h"
#include "Database2JSON.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"
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
int main(int argc, char *argv[])
{
   // Initialize MPI, SAMRAI, and enable logging.

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      std::string input_filename = argv[1];

      //-----------------------------------------------------------------------
      // Create input database and parse all data in input file.

      std::shared_ptr<tbox::MemoryDatabase> input_db(
          new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                       input_db);

      std::cout << "Compute PFM mobility for binary alloy" << std::endl;
      std::cout << "input_filename = " << input_filename << std::endl;

      std::shared_ptr<tbox::Database> model_db =
          input_db->getDatabase("ModelParameters");

      double phase_well_scale = model_db->getDouble("phi_well_scale");
      double epsilon = model_db->getDouble("epsilon_phi");

      EnergyInterpolationType energy_interp_func_type =
          EnergyInterpolationType::PBG;
      ConcInterpolationType conc_interp_func_type =
          ConcInterpolationType::LINEAR;

      std::string phase0 = model_db->getString("phaseL");
      std::string phase1 = model_db->getString("phaseS");

      double temperature = model_db->getDouble("temperature");

      std::shared_ptr<tbox::Database> conc_db(
          model_db->getDatabase("ConcentrationModel"));

      double mv = conc_db->getDouble("molar_volume");

      double DL = conc_db->getDouble("D_liquid");

      std::shared_ptr<tbox::Database> dcalphad_db =
          conc_db->getDatabase("Calphad");
      std::string calphad_filename = dcalphad_db->getString("filename");

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
                                                  conc_interp_func_type, phase0,
                                                  phase1);

      // initial guesses
      double init_guess[2];
      model_db->getDoubleArray("initial_guess", &init_guess[0], 2);

      double nominalc;
      model_db->getDoubleArray("concentration", &nominalc, 1);
      double lceq[5] = {init_guess[0], init_guess[1]};

      bool found_ceq = cafe.computeCeqT(temperature, &lceq[0], 25);
      if (lceq[0] > 1.) found_ceq = false;
      if (lceq[0] < 0.) found_ceq = false;
      if (lceq[1] > 1.) found_ceq = false;
      if (lceq[1] < 0.) found_ceq = false;

      if (found_ceq) {
         std::cout << "For nominal composition " << nominalc
                   << ", found equilibrium concentrations: " << std::endl;
         std::cout << "Liquid: " << lceq[0] << std::endl;
         std::cout << "Solid:  " << lceq[1] << std::endl;

         std::cout << "Interfacial energy: "
                   << epsilon * sqrt(16. * phase_well_scale) / (3. * sqrt(2.))
                   << " (J/m^2)" << std::endl;
         std::cout << "Delta: " << epsilon / sqrt(32. * phase_well_scale)
                   << " (um)" << std::endl;

         const PhaseIndex pi0 = PhaseIndex::phaseL;
         std::vector<double> d2fdc2(1);
         cafe.computeSecondDerivativeFreeEnergy(temperature, &lceq[0], pi0,
                                                &d2fdc2[0]);
         std::cout << "d2fdc2 = " << d2fdc2[0] << std::endl;

         double zeta = (lceq[0] - lceq[1]) * d2fdc2[0] * (lceq[0] - lceq[1]);
         zeta /= DL;
         zeta *= (1.e-6 / mv);  // convert from J/mol to pJ/um^3

         std::cout << "zeta = " << zeta << std::endl;

         double xi = epsilon / sqrt(32. * phase_well_scale);
         double a2 = 47. / 60.;
         double mobility = 1. / (3. * (2. * xi * xi) * a2 * zeta);

         std::cout << "mobility = " << mobility << std::endl;

      } else {
         std::cout << "ERROR: Equilibrium concentrations not found... "
                   << std::endl;
         std::cout << "Cannot compute mobility" << std::endl;
      }

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return (0);
}
#endif
