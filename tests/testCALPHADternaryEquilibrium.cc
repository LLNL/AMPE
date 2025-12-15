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
#include "CALPHADFreeEnergyFunctionsTernary.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <fstream>

using namespace SAMRAI;

#include "Database2JSON.h"
namespace pt = boost::property_tree;
using namespace Thermo4PFM;

int main(int argc, char *argv[])
{
   // feenableexcept(FE_INVALID | FE_OVERFLOW);

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   int ret = 0;
   {

      std::string input_filename(argv[1]);

      // Create input database and parse all data in input file.
      std::shared_ptr<tbox::MemoryDatabase> input_db(
          new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                       input_db);

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) std::cout << " AMPE: git version " << xstr(x) << std::endl;
      LOG(GITVERSION);
      std::cout << std::endl;
#endif

      std::shared_ptr<tbox::Database> model_db =
          input_db->getDatabase("ModelParameters");

      EnergyInterpolationType energy_interp_func_type =
          EnergyInterpolationType::PBG;
      ConcInterpolationType conc_interp_func_type = ConcInterpolationType::PBG;

      std::shared_ptr<tbox::Database> temperature_db =
          model_db->getDatabase("Temperature");
      double temperature = temperature_db->getDouble("temperature");

      std::shared_ptr<tbox::Database> conc_db(
          model_db->getDatabase("ConcentrationModel"));
      std::shared_ptr<tbox::Database> dcalphad_db =
          conc_db->getDatabase("Calphad");
      std::string calphad_filename = dcalphad_db->getString("filename");
      std::shared_ptr<tbox::MemoryDatabase> calphad_db(
          new tbox::MemoryDatabase("calphad_db"));
      tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                       calphad_db);

      int maxits = 20;
      std::shared_ptr<tbox::Database> newton_db;
      if (conc_db->isDatabase("NewtonSolver")) {
         newton_db = conc_db->getDatabase("NewtonSolver");
         maxits = newton_db->getIntegerWithDefault("max_its", 20);
      }

      pt::ptree calphad_pt;
      pt::ptree newton_pt;
      copyDatabase(calphad_db, calphad_pt);
      copyDatabase(newton_db, newton_pt);

      CALPHADFreeEnergyFunctionsTernary cafe(calphad_pt, newton_pt,
                                             energy_interp_func_type,
                                             conc_interp_func_type);


      // initial guesses
      double init_guess[5];
      model_db->getDoubleArray("initial_guess", &init_guess[0], 5);

      double nominalc[2];
      model_db->getDoubleArray("concentration", &nominalc[0], 2);
      double lceq[5] = {init_guess[0], init_guess[1],  // liquid
                        init_guess[2], init_guess[3],  // solid
                        init_guess[4]};

      bool found_ceq = cafe.computeTieLine(temperature, nominalc[0],
                                           nominalc[1], &lceq[0], maxits);

      if (lceq[0] != lceq[0]) found_ceq = false;
      if (lceq[0] > 1.) found_ceq = false;
      if (lceq[0] < 0.) found_ceq = false;
      if (lceq[1] > 1.) found_ceq = false;
      if (lceq[1] < 0.) found_ceq = false;

      std::cout << "Temperature = " << temperature << std::endl;
      if (found_ceq) {
         std::cout << "For nominal composition " << nominalc[0] << ","
                   << nominalc[1]
                   << ", found equilibrium concentrations: " << std::endl;
         std::cout << "Liquid: " << lceq[0] << "," << lceq[1] << std::endl;
         std::cout << "Solid:  " << lceq[2] << "," << lceq[3] << std::endl;
         std::cout << "Solid fraction: " << lceq[4] << std::endl;
      } else {
         std::cout << "TEST FAILED: Equilibrium concentrations not found!"
                   << std::endl;
         ret = 1;
      }

      double expected_cl[2];
      double expected_cs[2];
      std::shared_ptr<tbox::Database> result_db =
          input_db->getDatabase("ExpectedResults");
      result_db->getDoubleArray("cl", &expected_cl[0], 2);
      result_db->getDoubleArray("cs", &expected_cs[0], 2);
      double expected_fs = result_db->getDouble("fs");

      // test values
      const double tol = 1.e-6;
      if (fabs(expected_cl[0] - lceq[0]) > tol) {
         std::cout << "TEST FAILED: cleq[0] != " << expected_cl[0] << std::endl;
         ret = 1;
      }
      if (fabs(expected_cl[1] - lceq[1]) > tol) {
         std::cout << "TEST FAILED: cleq[1] != " << expected_cl[1] << std::endl;
         ret = 1;
      }
      if (fabs(expected_cs[0] - lceq[2]) > tol) {
         std::cout << "TEST FAILED: cseq[0] != " << expected_cs[0] << std::endl;
         ret = 1;
      }
      if (fabs(expected_cs[1] - lceq[3]) > tol) {
         std::cout << "TEST FAILED: cseq[1] != " << expected_cs[1] << std::endl;
         ret = 1;
      }
      if (fabs(expected_fs - lceq[4]) > tol) {
         std::cout << "TEST FAILED: fs != " << expected_fs << std::endl;
         ret = 1;
      }

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
