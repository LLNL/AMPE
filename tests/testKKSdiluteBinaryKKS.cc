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
#include "KKSFreeEnergyFunctionDiluteBinary.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>

using namespace SAMRAI;

#include "Database2JSON.h"
namespace pt = boost::property_tree;
using namespace Thermo4PFM;

int main(int argc, char* argv[])
{
   // Initialize MPI, SAMRAI

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /* This extra code block is used to scope some temporaries that are
    * created, it forces the destruction before the manager is
    * shutdown.
    */
   {
      std::string input_filename;
      input_filename = argv[1];

      //-----------------------------------------------------------------------
      // Create input database and parse all data in input file.

      std::shared_ptr<tbox::MemoryDatabase> input_db(
          new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                       input_db);

      // make from input file name
      std::string run_name =
          input_filename.substr(0, input_filename.rfind("."));

      // Logfile
      std::string log_file_name = run_name + ".log";

      tbox::PIO::logOnlyNodeZero(log_file_name);

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) tbox::plog << " AMPE: git version " << xstr(x) << std::endl;
      LOG(GITVERSION);
      tbox::plog << std::endl;
#endif

      tbox::plog << "input_filename = " << input_filename << std::endl;

      std::shared_ptr<tbox::Database> model_db =
          input_db->getDatabase("ModelParameters");

      EnergyInterpolationType energy_interp_func_type =
          EnergyInterpolationType::PBG;
      ConcInterpolationType conc_interp_func_type =
          ConcInterpolationType::LINEAR;

      std::shared_ptr<tbox::Database> temperature_db =
          model_db->getDatabase("Temperature");
      double temperature = temperature_db->getDouble("temperature");

      std::shared_ptr<tbox::Database> conc_db(
          model_db->getDatabase("ConcentrationModel"));

      pt::ptree conc_pt;
      copyDatabase(conc_db, conc_pt);

      KKSFreeEnergyFunctionDiluteBinary cafe(conc_pt, energy_interp_func_type,
                                             conc_interp_func_type);

      const double tol = 1.e-5;

      const PhaseIndex pi0 = PhaseIndex::phaseL;
      const PhaseIndex pi1 = PhaseIndex::phaseA;

      // initial guesses
      double c_init0 = 0.5;
      double c_init1 = 0.5;

      double sol[2] = {c_init0, c_init1};

      // solve KKS equations
      double conc = model_db->getDouble("concentration");
      double phi = model_db->getDouble("phi");
      cafe.computePhaseConcentrations(temperature, &conc, &phi, &sol[0]);

      tbox::pout << "-------------------------------" << std::endl;
      tbox::pout << "Temperature = " << temperature << std::endl;
      tbox::pout << "Result for c = " << conc << " and phi = " << phi
                 << std::endl;
      tbox::pout << "   cL = " << sol[0] << std::endl;
      tbox::pout << "   cS = " << sol[1] << std::endl;

      // verify concentrations satisfy equal chemical potentials
      tbox::pout << "Verification:" << std::endl;

      double derivL;
      double derivS;
      cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL);
      tbox::pout << "   dfL/dcL = " << derivL << std::endl;

      cafe.computeDerivFreeEnergy(temperature, &sol[1], pi1, &derivS);
      tbox::pout << "   dfS/dcS = " << derivS << std::endl;

      if (fabs(derivS - derivL) < tol) {
         tbox::pout << "TEST PASSED" << std::endl;
      } else {
         tbox::pout << "TEST FAILED\n!";
         tbox::pout << "Difference between derivatives: " << derivS - derivL
                    << std::endl;
         return 1;
      }

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
