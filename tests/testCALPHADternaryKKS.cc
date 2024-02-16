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

int main(int argc, char* argv[])
{
   // Initialize MPI, SAMRAI

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   /* This extra code block is used to scope some temporaries that are
    * created, it forces the destruction before the manager is
    * shutdown.
    */
   {

      //-----------------------------------------------------------------------
      /*
       * Process command line arguments and dump to log file.
       *
       *    executable <input file name>
       */

      std::string input_filename;
      input_filename = argv[1];

      //-----------------------------------------------------------------------
      // Create input database and parse all data in input file.

      std::shared_ptr<tbox::MemoryDatabase> input_db(
          new tbox::MemoryDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename,
                                                       input_db);

      //-----------------------------------------------------------------------
      // Read key input settings

      // make from input file name
      std::string run_name =
          input_filename.substr(0, input_filename.rfind("."));

      // Logfile
      std::string log_file_name = run_name + ".log";

      tbox::PIO::logOnlyNodeZero(log_file_name);

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) tbox::pout << " AMPE: git version " << xstr(x) << std::endl;
      LOG(GITVERSION);
      tbox::pout << std::endl;
#endif

      tbox::pout << "Run with " << mpi.getSize() << " MPI tasks" << std::endl;
      tbox::pout << "input_filename = " << input_filename << std::endl;

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
      std::string conc_avg_func_type =
          conc_db->getStringWithDefault("avg_func_type", "a");

      std::shared_ptr<tbox::Database> dcalphad_db =
          conc_db->getDatabase("Calphad");
      std::string calphad_filename = dcalphad_db->getString("filename");
      std::shared_ptr<tbox::MemoryDatabase> calphad_db(
          new tbox::MemoryDatabase("calphad_db"));
      tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                       calphad_db);

      std::shared_ptr<tbox::Database> newton_db;
      if (conc_db->isDatabase("NewtonSolver"))
         newton_db = conc_db->getDatabase("NewtonSolver");

      pt::ptree calphad_pt;
      pt::ptree newton_pt;
      copyDatabase(calphad_db, calphad_pt);
      copyDatabase(newton_db, newton_pt);

      CALPHADFreeEnergyFunctionsTernary cafe(calphad_pt, newton_pt,
                                             energy_interp_func_type,
                                             conc_interp_func_type);

      // initial guesses
      double sol[4];
      model_db->getDoubleArray("initial_guess", &sol[0], 4);

      // compute concentrations satisfying KKS equations
      double conc[2];
      model_db->getDoubleArray("concentration", &conc[0], 2);
      double phi = model_db->getDouble("phi");
      cafe.computePhaseConcentrations(temperature, &conc[0], &phi, &sol[0]);

      tbox::pout << "-------------------------------" << std::endl;
      tbox::pout << "Temperature: " << temperature << std::endl;
      tbox::pout << "Result for c = " << conc[0] << "," << conc[1]
                 << " and phi = " << phi << std::endl;
      tbox::pout << "   cL = " << sol[0] << ", " << sol[1] << std::endl;
      tbox::pout << "   cS = " << sol[2] << ", " << sol[3] << std::endl;

      const PhaseIndex pi0 = PhaseIndex::phaseL;
      const PhaseIndex pi1 = PhaseIndex::phaseA;

      tbox::pout << "Verification:" << std::endl;

      double derivL[2];
      cafe.computeDerivFreeEnergy(temperature, &sol[0], pi0, &derivL[0]);

      double derivS[2];
      cafe.computeDerivFreeEnergy(temperature, &sol[2], pi1, &derivS[0]);

      tbox::pout << "   dfL/dcL0 = " << derivL[0] << std::endl;
      tbox::pout << "   dfS/dcS0 = " << derivS[0] << std::endl;
      tbox::pout << "   dfL/dcL1 = " << derivL[1] << std::endl;
      tbox::pout << "   dfS/dcS1 = " << derivS[1] << std::endl;

      const double tol = 1.e-5;
      if (fabs(derivS[0] - derivL[0]) < tol &&
          fabs(derivS[0] - derivL[0]) < tol) {
         tbox::pout << "TEST PASSED" << std::endl;
      } else {
         tbox::pout << "TEST FAILED!\n";
         tbox::pout << "Difference between derivatives: "
                    << derivS[0] - derivL[0] << "," << derivS[1] - derivL[1]
                    << std::endl;
      }

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
