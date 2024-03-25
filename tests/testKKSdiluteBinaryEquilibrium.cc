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

      tbox::pout << "-------------------------------" << std::endl;
      tbox::pout << "Temperature = " << temperature << std::endl;

      // compute equilibrium compositions
      double ceq[2];
      cafe.computeCeqT(temperature, &ceq[0]);
      tbox::pout << "   ceL = " << ceq[0] << std::endl;
      tbox::pout << "   ceS = " << ceq[1] << std::endl;

      // test if chemical potentials are equal for equilibrium compositions
      double derivL;
      double derivS;
      cafe.computeDerivFreeEnergy(temperature, &ceq[0], pi0, &derivL);
      cafe.computeDerivFreeEnergy(temperature, &ceq[1], pi1, &derivS);
      tbox::pout << "   dfL/dcL = " << derivL << std::endl;
      tbox::pout << "   dfS/dcS = " << derivS << std::endl;
      if (fabs(derivS - derivL) < tol) {
         tbox::pout << "TEST PASSED" << std::endl;
      } else {
         tbox::pout << "TEST FAILED\n!";
         tbox::pout << "Difference between derivatives: " << derivS - derivL
                    << std::endl;
         return 1;
      }

      // test if driving force is 0 for equilibrium compositions
      double fl = cafe.computeFreeEnergy(temperature, &ceq[0], pi0, false);
      double fa = cafe.computeFreeEnergy(temperature, &ceq[1], pi1, false);
      tbox::pout << "   fL = " << fl << std::endl;
      tbox::pout << "   fS = " << fa << std::endl;
      double diff = fa - fl - derivS * (ceq[1] - ceq[0]);
      if (fabs(diff) < tol) {
         tbox::pout << "TEST PASSED" << std::endl;
      } else {
         tbox::pout << "TEST FAILED\n!";
         tbox::pout << "Driving force not zero: " << diff << std::endl;
         return 1;
      }

      // compute second derivatives for info only
      std::vector<double> d2fdc2(1);
      cafe.computeSecondDerivativeFreeEnergy(temperature, &ceq[0], pi0,
                                             d2fdc2.data());
      tbox::pout << "-------------------------------" << std::endl;
      tbox::pout << "Second derivatives" << std::endl;
      tbox::pout << "At ceL: " << d2fdc2[0] << std::endl;
      cafe.computeSecondDerivativeFreeEnergy(temperature, &ceq[1], pi1,
                                             d2fdc2.data());
      tbox::pout << "At ceS: " << d2fdc2[0] << std::endl;

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
