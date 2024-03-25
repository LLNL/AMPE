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
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "InterpolationType.h"

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>

using namespace SAMRAI;
#include "Database2JSON.h"
namespace pt = boost::property_tree;
using namespace Thermo4PFM;


int main(int argc, char *argv[])
{
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

      std::cout << "input_filename = " << input_filename << std::endl;

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

      std::shared_ptr<tbox::Database> newton_db;
      if (conc_db->isDatabase("NewtonSolver"))
         newton_db = conc_db->getDatabase("NewtonSolver");

      pt::ptree calphad_pt;
      pt::ptree newton_pt;
      copyDatabase(calphad_db, calphad_pt);
      copyDatabase(newton_db, newton_pt);

      CALPHADFreeEnergyFunctionsBinary cafe(calphad_pt, newton_pt,
                                            energy_interp_func_type,
                                            conc_interp_func_type);

      // initial guesses
      double init_guess[2];
      model_db->getDoubleArray("initial_guess", &init_guess[0], 2);

      double lceq[2] = {init_guess[0], init_guess[1]};

      // compute equilibrium concentrations in each phase
      bool found_ceq = cafe.computeCeqT(temperature, &lceq[0]);
      if (lceq[0] > 1.) found_ceq = false;
      if (lceq[0] < 0.) found_ceq = false;
      if (lceq[1] > 1.) found_ceq = false;
      if (lceq[1] < 0.) found_ceq = false;

      std::cout << "Temperature = " << temperature << std::endl;
      if (found_ceq) {
         std::cout << "Found equilibrium concentrations: " << lceq[0] << " and "
                   << lceq[1] << "..." << std::endl;
         ret = 0;
      } else {
         std::cout << "Equilibrium concentrations not found!" << std::endl;
         ret = 1;
      }

      std::vector<double> d2fdc2(1, 0.);
      cafe.computeSecondDerivativeFreeEnergy(temperature, lceq,
                                             PhaseIndex::phaseL, d2fdc2.data());
      std::cout << "2nd derivative of fL [J/mol]: " << d2fdc2[0] << std::endl;

      cafe.computeSecondDerivativeFreeEnergy(temperature, lceq,
                                             PhaseIndex::phaseA, d2fdc2.data());
      std::cout << "2nd derivative of fS [J/mol]: " << d2fdc2[0] << std::endl;

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
