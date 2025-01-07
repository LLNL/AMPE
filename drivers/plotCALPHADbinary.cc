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
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>

using namespace SAMRAI;

#include "Database2JSON.h"
namespace pt = boost::property_tree;

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

int main(int argc, char *argv[])
{
   // Initialize MPI, SAMRAI

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      std::string input_filename(argv[1]);

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

      std::shared_ptr<tbox::Database> model_db =
          input_db->getDatabase("ModelParameters");

      Thermo4PFM::EnergyInterpolationType energy_interp_func_type =
          Thermo4PFM::EnergyInterpolationType::PBG;
      Thermo4PFM::ConcInterpolationType conc_interp_func_type =
          Thermo4PFM::ConcInterpolationType::PBG;

      std::shared_ptr<tbox::Database> temperature_db =
          model_db->getDatabase("Temperature");
      double temperature = temperature_db->getDouble("temperature");

      std::shared_ptr<tbox::Database> conc_db(
          model_db->getDatabase("ConcentrationModel"));

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
      if (conc_db->isDatabase("NewtonSolver")) {
         std::shared_ptr<tbox::Database> newton_db;
         newton_db = conc_db->getDatabase("NewtonSolver");
         copyDatabase(newton_db, newton_pt);
      }

      Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(calphad_pt, newton_pt,
                                                        energy_interp_func_type,
                                                        conc_interp_func_type);

      std::ofstream os("FvsC.dat", std::ios::out);
      cafe.printEnergyVsComposition(temperature, os);

      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
