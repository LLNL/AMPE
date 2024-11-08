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

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <map>
#include <iostream>
#include <fstream>

using namespace SAMRAI;
#include "Database2JSON.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

int main(int argc, char *argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      std::string input_filename(argv[1]);

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

      Thermo4PFM::EnergyInterpolationType energy_interp_func_type =
          Thermo4PFM::EnergyInterpolationType::PBG;
      Thermo4PFM::ConcInterpolationType conc_interp_func_type =
          Thermo4PFM::ConcInterpolationType::PBG;

      std::shared_ptr<tbox::Database> temperature_db =
          model_db->getDatabase("Temperature");
      double temperature_low = temperature_db->getDouble("low");
      double temperature_high = temperature_db->getDouble("high");

      std::shared_ptr<tbox::Database> conc_db(
          model_db->getDatabase("ConcentrationModel"));
      std::string conc_avg_func_type =
          conc_db->getStringWithDefault("avg_func_type", "a");

      std::shared_ptr<tbox::Database> dcalphad_db =
          conc_db->getDatabase("Calphad");
      std::string calphad_filename = dcalphad_db->getString("filename");

      boost::property_tree::ptree calphad_pt;
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

      std::shared_ptr<tbox::Database> newton_db;
      if (conc_db->isDatabase("NewtonSolver"))
         newton_db = conc_db->getDatabase("NewtonSolver");

      pt::ptree newton_pt;
      copyDatabase(newton_db, newton_pt);

      Thermo4PFM::CALPHADFreeEnergyFunctionsBinary cafe(calphad_pt, newton_pt,
                                                        energy_interp_func_type,
                                                        conc_interp_func_type);


      // initial guesses
      double init_guess[2];
      model_db->getDoubleArray("initial_guess", &init_guess[0], 2);

      double lceq[2] = {init_guess[0], init_guess[1]};

      std::map<double, double> cseq;
      std::map<double, double> cleq;

      double dT = (temperature_high - temperature_low) / 50;

      model_db->printClassData(tbox::plog);

      const double tol = 1.e-4;

      // loop over temperature range
      for (int iT = 0; iT < 50; iT++) {

         double temperature = temperature_low + iT * dT;

         // compute equilibrium concentrations
         bool found_ceq = cafe.computeCeqT(temperature, &lceq[0]);
         if (lceq[0] > 1. + tol) found_ceq = false;
         if (lceq[0] < 0. - tol) found_ceq = false;
         if (lceq[1] > 1. + tol) found_ceq = false;
         if (lceq[1] < 0. - tol) found_ceq = false;

         if (found_ceq) {
            // tbox::pout<<"Found equilibrium concentrations: "
            //          <<lceq[0]<<" and "<<lceq[1]<<"..."<<endl;
            cleq.insert(std::pair<double, double>(temperature, lceq[0]));
            cseq.insert(std::pair<double, double>(temperature, lceq[1]));

         } else {
            tbox::pout << "Temperature = " << temperature << std::endl;
            tbox::pout << "ERROR: Equilibrium concentrations not found... "
                       << std::endl;
            return 1;
         }
      }

      {
         std::ofstream os("CvsTliquid.csv");
         os << "T, ceq\n";
         {
            std::map<double, double>::iterator it = cleq.begin();
            while (it != cleq.end()) {
               os << it->first << ", " << it->second << std::endl;
               ++it;
            }
         }
      }

      {
         std::ofstream os("CvsTsolid.csv");
         os << "T, ceq\n";
         {
            std::map<double, double>::iterator it = cseq.begin();
            while (it != cseq.end()) {
               os << it->first << ", " << it->second << std::endl;
               ++it;
            }
         }
      }
      input_db.reset();
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
