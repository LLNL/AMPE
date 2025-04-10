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
#include "ParabolicFreeEnergyFunctionsBinary.h"
#include "ParabolicTools.h"

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
#include <iomanip>

using namespace SAMRAI;

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

      Thermo4PFM::EnergyInterpolationType energy_interp_func_type =
          Thermo4PFM::EnergyInterpolationType::PBG;
      Thermo4PFM::ConcInterpolationType conc_interp_func_type =
          Thermo4PFM::ConcInterpolationType::PBG;

      std::shared_ptr<tbox::Database> temperature_db =
          input_db->getDatabase("Temperature");
      double temperature_low = temperature_db->getDouble("low");
      double temperature_high = temperature_db->getDouble("high");


      double coeffL[3][2];
      readParabolicData(input_db, "Liquid", coeffL);

      double coeffA[3][2];
      readParabolicData(input_db, "PhaseA", coeffA);

      double Tref = input_db->getDouble("Tref");

      Thermo4PFM::ParabolicFreeEnergyFunctionsBinary fe(Tref, coeffL, coeffA,
                                                        energy_interp_func_type,
                                                        conc_interp_func_type);

      // initial guesses
      double init_guess[2];
      input_db->getDoubleArray("initial_guess", &init_guess[0], 2);

      double lceq[2] = {init_guess[0], init_guess[1]};

      std::map<double, double> cseq;
      std::map<double, double> cleq;

      double dT = (temperature_high - temperature_low) / 50;

      input_db->printClassData(tbox::plog);

      const double tol = 1.e-4;

      // loop over temperature range
      for (int iT = 0; iT < 50; iT++) {

         double temperature = temperature_low + iT * dT;

         // compute equilibrium concentrations
         bool found_ceq = fe.computeCeqT(temperature, &lceq[0], 50, true);
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
         std::ofstream os("CvsT.csv");
         os << std::setprecision(9);
         os << "T, ceqL, ceqS\n";
         {
            std::map<double, double>::iterator itl = cleq.begin();
            std::map<double, double>::iterator its = cseq.begin();
            while (itl != cleq.end()) {
               os << itl->first << ", " << itl->second << ", " << its->second
                  << std::endl;
               ++itl;
               ++its;
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
