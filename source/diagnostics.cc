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
#include "CALPHADMobility.h"
#include "QuatModelParameters.h"

#include "SAMRAI/tbox/InputManager.h"

#include <string>
#include <memory>

void preRunDiagnosticsMobilityInPhases(const double temperature,
                                       QuatModelParameters& model_parameters,
                                       std::shared_ptr<tbox::Database> db)
{
   assert(db);
   std::string calphad_filename = db->getString("filename");

   // JSON databases do not contain mobility parameters at the moment
   // jlf 08/27/2021
   if (calphad_filename.compare(calphad_filename.size() - 4, 4, "json") == 0)
      return;

   tbox::pout << "preRunDiagnosticsMobilityInPhases, open " << calphad_filename
              << std::endl;

   std::shared_ptr<tbox::MemoryDatabase> calphad_db(
       new tbox::MemoryDatabase("calphad_db"));
   // WARNING: needs to be called from all MPI tasks!!!
   tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                    calphad_db);

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getRank() == 0)
      if (calphad_db->isDatabase("MobilityParameters")) {
         std::shared_ptr<tbox::Database> mobility_db(
             calphad_db->getDatabase("MobilityParameters"));
         std::shared_ptr<tbox::Database> species0_db(
             mobility_db->getDatabase("Species0"));
         std::shared_ptr<tbox::Database> species1_db(
             mobility_db->getDatabase("Species1"));

         CALPHADMobility calphad_mobility0_phaseL("MobilitySpecies0");
         calphad_mobility0_phaseL.initialize(
             species0_db->getDatabase("PhaseL"));

         CALPHADMobility calphad_mobility1_phaseL("MobilitySpecies1");
         calphad_mobility1_phaseL.initialize(
             species1_db->getDatabase("PhaseL"));

         CALPHADMobility calphad_mobility0_phaseA("MobilitySpecies0");
         calphad_mobility0_phaseA.initialize(
             species0_db->getDatabase("PhaseA"));

         CALPHADMobility calphad_mobility1_phaseA("MobilitySpecies1");
         calphad_mobility1_phaseA.initialize(
             species1_db->getDatabase("PhaseA"));

         const double tempmin = temperature * 0.5;
         const double tempmax = temperature * 2.;

         std::ofstream tfile("D.dat", std::ios::out);
         tfile << "#Diffusion in liquid phase for species 0 [m^2/s] vs. "
                  "10000./T"
               << std::endl;
         calphad_mobility0_phaseL.printDiffusionVsTemperature(tempmin, tempmax,
                                                              tfile);

         tfile << std::endl
               << "#Diffusion in liquid phase for species 1 [m^2/s] vs. "
                  "10000./T"
               << std::endl;
         calphad_mobility0_phaseL.printDiffusionVsTemperature(tempmin, tempmax,
                                                              tfile);

         tfile << std::endl
               << "#Diffusion in phase A for species 0 [m^2/s] vs. "
                  "10000./T"
               << std::endl;
         calphad_mobility0_phaseA.printDiffusionVsTemperature(tempmin, tempmax,
                                                              tfile);

         tfile << std::endl
               << "#Diffusion in phase A for species 1 [m^2/s] vs. "
                  "10000./T"
               << std::endl;
         calphad_mobility1_phaseA.printDiffusionVsTemperature(tempmin, tempmax,
                                                              tfile);
      }
}
