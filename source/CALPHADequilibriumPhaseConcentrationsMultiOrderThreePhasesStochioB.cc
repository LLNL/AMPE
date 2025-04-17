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
#include "Database2JSON.h"

#include <boost/property_tree/json_parser.hpp>
namespace pt = boost::property_tree;

#include "CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

// Thermo4PFM
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB ::
    CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB(
        const short norderp_A, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const QuatModelParameters& model_parameters,
        std::shared_ptr<tbox::Database> conc_db,
        std::shared_ptr<tbox::Database> newton_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases(
          norderp_A, conc_l_id, conc_a_id, conc_b_id, model_parameters,
          conc_db),
      d_model_parameters(model_parameters)
{
   tbox::plog << "CALPHADequilPhaseConcMultiOrderThreePhasesStochioB..."
              << std::endl;
   std::shared_ptr<tbox::Database> conc_calphad_db =
       conc_db->getDatabase("Calphad");
   std::string calphad_filename = conc_calphad_db->getString("filename");

   std::shared_ptr<tbox::MemoryDatabase> calphad_db;
   boost::property_tree::ptree calphad_pt;

   if (calphad_filename.compare(calphad_filename.size() - 4, 4, "json") == 0) {
      boost::property_tree::read_json(calphad_filename, calphad_pt);
   } else {
      calphad_db.reset(new tbox::MemoryDatabase("calphad_db"));
      tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                       calphad_db);
      copyDatabase(calphad_db, calphad_pt);
   }

   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);

   d_calphad_fenergy = std::unique_ptr<
       Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>(
       new Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB(
           model_parameters.getStochioB(), calphad_pt, newton_pt,
           Thermo4PFM::EnergyInterpolationType::LINEAR,
           Thermo4PFM::ConcInterpolationType::LINEAR));
}

int CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB ::
    computePhaseConcentrations(const double temp, double* c, double* hphi,
                               double* x)
{
   assert(!std::isnan(x[0]));
   assert(!std::isnan(x[1]));

   const double epsilon1 = 1.e-4;
   const double epsilon2 = 2.e-4;

   // initialize to NaN to trigger error if used when not set
   double xkks[2] = {tbox::IEEE::getSignalingNaN(),
                     tbox::IEEE::getSignalingNaN()};
   double xeq[2] = {tbox::IEEE::getSignalingNaN(),
                    tbox::IEEE::getSignalingNaN()};

   if (hphi[2] < 1. - epsilon1) {
      // solve KKS problem
      xkks[0] = x[0];
      xkks[1] = x[1];
      int status =
          d_calphad_fenergy->computePhaseConcentrations(temp, c, hphi, xkks);
#ifndef GPU_OFFLOAD
      if (status < 0) {
         std::cerr << "computePhaseConcentrations failed for T=" << temp
                   << ", hphi=";
         for (short i = 0; i < 3; i++)
            std::cerr << hphi[i] << ", ";
         std::cerr << ", c=" << c[0] << std::endl;
         std::cerr << "x = " << x[0] << ", " << x[1] << std::endl;
         std::cerr << "xkks = " << xkks[0] << ", " << xkks[1] << std::endl;
         std::cerr << "conc[0] - hphi2 * cB_ = "
                   << c[0] - hphi[2] * d_model_parameters.getStochioB()
                   << std::endl;
         std::cerr << "cB = " << d_model_parameters.getStochioB() << std::endl;
         const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
         MPI_Abort(mpi.getCommunicator(), EXIT_FAILURE);
      }
#endif
   }

   if (hphi[2] > 1. - epsilon2) {
      xeq[0] = d_model_parameters.ceq_liquid(temp);
      xeq[1] = d_model_parameters.ceq_solidA(temp);
      // std::cout << "xeq[0] = " << xeq[0] << std::endl;
      // std::cout << "xeq[1] = " << xeq[1] << std::endl;
      assert(!std::isnan(xeq[0]));
      assert(!std::isnan(xeq[1]));
   }

   if (hphi[2] < 1. - epsilon2) {
      x[0] = xkks[0];
      x[1] = xkks[1];
   } else if (hphi[2] > 1. - epsilon1) {
      x[0] = xeq[0];
      x[1] = xeq[1];
   } else {
      // map 1.-epsilon2 < hphi < 1-epsilon1 to (0,1)
      double h = (hphi[2] - (1. - epsilon2)) / (epsilon2 - epsilon1);
      // mix equilibrium compositions and KKS solution
      double f =
          Thermo4PFM::interp_func(Thermo4PFM::EnergyInterpolationType::PBG, h);
      x[0] = f * xeq[0] + (1. - f) * xkks[0];
      x[1] = f * xeq[1] + (1. - f) * xkks[1];
   }
   return 0;
}
