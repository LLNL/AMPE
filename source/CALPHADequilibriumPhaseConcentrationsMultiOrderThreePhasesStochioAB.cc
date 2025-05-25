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

#include "CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioAB.h"

CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioAB ::
    CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioAB(
        const short norderp_A, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const QuatModelParameters& model_parameters,
        std::shared_ptr<tbox::Database> conc_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases(
          norderp_A, conc_l_id, conc_a_id, conc_b_id, model_parameters,
          conc_db),
      d_model_parameters(model_parameters)
{
   tbox::plog << "CALPHADequilPhaseConcMultiOrderThreePhasesStochioAB..."
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
}

int CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioAB ::
    computePhaseConcentrations(const double temp, double* c, double* hphi,
                               double* x)
{
   assert(!std::isnan(x[0]));

   const double epsilon = 1e-1;

   const double cA = d_model_parameters.getStochioA();
   const double cB = d_model_parameters.getStochioB();

   const double xeq = d_model_parameters.ceq_liquid(temp);
   // std::cout << "xeq = " << xeq << std::endl;
   assert(!std::isnan(xeq));

   // solve explicit KKS problem
   const double a = hphi[0] + epsilon * std::exp(-hphi[0] / epsilon);
   const double dkks = (c[0] - hphi[0] * xeq - hphi[1] * cA - hphi[2] * cB) / a;

   const double xmin = 0.7857;
   const double xmax = 0.999;
   const double d = dkks > 0. ? (xmax - xeq) : xeq - xmin;
   const double t = std::tanh(dkks / d);
   const double xkks = xeq + d * t;

#if 0
   for (short i = 0; i < 3; i++)
      std::cerr << hphi[i] << ", ";
   std::cerr << ", c=" << c[0] << std::endl;
   std::cerr << "x = " << x[0] << std::endl;
   std::cerr << "xkks = " << xkks << std::endl;
   std::cerr << "conc[0] - hphi1 * cA - hphi2 * cB = "
             << c[0] - hphi[1] * cA - hphi[2] * cB << std::endl;
#endif

   x[0] = xkks;
   x[1] = cA;
   x[2] = cB;

   return 0;
}
