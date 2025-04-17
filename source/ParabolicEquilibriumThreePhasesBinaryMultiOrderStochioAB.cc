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
#include "ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB.h"
#include "FuncFort.h"
#include "ParabolicTools.h"

ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB::
    ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB(
        const short norderp_A, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const QuatModelParameters& model_parameters,
        std::shared_ptr<tbox::Database> conc_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases(
          norderp_A, conc_l_id, conc_a_id, conc_b_id, model_parameters,
          conc_db),
      d_model_parameters(model_parameters)
{
}

ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB::
    ~ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB(){};

int ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB::
    computePhaseConcentrations(const double temp, double* c, double* hphi,
                               double* x)
{
   assert(!std::isnan(x[0]));

   const double epsilon1 = 1.e-4;
   const double epsilon2 = 2.e-4;

   // initialize to NaN to trigger error if used when not set
   double xkks = tbox::IEEE::getSignalingNaN();
   double xeq = tbox::IEEE::getSignalingNaN();

   const double cA = d_model_parameters.getStochioA();
   const double cB = d_model_parameters.getStochioB();

   if (hphi[0] >= epsilon1) {
      // solve explicit KKS problem
      xkks = (c[0] - hphi[1] * cA - hphi[2] * cB) / hphi[0];

#if 0
      for (short i = 0; i < 3; i++)
         std::cerr << hphi[i] << ", ";
      std::cerr << ", c=" << c[0] << std::endl;
      std::cerr << "x = " << x[0] << std::endl;
      std::cerr << "xkks = " << xkks << std::endl;
      std::cerr << "conc[0] - hphi1 * cA - hphi2 * cB = "
                << c[0] - hphi[1] * cA - hphi[2] * cB << std::endl;
#endif
   }

   if (hphi[0] < epsilon2) {
      xeq = d_model_parameters.ceq_liquid(temp);
      // std::cout << "xeq = " << xeq << std::endl;
      assert(!std::isnan(xeq));
   }

   if (hphi[0] >= epsilon2) {
      x[0] = xkks;
   } else if (hphi[0] <= epsilon1) {
      x[0] = xeq;
   } else {
      // map epsilon1 < hphi < epsilon2 to (0,1)
      double h = (hphi[0] - epsilon1) / (epsilon2 - epsilon1);
      // mix equilibrium compositions and KKS solution
      double f =
          Thermo4PFM::interp_func(Thermo4PFM::EnergyInterpolationType::PBG, h);
      x[0] = (1. - f) * xeq + f * xkks;
   }
   x[1] = cA;
   x[2] = cB;

   return 0;
}
