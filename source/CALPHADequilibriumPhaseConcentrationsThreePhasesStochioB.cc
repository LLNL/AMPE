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
namespace pt = boost::property_tree;

#include "CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

// Thermo4PFM
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB ::
    CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB(
        QuatModelParameters& model_parameters, const int conc_l_scratch_id,
        const int conc_a_scratch_id, const int conc_b_scratch_id,
        const int conc_l_ref_id, const int conc_a_ref_id,
        const int conc_b_ref_id,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db,
        const unsigned ncompositions)
    : EquilibriumPhaseConcentrationsThreePhases<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>(
          conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
          conc_l_ref_id, conc_a_ref_id, conc_b_ref_id, energy_interp_func_type,
          Thermo4PFM::ConcInterpolationType::LINEAR, calphad_pt, newton_db,
          ncompositions),
      d_model_parameters(model_parameters)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   d_calphad_fenergy = std::unique_ptr<
       Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>(
       new Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB(
           model_parameters.getStochioB(), calphad_pt, newton_pt,
           energy_interp_func_type, Thermo4PFM::ConcInterpolationType::LINEAR));
}

int CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB ::
    computeAuxilliaryConcentrations(const double temp, double* c, double* hphi,
                                    double* x)
{
   const double epsilon1 = 1.e-2;
   const double epsilon2 = 2.e-2;

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
         std::cerr << ", c=" << c[0] << ", " << c[1] << ", " << c[2]
                   << std::endl;
         std::cerr << "xkks = " << xkks[0] << ", " << xkks[1] << std::endl;
         const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
         MPI_Abort(mpi.getCommunicator(), -1);
      }
#endif
   }

   if (hphi[2] > 1. - epsilon2) {
      xeq[0] = d_model_parameters.ceq_liquid(temp);
      xeq[1] = d_model_parameters.ceq_solidA(temp);
      // std::cout << "xeq[0] = " << xeq[0] << std::endl;
      // std::cout << "xeq[1] = " << xeq[1] << std::endl;
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
