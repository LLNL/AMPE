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
#ifndef included_ParabolicFreeEnergyMultiOrderBinary
#define included_ParabolicFreeEnergyMultiOrderBinary

#include "ParabolicFreeEnergyFunctionsBinary.h"
#include "FreeEnergyStrategyBinary.h"
#include "InterpolationType.h"
#include "MultiOrderBinaryDrivingForce.h"
#include "MolarVolumeStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

class ParabolicFreeEnergyMultiOrderBinary : public FreeEnergyStrategyBinary
{
 public:
   ParabolicFreeEnergyMultiOrderBinary(
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       MolarVolumeStrategy* mvstrategy, const int conc_l_id,
       const int conc_a_id, std::shared_ptr<tbox::Database> conc_db);

   ~ParabolicFreeEnergyMultiOrderBinary(){};

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int f_l_id, const int f_a_id,
                        const int f_b_id, const int rhs_id) override;

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

   void computeSecondDerivativeEnergyPhase(const char phase,
                                           const std::vector<double>& c,
                                           std::vector<double>& d2fdc2,
                                           const bool use_internal_units)
   {
      switch (phase) {
         case 'l':
            computeSecondDerivativeEnergyPhaseL(-1., c, d2fdc2,
                                                use_internal_units);
            break;

         case 'a':
            computeSecondDerivativeEnergyPhaseA(-1., c, d2fdc2,
                                                use_internal_units);
            break;

         default:
            tbox::pout << "undefined phase=" << phase << "!!!" << std::endl;
            tbox::SAMRAI_MPI::abort();
      }
   }

   void preRunDiagnostics(const double temperature){};


 private:
   MolarVolumeStrategy* d_mv_strategy;

   std::shared_ptr<Thermo4PFM::ParabolicFreeEnergyFunctionsBinary>
       d_parabolic_fenergy;

   std::shared_ptr<MultiOrderBinaryDrivingForce> d_multiorder_driving_force;

   double computeFreeEnergy(const double temperature, double* c_i,
                            const Thermo4PFM::PhaseIndex pi, const bool gp);

   double computeDerivFreeEnergy(const double temperature, double* c_i,
                                 const Thermo4PFM::PhaseIndex pi);

   double computeMuL(const double t, const double c);
   double computeMuA(const double t, const double c);
};

#endif
