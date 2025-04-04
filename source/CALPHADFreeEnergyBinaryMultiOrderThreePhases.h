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
#ifndef included_CALPHADFreeEnergyBinaryMultiOrderThreePhases
#define included_CALPHADFreeEnergyBinaryMultiOrderThreePhases

#include "FreeEnergyStrategyBinary.h"
#include "MultiOrderBinaryThreePhasesDrivingForce.h"
#include "MolarVolumeStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

template <class FreeEnergyFunctionType>
class CALPHADFreeEnergyBinaryMultiOrderThreePhases
    : public FreeEnergyStrategyBinary
{
 public:
   CALPHADFreeEnergyBinaryMultiOrderThreePhases(
       boost::property_tree::ptree calphad_pt,
       std::shared_ptr<tbox::Database> newton_db,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       const short norderp_A, MolarVolumeStrategy* mvstrategy,
       const int conc_l_id, const int conc_a_id, const int conc_b_id);

   ~CALPHADFreeEnergyBinaryMultiOrderThreePhases(){};

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int f_l_id, const int f_a_id,
                        const int f_b_id, const int rhs_id) override;

   void preRunDiagnostics(const double temperature)
   {
      d_calphad_fenergy->preRunDiagnostics(temperature);
   }

   bool computeCeqT(const double temperature, const Thermo4PFM::PhaseIndex pi0,
                    const Thermo4PFM::PhaseIndex pi1, double* ceq);

 protected:
   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

   double computeFreeEnergy(const double temperature, double* c_i,
                            const Thermo4PFM::PhaseIndex pi,
                            const bool gp) override;

   double computeDerivFreeEnergy(const double temperature, double* c_i,
                                 const Thermo4PFM::PhaseIndex pi) override;

   double computeMuA(const double t, const double c);
   double computeMuL(const double t, const double c);

 private:
   MolarVolumeStrategy* d_mv_strategy;

   std::shared_ptr<FreeEnergyFunctionType> d_calphad_fenergy;

   std::shared_ptr<MultiOrderBinaryThreePhasesDrivingForce>
       d_multiorder_driving_force;

   void setup(boost::property_tree::ptree calphad_db,
              std::shared_ptr<tbox::Database> newton_db);
};

#endif
