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
#ifndef included_CALPHADFreeEnergyStrategyBinary
#define included_CALPHADFreeEnergyStrategyBinary

#include "FreeEnergyStrategyBinary.h"
#include "FuncFort.h"
#include "InterpolationType.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
class MolarVolumeStrategy;

#include <boost/property_tree/ptree.hpp>

#include <string>
#include <vector>

template <class FreeEnergyFunctionType>
class CALPHADFreeEnergyStrategyBinary : public FreeEnergyStrategyBinary
{
 public:
   CALPHADFreeEnergyStrategyBinary(
       boost::property_tree::ptree calphad_db,
       std::shared_ptr<tbox::Database> newton_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       MolarVolumeStrategy* mvstrategy, const int conc_l_id,
       const int conc_a_id, const int conc_b_id, const bool with_third_phase);

   virtual ~CALPHADFreeEnergyStrategyBinary(){};

   virtual void setup(boost::property_tree::ptree calphad_db,
                      std::shared_ptr<tbox::Database> newton_db);

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

   MolarVolumeStrategy* d_mv_strategy;

   std::shared_ptr<FreeEnergyFunctionType> d_calphad_fenergy;

   double computeMuA(const double t, const double c);
   double computeMuL(const double t, const double c);
   double computeMuB(const double t, const double c);

 private:
   double computeFreeEnergy(const double temperature, double* c_i,
                            const Thermo4PFM::PhaseIndex pi, const bool gp);

   double computeDerivFreeEnergy(const double temperature, double* c_i,
                                 const Thermo4PFM::PhaseIndex pi);
};

#endif
