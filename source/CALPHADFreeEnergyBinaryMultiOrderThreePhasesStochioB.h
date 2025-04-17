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
#ifndef included_CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB
#define included_CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB

#include "CALPHADFreeEnergyBinaryMultiOrderThreePhases.h"
#include "InterpolationType.h"
#include "MultiOrderBinaryThreePhasesDrivingForceStochioB.h"

#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
class MolarVolumeStrategy;

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

class CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB
    : public CALPHADFreeEnergyBinaryMultiOrderThreePhases<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>
{
 public:
   CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB(
       const int norderp_A, boost::property_tree::ptree calphad_db,
       std::shared_ptr<tbox::Database> newton_db,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       MolarVolumeStrategy* mvstrategy, const int conc_l_id,
       const int conc_a_id, const int conc_b_id);

   bool computeCeqT(const double temperature, const Thermo4PFM::PhaseIndex pi0,
                    const Thermo4PFM::PhaseIndex pi1, double* ceq) override;

 private:
   void setup(boost::property_tree::ptree calphad_pt,
              std::shared_ptr<tbox::Database> newton_db);

   std::shared_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsBinary> d_ceq_fenergy;


   std::shared_ptr<MultiOrderBinaryThreePhasesDrivingForceStochioB>
       d_multiorder_driving_force;
};

#endif
