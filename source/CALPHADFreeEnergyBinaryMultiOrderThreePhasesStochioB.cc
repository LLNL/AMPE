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
#include "CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB.h"

#include <boost/property_tree/json_parser.hpp>
#include "Database2JSON.h"
namespace pt = boost::property_tree;

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

using namespace SAMRAI;

#include <cassert>

//=======================================================================

CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB::
    CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB(
        const int norderp_A, boost::property_tree::ptree calphad_pt,
        std::shared_ptr<tbox::Database> newton_db,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        MolarVolumeStrategy* mvstrategy, const int conc_l_id,
        const int conc_a_id, const int conc_b_id)
    : CALPHADFreeEnergyBinaryMultiOrderThreePhases<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>(
          calphad_pt, newton_db, conc_interp_func_type, norderp_A, mvstrategy,
          conc_l_id, conc_a_id, conc_b_id)
{
   tbox::plog << "CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB..."
              << std::endl;
   setup(calphad_pt, newton_db);

   d_multiorder_driving_force.reset(
       new MultiOrderBinaryThreePhasesDrivingForceStochioB(this, norderp_A));
}

//=======================================================================

void CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB::setup(
    pt::ptree calphad_pt, std::shared_ptr<tbox::Database> newton_db)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);
   // set looser tol in solver
   newton_pt.put("tol", 1.e-6);
   d_ceq_fenergy.reset(new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(
       calphad_pt, newton_pt, d_energy_interp_func_type,
       d_conc_interp_func_type));
}


bool CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB::computeCeqT(
    const double temperature, const Thermo4PFM::PhaseIndex pi0,
    const Thermo4PFM::PhaseIndex pi1, double* ceq)
{
   std::cout << "CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB::"
                "computeCeqT..."
             << std::endl;
   return d_ceq_fenergy->computeCeqT(temperature, &ceq[0], 50, true);
}
