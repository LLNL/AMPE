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
#include "CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB.h"

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

CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB::
    CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB(
        boost::property_tree::ptree calphad_db,
        std::shared_ptr<tbox::Database> newton_db,
        const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
        const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
        MolarVolumeStrategy* mvstrategy, const int conc_l_id,
        const int conc_a_id, const int conc_b_id)
    : CALPHADFreeEnergyStrategyBinaryFolchPlapp<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>(
          calphad_db, newton_db, energy_interp_func_type, conc_interp_func_type,
          mvstrategy, conc_l_id, conc_a_id, conc_b_id)
{
   setup(calphad_db, newton_db);
}

//=======================================================================

void CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB::setup(
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


bool CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB::computeCeqT(
    const double temperature, const Thermo4PFM::PhaseIndex pi0,
    const Thermo4PFM::PhaseIndex pi1, double* ceq)
{
   std::cout << "CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB::"
                "computeCeqT..."
             << std::endl;
   return d_ceq_fenergy->computeCeqT(temperature, &ceq[0], 50, true);
}
