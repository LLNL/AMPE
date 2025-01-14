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
#include "CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder.h"
#include "Database2JSON.h"
#include "FuncFort.h"

#include "SAMRAI/tbox/InputManager.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder::
    CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder(
        const int conc_l_id, const int conc_a_id,
        std::shared_ptr<tbox::Database> conc_db,
        std::shared_ptr<tbox::Database> newton_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrder(conc_l_id, conc_a_id,
                                                     conc_db)
{
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

   d_fenergy.reset(new Thermo4PFM::CALPHADFreeEnergyFunctionsBinary(
       calphad_pt, newton_pt, Thermo4PFM::EnergyInterpolationType::LINEAR,
       Thermo4PFM::ConcInterpolationType::LINEAR));
}
