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
#include "CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "Database2JSON.h"

#include "SAMRAI/tbox/InputManager.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

template <class FreeEnergyType>
CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases<
    FreeEnergyType>::
    CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases(
        const short norderp_A, const int conc_l_id, const int conc_a_id,
        const int conc_b_id, const QuatModelParameters& model_parameters,
        std::shared_ptr<tbox::Database> conc_db,
        std::shared_ptr<tbox::Database> newton_db)
    : EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases(
          norderp_A, conc_l_id, conc_a_id, conc_b_id, model_parameters, conc_db)
{
   tbox::plog << "CALPHADequilPhaseConcMultiOrderThreePhases..." << std::endl;
   tbox::plog << "Type: " << typeid(FreeEnergyType).name() << std::endl;

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

   d_fenergy.reset(
       new FreeEnergyType(calphad_pt, newton_pt,
                          Thermo4PFM::EnergyInterpolationType::LINEAR,
                          model_parameters.conc_interp_func_type()));
}

template class
    CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases<
        Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhase>;
template class
    CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases<
        Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>;
