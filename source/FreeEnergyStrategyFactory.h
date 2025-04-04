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
#ifndef included_FreeEnergyStrategyFactory
#define included_FreeEnergyStrategyFactory

#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "CALPHADFunctions.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFreeEnergyStrategyTernary.h"
#include "CALPHADFreeEnergyBinaryMultiOrder.h"
#include "CALPHADFreeEnergyStrategyBinaryFolchPlapp.h"
#include "CALPHADFreeEnergyBinaryMultiOrderThreePhases.h"
#include "CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB.h"
#include "CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB.h"
#include "ParabolicFreeEnergyBinary.h"
#include "ParabolicFreeEnergyMultiOrderBinary.h"
#include "QuadraticFreeEnergyBinary.h"
#include "QuadraticFreeEnergyMultiOrderBinary.h"
#include "QuadraticFreeEnergyMultiOrderTernaryThreePhase.h"
#include "ParabolicFreeEnergyMultiOrderBinaryThreePhase.h"
#include "KKSdiluteBinary.h"
#include "BiasDoubleWellBeckermannFreeEnergyStrategy.h"
#include "BiasDoubleWellUTRCFreeEnergyStrategy.h"
#include "DeltaTemperatureFreeEnergyStrategy.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

class FreeEnergyStrategyFactory
{
 public:
   static std::shared_ptr<FreeEnergyStrategy> create(
       QuatModelParameters& model_parameters, const int ncompositions,
       const int conc_l_scratch_id, const int conc_a_scratch_id,
       const int conc_b_scratch_id, MolarVolumeStrategy* mvstrategy,
       MeltingTemperatureStrategy* meltingT_strategy, const double Tref,
       std::shared_ptr<tbox::Database> conc_db)
   {
      tbox::plog << "FreeEnergyStrategyFactory: ncompositions = "
                 << ncompositions << std::endl;

      std::shared_ptr<FreeEnergyStrategy> free_energy_strategy;

      if (model_parameters.with_concentration()) {

         if (model_parameters.isConcentrationModelCALPHAD()) {
            std::shared_ptr<tbox::MemoryDatabase> calphad_db;
            boost::property_tree::ptree calphad_pt;

            tbox::pout << "QuatModel: "
                       << "Using CALPHAD model for concentration" << std::endl;
            std::shared_ptr<tbox::Database> db(conc_db->getDatabase("Calphad"));
            std::string calphad_filename = db->getString("filename");
            if (calphad_filename.compare(calphad_filename.size() - 4, 4,
                                         "json") == 0) {
               boost::property_tree::read_json(calphad_filename, calphad_pt);
            } else {
               calphad_db.reset(new tbox::MemoryDatabase("calphad_db"));
               tbox::pout << "FreeEnergyStrategyFactory: Read "
                          << calphad_filename << std::endl;
               tbox::InputManager::getManager()->parseInputFile(
                   calphad_filename, calphad_db);
               copyDatabase(calphad_db, calphad_pt);
            }

            std::shared_ptr<tbox::Database> newton_db;
            if (conc_db->isDatabase("NewtonSolver")) {
               newton_db = conc_db->getDatabase("NewtonSolver");
            }

            if (ncompositions == 1) {  // binary case
               tbox::plog << "ncompositions: 1" << std::endl;
               if (model_parameters.withMultipleOrderP()) {
                  tbox::plog << "MultiOrder..." << std::endl;
                  if (conc_b_scratch_id >= 0) {
                     if (model_parameters.getStochioB() >= 0.) {
                        tbox::plog << "StochioB..." << std::endl;
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB(
                                model_parameters.norderpA(), calphad_pt,
                                newton_db,
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id));
                     } else {
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyBinaryMultiOrderThreePhases<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinaryThreePhase>(
                                calphad_pt, newton_db,
                                model_parameters.conc_interp_func_type(),
                                model_parameters.norderpA(), mvstrategy,
                                conc_l_scratch_id, conc_a_scratch_id,
                                conc_b_scratch_id));
                     }
                  } else {
                     free_energy_strategy.reset(
                         new CALPHADFreeEnergyBinaryMultiOrder(
                             calphad_pt, newton_db,
                             model_parameters.energy_interp_func_type(),
                             model_parameters.conc_interp_func_type(),
                             mvstrategy, conc_l_scratch_id, conc_a_scratch_id));
                  }
               } else {
                  // check if sublattice parameters are in CALPHAD database
                  bool subl = Thermo4PFM::checkSublattice(calphad_pt);
                  if (subl) tbox::plog << "CALPHAD sublattice..." << std::endl;
                  if (conc_b_scratch_id >= 0) {
                     if (subl) {
                        tbox::plog << "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.."
                                      "."
                                   << std::endl;
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyStrategyBinaryFolchPlapp<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinary3Ph2Sl>(
                                calphad_pt, newton_db,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id));
                     } else {
                        tbox::plog << "CALPHADFreeEnergyFunctionsBinaryThreePha"
                                      "se"
                                   << std::endl;
                        if (model_parameters.getStochioB()) {
                           tbox::plog << "Stochio..." << std::endl;
                           free_energy_strategy.reset(
                               new CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB(
                                   calphad_pt, newton_db,
                                   model_parameters.energy_interp_func_type(),
                                   model_parameters.conc_interp_func_type(),
                                   mvstrategy, conc_l_scratch_id,
                                   conc_a_scratch_id, conc_b_scratch_id));
                        } else {
                           free_energy_strategy.reset(
                               new CALPHADFreeEnergyStrategyBinaryFolchPlapp<
                                   Thermo4PFM::
                                       CALPHADFreeEnergyFunctionsBinaryThreePhase>(
                                   calphad_pt, newton_db,
                                   model_parameters.energy_interp_func_type(),
                                   model_parameters.conc_interp_func_type(),
                                   mvstrategy, conc_l_scratch_id,
                                   conc_a_scratch_id, conc_b_scratch_id));
                        }
                     }
                     // conc_b_scratch_id<0
                  } else {
                     if (subl) {
                        tbox::plog << "CALPHADFreeEnergyFunctionsBinary2Ph1Sl"
                                   << std::endl;
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyStrategyBinary<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinary2Ph1Sl>(
                                calphad_pt, newton_db,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id, false));
                     } else {
                        tbox::plog << "CALPHADFreeEnergyStrategyBinary..."
                                   << std::endl;
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyStrategyBinary<
                                Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>(
                                calphad_pt, newton_db,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id,
                                model_parameters.with_third_phase()));
                     }
                  }
               }
            } else {  // ncompositions!=1
               tbox::plog << "CALPHADFreeEnergyStrategyTernary..." << std::endl;
               assert(ncompositions == 2);
               free_energy_strategy.reset(new CALPHADFreeEnergyStrategyTernary(
                   calphad_db, newton_db,
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(), mvstrategy,
                   conc_l_scratch_id, conc_a_scratch_id));
            }
         }
         // not CALPHAD
         else if (model_parameters.isConcentrationModelKKSdilute()) {
            tbox::plog << "Using KKS dilute model for concentration"
                       << std::endl;
            free_energy_strategy.reset(new KKSdiluteBinary(
                conc_db, model_parameters.energy_interp_func_type(),
                model_parameters.conc_interp_func_type(), mvstrategy,
                conc_l_scratch_id, conc_a_scratch_id));
         } else if (model_parameters.isConcentrationModelParabolic()) {
            tbox::plog << "Using Parabolic model for concentration"
                       << std::endl;
            if (model_parameters.norderp() > 1) {
               if (conc_b_scratch_id > -1) {
                  tbox::plog << "ParabolicFreeEnergyMultiOrderBinaryThreePhase."
                                ".."
                             << std::endl;
                  free_energy_strategy.reset(
                      new ParabolicFreeEnergyMultiOrderBinaryThreePhase(
                          conc_db->getDatabase("Parabolic"),
                          model_parameters.conc_interp_func_type(),
                          model_parameters.norderpA(), mvstrategy,
                          conc_l_scratch_id, conc_a_scratch_id,
                          conc_b_scratch_id));
               } else {
                  tbox::plog << "ParabolicFreeEnergyMultiOrderBinary"
                             << std::endl;
                  free_energy_strategy.reset(
                      new ParabolicFreeEnergyMultiOrderBinary(
                          model_parameters.energy_interp_func_type(),
                          model_parameters.conc_interp_func_type(), mvstrategy,
                          conc_l_scratch_id, conc_a_scratch_id, conc_db));
               }
            } else {
               tbox::plog << "ParabolicFreeEnergyBinary" << std::endl;
               free_energy_strategy.reset(new ParabolicFreeEnergyBinary(
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(), mvstrategy,
                   conc_l_scratch_id, conc_a_scratch_id, conc_db));
            }
         } else if (model_parameters.isConcentrationModelQuadratic()) {
            tbox::plog << "Using Quadratic model for concentration"
                       << std::endl;
            if (model_parameters.norderp() > 1) {
               if (conc_b_scratch_id > -1) {
                  tbox::plog << "QuadraticFreeEnergyStrategyMultiOrderTernaryTh"
                                "reePhase..."
                             << std::endl;
                  free_energy_strategy.reset(
                      new QuadraticFreeEnergyMultiOrderTernaryThreePhase(
                          conc_db->getDatabase("Quadratic"),
                          model_parameters.energy_interp_func_type(),
                          model_parameters.norderpA(),
                          model_parameters.molar_volume_liquid(),
                          model_parameters.molar_volume_solid_A(),
                          model_parameters.molar_volume_solid_B(),
                          conc_l_scratch_id, conc_a_scratch_id,
                          conc_b_scratch_id));
               } else {
                  tbox::plog << "QuadraticFreeEnergyStrategyMultiOrder..."
                             << std::endl;
                  free_energy_strategy.reset(
                      new QuadraticFreeEnergyMultiOrderBinary(
                          conc_db->getDatabase("Quadratic"),
                          model_parameters.energy_interp_func_type(),
                          model_parameters.conc_interp_func_type(),
                          model_parameters.molar_volume_liquid(),
                          model_parameters.molar_volume_solid_A(),
                          conc_l_scratch_id, conc_a_scratch_id));
               }
            } else {
               tbox::plog << "QuadraticFreeEnergyStrategy" << std::endl;
               free_energy_strategy.reset(new QuadraticFreeEnergyBinary(
                   conc_db->getDatabase("Quadratic"),
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(),
                   model_parameters.molar_volume_liquid(),
                   model_parameters.molar_volume_solid_A(), conc_l_scratch_id,
                   conc_a_scratch_id));
            }
         } else if (model_parameters.with_bias_well()) {
            if (model_parameters.wellBiasBeckermann()) {
               free_energy_strategy.reset(
                   new BiasDoubleWellBeckermannFreeEnergyStrategy(
                       model_parameters.well_bias_alpha(), meltingT_strategy));
            } else {
               free_energy_strategy.reset(
                   new BiasDoubleWellUTRCFreeEnergyStrategy(
                       model_parameters.well_bias_alpha(),
                       model_parameters.well_bias_gamma(), meltingT_strategy));
            }
         }
      } else if (model_parameters.with_heat_equation()) {
         if (model_parameters.with_bias_well()) {
            free_energy_strategy.reset(new BiasDoubleWellUTRCFreeEnergyStrategy(
                model_parameters.well_bias_alpha(),
                model_parameters.well_bias_gamma(), meltingT_strategy));
         } else if (model_parameters.free_energy_type()[0] == 'l') {
            free_energy_strategy.reset(new DeltaTemperatureFreeEnergyStrategy(
                Tref, model_parameters.latent_heat(),
                model_parameters.energy_interp_func_type()));
         } else
            free_energy_strategy.reset(new TemperatureFreeEnergyStrategy(
                model_parameters.energy_interp_func_type(),
                model_parameters.free_energy_solid_A(),
                model_parameters.molar_volume_solid_A(),
                model_parameters.latent_heat(), Tref));
      } else {  // no composition, no heat equation
         if (model_parameters.free_energy_type()[0] == 's') {
            tbox::plog << "PhaseFreeEnergyStrategy..." << std::endl;
            free_energy_strategy.reset(new PhaseFreeEnergyStrategy(
                model_parameters.energy_interp_func_type(),
                model_parameters.free_energy_liquid(),
                model_parameters.free_energy_solid_A(),
                model_parameters.molar_volume_liquid(),
                model_parameters.molar_volume_solid_A()));
         }
      }
      // pure element free energy
      if (model_parameters.free_energy_type()[0] == 'l') {
         free_energy_strategy.reset(new DeltaTemperatureFreeEnergyStrategy(
             Tref, model_parameters.latent_heat(),
             model_parameters.energy_interp_func_type()));
      }

      return free_energy_strategy;
   }
};
#endif
