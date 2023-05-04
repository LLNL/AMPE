// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#ifndef included_FreeEnergyStrategyFactory
#define included_FreeEnergyStrategyFactory

#ifdef HAVE_THERMO4PFM
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "CALPHADFunctions.h"
#endif

#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFreeEnergyStrategyTernary.h"
#include "CALPHADFreeEnergyStrategyWithPenalty.h"
#include "CALPHADFreeEnergyStrategyBinaryThreePhase.h"
#include "KKSdiluteBinary.h"
#include "QuadraticFreeEnergyStrategy.h"
#include "BiasDoubleWellBeckermannFreeEnergyStrategy.h"
#include "BiasDoubleWellUTRCFreeEnergyStrategy.h"
#include "DeltaTemperatureFreeEnergyStrategy.h"
#include "TiltingFolchPlapp2005.h"
#include "TiltingMoelans2011.h"
#include "thermo/InterpolationType.h"

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
      std::shared_ptr<FreeEnergyStrategy> free_energy_strategy;

      if (model_parameters.with_concentration()) {

         if (model_parameters.isConcentrationModelCALPHAD()) {
            std::shared_ptr<tbox::MemoryDatabase> calphad_db;
            boost::property_tree::ptree calphad_pt;

            tbox::pout << "QuatModel: "
                       << "Using CALPHAD model for concentration" << std::endl;
            std::shared_ptr<tbox::Database> db(conc_db->getDatabase("Calphad"));
            std::string calphad_filename = db->getString("filename");
            bool calphad_file_is_json = false;
            if (calphad_filename.compare(calphad_filename.size() - 4, 4,
                                         "json") == 0) {
               boost::property_tree::read_json(calphad_filename, calphad_pt);
               calphad_file_is_json = true;
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

            if (ncompositions == 1) {
#ifdef HAVE_THERMO4PFM
               if (conc_b_scratch_id >= 0) {
                  tbox::pout << "Use CALPHADFreeEnergyStrategyBinaryThreePhase"
                             << std::endl;
                  // check if sublattice parameters are in CALPHAD database
                  bool subl = checkSublattice(calphad_pt);
                  if (subl) {
                     tbox::plog << "CALPHADFreeEnergyFunctionsBinary3Ph2Sl"
                                << std::endl;
                     if (model_parameters
                             .energy_three_args_interp_func_type() ==
                         EnergyThreeArgsInterpolationType::MOELANS2011)
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyStrategyBinaryThreePhase<
                                CALPHADFreeEnergyFunctionsBinary3Ph2Sl,
                                TiltingMoelans2011>(
                                calphad_pt, newton_db,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id));
                     else
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyStrategyBinaryThreePhase<
                                CALPHADFreeEnergyFunctionsBinary3Ph2Sl,
                                TiltingFolchPlapp2005>(
                                calphad_pt, newton_db,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id));
                  } else {
                     tbox::plog << "CALPHADFreeEnergyFunctionsBinaryThreePhase"
                                << std::endl;
                     if (model_parameters
                             .energy_three_args_interp_func_type() ==
                         EnergyThreeArgsInterpolationType::MOELANS2011)
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyStrategyBinaryThreePhase<
                                CALPHADFreeEnergyFunctionsBinaryThreePhase,
                                TiltingMoelans2011>(
                                calphad_pt, newton_db,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id));
                     else  // default
                        free_energy_strategy.reset(
                            new CALPHADFreeEnergyStrategyBinaryThreePhase<
                                CALPHADFreeEnergyFunctionsBinaryThreePhase,
                                TiltingFolchPlapp2005>(
                                calphad_pt, newton_db,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                mvstrategy, conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id));
                  }
                  // conc_b_scratch_id<0
               } else
#endif
               {
#ifdef HAVE_THERMO4PFM
                  // check if sublattice parameters are in CALPHAD database
                  bool subl = checkSublattice(calphad_pt);
                  if (subl) {
                     tbox::plog << "CALPHADFreeEnergyFunctionsBinary2Ph1Sl"
                                << std::endl;
                     free_energy_strategy.reset(
                         new CALPHADFreeEnergyStrategyBinary<
                             CALPHADFreeEnergyFunctionsBinary2Ph1Sl>(
                             calphad_pt, newton_db,
                             model_parameters.energy_interp_func_type(),
                             model_parameters.conc_interp_func_type(),
                             mvstrategy, conc_l_scratch_id, conc_a_scratch_id,
                             conc_b_scratch_id, false));
                  } else
#endif
                     free_energy_strategy.reset(
                         new CALPHADFreeEnergyStrategyBinary<
                             CALPHADFreeEnergyFunctionsBinary>(
#ifdef HAVE_THERMO4PFM
                             calphad_pt,
#else
                          calphad_db,
#endif
                             newton_db,
                             model_parameters.energy_interp_func_type(),
                             model_parameters.conc_interp_func_type(),
                             mvstrategy, conc_l_scratch_id, conc_a_scratch_id,
                             conc_b_scratch_id,
                             model_parameters.with_third_phase()));
               }
            } else {  // ncompositions!=1
               assert(ncompositions == 2);
               free_energy_strategy.reset(new CALPHADFreeEnergyStrategyTernary(
                   calphad_db, newton_db,
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(), mvstrategy,
                   conc_l_scratch_id, conc_a_scratch_id));
            }
#ifndef HAVE_THERMO4PFM
            if (!calphad_file_is_json) {
               if (calphad_db->keyExists("PenaltyPhaseL")) {
                  tbox::plog << "QuatModel: "
                             << "Adding penalty to CALPHAD energy" << std::endl;

                  assert(ncompositions == 1);

                  free_energy_strategy.reset(
                      new CALPHADFreeEnergyStrategyWithPenalty(
                          calphad_db, newton_db,
                          model_parameters.energy_interp_func_type(),
                          model_parameters.conc_interp_func_type(), mvstrategy,
                          conc_l_scratch_id, conc_a_scratch_id,
                          conc_b_scratch_id, ncompositions,
                          model_parameters.with_third_phase()));
               }
            }
#endif
         }
         // not CALPHAD
         else if (model_parameters.isConcentrationModelKKSdilute()) {
            tbox::pout << "QuatModel: "
                       << "Using KKS dilute model for concentration"
                       << std::endl;
            free_energy_strategy.reset(new KKSdiluteBinary(
                conc_db, model_parameters.energy_interp_func_type(),
                model_parameters.conc_interp_func_type(), mvstrategy,
                conc_l_scratch_id, conc_a_scratch_id));
         } else if (model_parameters.isConcentrationModelQuadratic()) {
            tbox::pout << "QuatModel: "
                       << "Using Quadratic model for concentration"
                       << std::endl;
            free_energy_strategy.reset(new QuadraticFreeEnergyStrategy(
                conc_db->getDatabase("Quadratic"),
                model_parameters.energy_interp_func_type(),
                model_parameters.molar_volume_liquid(),
                model_parameters.molar_volume_solid_A(),
                model_parameters.molar_volume_solid_B(), conc_l_scratch_id,
                conc_a_scratch_id, conc_b_scratch_id,
                model_parameters.with_third_phase()));
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
                model_parameters.eta_interp_func_type(),
                model_parameters.free_energy_solid_A(),
                model_parameters.free_energy_solid_B(),
                model_parameters.molar_volume_solid_A(),
                model_parameters.molar_volume_solid_B(),
                model_parameters.latent_heat(), Tref,
                model_parameters.with_third_phase()));

      } else {  // no composition, no heat equation
         if (model_parameters.free_energy_type()[0] == 's') {
            free_energy_strategy.reset(new PhaseFreeEnergyStrategy(
                model_parameters.energy_interp_func_type(),
                model_parameters.eta_interp_func_type(),
                model_parameters.free_energy_liquid(),
                model_parameters.free_energy_solid_A(),
                model_parameters.free_energy_solid_B(),
                model_parameters.molar_volume_liquid(),
                model_parameters.molar_volume_solid_A(),
                model_parameters.molar_volume_solid_B(),
                model_parameters.with_third_phase()));
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
