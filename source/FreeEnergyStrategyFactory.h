#ifndef included_FreeEnergyStrategyFactory
#define included_FreeEnergyStrategyFactory

#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFreeEnergyStrategyTernary.h"
#include "CALPHADFreeEnergyStrategyWithPenalty.h"
#include "KKSdiluteBinary.h"
#include "HBSMFreeEnergyStrategy.h"
#include "BiasDoubleWellBeckermannFreeEnergyStrategy.h"
#include "BiasDoubleWellUTRCFreeEnergyStrategy.h"

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
            std::shared_ptr<tbox::MemoryDatabase> newton_db;

            tbox::pout << "QuatModel: "
                       << "Using CALPHAD model for concentration" << std::endl;
            std::shared_ptr<tbox::Database> db(conc_db->getDatabase("Calphad"));
            std::string calphad_filename = db->getString("filename");
            calphad_db.reset(new tbox::MemoryDatabase("calphad_db"));
            tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                             calphad_db);

            if (conc_db->isDatabase("NewtonSolver")) {
               db = conc_db->getDatabase("NewtonSolver");
               newton_db.reset(new tbox::MemoryDatabase("newton_db"));
            }

            if (ncompositions == 1) {
               free_energy_strategy.reset(new CALPHADFreeEnergyStrategyBinary(
                   calphad_db, newton_db,
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(), mvstrategy,
                   conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                   model_parameters.with_third_phase()));
            } else {
               assert(ncompositions == 2);
               free_energy_strategy.reset(new CALPHADFreeEnergyStrategyTernary(
                   calphad_db, newton_db,
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(), mvstrategy,
                   conc_l_scratch_id, conc_a_scratch_id));
            }
            if (!calphad_db->keyExists("PenaltyPhaseL")) {

#ifndef HAVE_THERMO4PFM
            } else {
               tbox::plog << "QuatModel: "
                          << "Adding penalty to CALPHAD energy" << std::endl;

               assert(ncompositions == 1);

               free_energy_strategy.reset(
                   new CALPHADFreeEnergyStrategyWithPenalty(
                       calphad_db, newton_db,
                       model_parameters.energy_interp_func_type(),
                       model_parameters.conc_interp_func_type(), mvstrategy,
                       conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                       ncompositions, model_parameters.with_third_phase()));
#endif
            }
         }

         else if (model_parameters.isConcentrationModelKKSdilute()) {
            tbox::pout << "QuatModel: "
                       << "Using KKS dilute model for concentration"
                       << std::endl;
            free_energy_strategy.reset(new KKSdiluteBinary(
                conc_db, model_parameters.energy_interp_func_type(),
                model_parameters.conc_interp_func_type(), mvstrategy,
                conc_l_scratch_id, conc_a_scratch_id));
         } else if (model_parameters.isConcentrationModelHBSM()) {
            tbox::pout << "QuatModel: "
                       << "Using HBSM model for concentration" << std::endl;
            free_energy_strategy.reset(new HBSMFreeEnergyStrategy(
                conc_db->getDatabase("HBSM"),
                model_parameters.energy_interp_func_type(),
                model_parameters.molar_volume_liquid(),
                model_parameters.molar_volume_solid_A(),
                model_parameters.molar_volume_solid_B(),
                model_parameters.D_liquid(), model_parameters.D_solid_A(),
                model_parameters.D_solid_B(), model_parameters.Q0_liquid(),
                model_parameters.Q0_solid_A(), model_parameters.Q0_solid_B(),
                conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
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
