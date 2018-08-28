#ifndef included_MobilityFactory
#define included_MobilityFactory

#include "KimMobilityStrategyInfMob.h"
#include "KimMobilityStrategyFiniteMob.h"
#include "QuatModelParameters.h"

class QuatModel;

class MobilityFactory
{
public:
   static boost::shared_ptr<QuatMobilityStrategy> create(
             QuatModel* model,
             QuatModelParameters& model_parameters,
             const int conc_l_scratch_id,
             const int conc_a_scratch_id,
             const int temperature_scratch_id,
             const int ncompositions,
             boost::shared_ptr<tbox::Database> conc_db)
   {
      boost::shared_ptr<QuatMobilityStrategy> mobility_strategy;

      if( model_parameters.isPhaseMobilityScalar() )
         mobility_strategy.reset( new SimpleQuatMobilityStrategy( model ) );
      else{
         //no support for composition dependent diffusion for now
         assert( !model_parameters.conDiffusionStrategyIsCTD() );

         boost::shared_ptr<tbox::Database> conc_calphad_db=
            conc_db->getDatabase( "Calphad" );
         string calphad_filename = conc_calphad_db->getString( "filename" );
         boost::shared_ptr<tbox::MemoryDatabase> calphad_db
            ( new tbox::MemoryDatabase( "calphad_db" ) );
         tbox::InputManager::getManager()->parseInputFile(
            calphad_filename, calphad_db );

         boost::shared_ptr<tbox::Database> newton_db=
            conc_db->getDatabase( "NewtonSolver" );

         if( model_parameters.interfaceMobility()>0. ){
            mobility_strategy.reset(
               new KimMobilityStrategyFiniteMob(
                      model,
                      conc_l_scratch_id, conc_a_scratch_id,
                      temperature_scratch_id,
                      model_parameters.interfaceMobility(),
                      model_parameters.epsilon_phase(),
                      model_parameters.phase_well_scale(),
                      model_parameters.energy_interp_func_type(),
                      model_parameters.conc_interp_func_type(),
                      calphad_db,
                      newton_db,
                      ncompositions,
                      model_parameters.D_liquid(),
                      model_parameters.Q0_liquid(),
                      model_parameters.molar_volume_liquid()) );
         }else{
            mobility_strategy.reset(
               new KimMobilityStrategyInfMob(
                      model,
                      conc_l_scratch_id, conc_a_scratch_id,
                      temperature_scratch_id,
                      model_parameters.epsilon_phase(),
                      model_parameters.phase_well_scale(),
                      model_parameters.energy_interp_func_type(),
                      model_parameters.conc_interp_func_type(),
                      calphad_db,
                      newton_db,
                      ncompositions,
                      model_parameters.D_liquid(),
                      model_parameters.Q0_liquid(),
                      model_parameters.molar_volume_liquid()) );
         }
      }
      return mobility_strategy;
   }

};

#endif

