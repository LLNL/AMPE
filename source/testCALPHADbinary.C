#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADMobility.h"
#include "CompositionStrategyMobilities.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "ConstantMolarVolumeStrategy.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include <boost/make_shared.hpp>
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>

using namespace SAMRAI;
using namespace std;


int main( int argc, char *argv[] )
{
   // Initialize MPI, SAMRAI

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /* This extra code block is used to scope some temporaries that are
    * created, it forces the destruction before the manager is
    * shutdown.
    */
   {

   //-----------------------------------------------------------------------
   /*
    * Process command line arguments and dump to log file.
    * For non-restarted case, command line is:
    *
    *    executable <input file name>
    *
    */

   std::string input_filename;
   input_filename = argv[1];

   //-----------------------------------------------------------------------
   // Create input database and parse all data in input file.

   boost::shared_ptr<tbox::MemoryDatabase> input_db(
      new tbox::MemoryDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   // make from input file name
   std::string run_name =
      input_filename.substr( 0, input_filename.rfind( "." ) );

   // Logfile
   std::string log_file_name = run_name + ".log";
   tbox::PIO::logOnlyNodeZero( log_file_name );

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) tbox::plog<<" AMPE: git version "<<xstr(x)<<endl;
    LOG(GITVERSION);
    tbox::plog<<endl;
#endif

   tbox::plog << "input_filename = " << input_filename << endl;

   boost::shared_ptr<tbox::Database> model_db =
      input_db->getDatabase("ModelParameters");

   double phase_well_scale = model_db->getDouble( "phi_well_scale" );

   string phase_well_func_type =
         model_db->getString( "phi_well_func_type" );

   string energy_interp_func_type = "pbg";
   string conc_interp_func_type = "pbg";
   string eta_interp_func_type   ="pbg";
   
   boost::shared_ptr<tbox::Database> temperature_db =
      model_db->getDatabase( "Temperature" );
   double temperature = temperature_db->getDouble( "temperature" );

   boost::shared_ptr<tbox::Database> conc_db(
      model_db->getDatabase( "ConcentrationModel" ));
   string conc_avg_func_type =
      conc_db->getStringWithDefault( "avg_func_type", "a" );

   boost::shared_ptr<tbox::Database> dcalphad_db=
      conc_db->getDatabase( "Calphad" );
   std::string calphad_filename = dcalphad_db->getString( "filename" );
   boost::shared_ptr<tbox::MemoryDatabase> calphad_db (
      new tbox::MemoryDatabase( "calphad_db" ) );
   tbox::InputManager::getManager()->parseInputFile(
      calphad_filename, calphad_db );
   
   boost::shared_ptr<tbox::Database> newton_db;
   if ( conc_db->isDatabase( "NewtonSolver" ) )
      newton_db = conc_db->getDatabase( "NewtonSolver" );

   bool with_third_phase=false;
   
   CALPHADFreeEnergyFunctionsBinary
      cafe(calphad_db, newton_db,
           energy_interp_func_type,
           conc_interp_func_type,
           eta_interp_func_type,
           conc_avg_func_type,
           with_third_phase,
           phase_well_scale,
           0.,
           phase_well_func_type,
           "");
   
   cafe.printEnergyVsComposition(temperature);

   cafe.preRunDiagnostics(303., 2899. );

   // choose pair of phases: phaseL, phaseA, phaseB
   const PHASE_INDEX pi0=phaseL;
   const PHASE_INDEX pi1=phaseA;
   
   // initial guesses
   double init_guess[2];
   model_db->getDoubleArray("initial_guess", &init_guess[0], 2);

   double lceq[2]={init_guess[0], init_guess[1]};
   
   // compute equilibrium concentrations
   bool found_ceq =
      cafe.computeCeqT(temperature,pi0,pi1,&lceq[0]);
   if( lceq[0]>1. )found_ceq = false;
   if( lceq[0]<0. )found_ceq = false;
   if( lceq[1]>1. )found_ceq = false;
   if( lceq[1]<0. )found_ceq = false;
   
   if( found_ceq ){
      tbox::pout<<"Found equilibrium concentrations: "
                <<lceq[0]<<" and "<<lceq[1]<<"..."<<endl;
   }else{
      tbox::pout<<"WARNING: Equilibrium concentrations not found... "<<endl;
   }

   double molar_volume = conc_db->getDouble( "molar_volume" );

   tbox::plog<<"ConstantMolarVolumeStrategy... "<<endl;
   ConstantMolarVolumeStrategy mvstrategy(
      molar_volume,
      molar_volume,
      molar_volume);
   tbox::plog<<"CALPHADFreeEnergyStrategy... "<<endl;
   CALPHADFreeEnergyStrategyBinary free_energy_strategy(
               calphad_db, newton_db,
               energy_interp_func_type,
               conc_interp_func_type,
               "",
               conc_avg_func_type,
               &mvstrategy,
               -1,-1,-1,
               false,
               phase_well_scale,
               0.,
               phase_well_func_type,
               "" );

   if( calphad_db->keyExists( "MobilityParameters" ) ){
      tbox::plog<<"CompositionStrategyMobilities... "<<endl;
      int ncompositions=1;
      bool with_third_phase=false;
      CompositionStrategyMobilities composition_strategy_mobilities(
         dcalphad_db, with_third_phase, ncompositions, &free_energy_strategy );

      composition_strategy_mobilities.printDiagnostics(temperature,temperature);
   }
   cafe.energyVsPhiAndC(temperature, &lceq[0], found_ceq, with_third_phase,
                        101, 100);

   input_db.reset();

   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return(0);
}
