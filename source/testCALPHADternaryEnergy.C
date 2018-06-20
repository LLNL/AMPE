#include "CALPHADFreeEnergyFunctionsTernary.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Database.h"

#include <boost/make_shared.hpp>

#include <string>
#include <fstream>

using namespace SAMRAI;
using namespace std;


int main( int argc, char *argv[] )
{
   // Initialize MPI, SAMRAI, and enable logging.

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

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

   boost::shared_ptr<tbox::MemoryDatabase> input_db(new tbox::MemoryDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   std::string run_name = input_filename.substr( 0, input_filename.rfind( "." ) );

   std::string log_file_name = run_name + ".log";
   tbox::PIO::logOnlyNodeZero( log_file_name );

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) tbox::plog<<" AMPE git version "<<xstr(x)<<endl;
    LOG(GITVERSION);
    tbox::plog<<endl;
#endif

   tbox::plog << "Run with "<<mpi.getSize()<<" MPI tasks"<<endl; 	
   tbox::plog << "input_filename = " << input_filename << endl;

   boost::shared_ptr<tbox::Database> model_db =
      input_db->getDatabase("ModelParameters");

   double phase_well_scale = 0.;

   string phase_well_func_type = "double";

   string energy_interp_func_type = "pbg";
   string conc_interp_func_type = "pbg";

   boost::shared_ptr<tbox::Database> temperature_db =
      model_db->getDatabase( "Temperature" );
   double temperature = temperature_db->getDouble( "temperature" );

   boost::shared_ptr<tbox::Database> conc_db(
      model_db->getDatabase( "ConcentrationModel" ));
   string conc_avg_func_type = "a";

   boost::shared_ptr<tbox::Database> dcalphad_db=
      conc_db->getDatabase( "Calphad" );
   std::string calphad_filename = dcalphad_db->getString( "filename" );
   boost::shared_ptr<tbox::MemoryDatabase> calphad_db (
      new tbox::MemoryDatabase( "calphad_db" ) );
   tbox::InputManager::getManager()->parseInputFile(
      calphad_filename, calphad_db );
   boost::shared_ptr<tbox::Database> newton_db;  
 
   CALPHADFreeEnergyFunctionsTernary
      cafe(calphad_db, newton_db,
           energy_interp_func_type,
           conc_interp_func_type,
           conc_avg_func_type,
           phase_well_scale,
           phase_well_func_type);
  

   double c[2];
   model_db->getDoubleArray("concentration",&c[0],2);

   // choose pair of phases: phaseL, phaseA
   PHASE_INDEX pindex;
   string phase = model_db->getString("phase");
   if( phase=="solid" )pindex=phaseA;
   else if(phase=="liquid" )pindex=phaseL;
   else{
      cerr<<"ERROR: Phase needs to be 'solid' or liquid'"<<endl;
      tbox::SAMRAI_MPI::abort();
   }

   double energy = cafe.computeFreeEnergy(temperature,&c[0],pindex);

   cout<<setprecision(12);
   cout<<"Energy: "<<energy<<endl;

   input_db.reset();

   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return(0);
}
