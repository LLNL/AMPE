#include "AMPE.h"
#include "QuatModel.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"

using namespace std;

AMPE::AMPE(MPI_Comm comm)
{
   tbox::SAMRAI_MPI::init(comm);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   tbox::plog<<"AMPE: git_version "<<gitCommitID()<<endl;

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   tbox::plog << "Run with "<<mpi.getSize()<<" MPI tasks"<<endl;
}

AMPE::~AMPE()
{
   delete d_pfm;

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
}

void AMPE::initialize(const string input_filename,
                      const string restart_read_dirname,
                      const int restore_num)
{
   bool is_from_restart = (restore_num>=0);

   boost::shared_ptr<tbox::MemoryDatabase> input_db(
      new tbox::MemoryDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   string model_type = input_db->getStringWithDefault( "model_type", "Quat" );

   string run_name;
   if ( input_db->keyExists( "run_name" ) ) {
      run_name = input_db->getString( "run_name" );
   }
   else {
      // make from input file name
      run_name = input_filename.substr( 0, input_filename.rfind( "." ) );
   }

   bool log_all_nodes = false;
   string log_file_name = run_name + ".log";

   if ( input_db->isDatabase( "Logging" ) ) {
      boost::shared_ptr<tbox::Database> log_db = 
         input_db->getDatabase( "Logging" );
      if ( log_db->keyExists( "filename" ) ) {
         log_file_name = log_db->getString( "filename" );
      }

      if ( log_db->keyExists( "log_all_nodes" ) ) {
         log_all_nodes = log_db->getBool( "log_all_nodes" );
      }
   }

   if ( log_all_nodes ) {
      tbox::PIO::logAllNodes( log_file_name );
   }
   else {
      tbox::PIO::logOnlyNodeZero( log_file_name );
   }

   if ( is_from_restart ) {
      tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
      tbox::plog << "restore_num = " << restore_num << endl;
   }

   // Create timers
   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

   d_time_man = tbox::TimerManager::getManager();
   d_time_man->resetAllTimers();

   // Create a PFModel object
   if ( model_type == "Quat" ) {
      d_pfm = new QuatModel( 4 );
   }
#if NDIM==2
   else if ( model_type == "KWC" ) {
      d_pfm = new QuatModel( 1 );
   }
   else if ( model_type == "KWCcomplex" ) {
      d_pfm = new QuatModel( 2 );
   }
#endif
   else {
      TBOX_ERROR( "Invalid model_type" << endl );
   }

   d_pfm->Initialize(
      input_db, run_name,
      is_from_restart, restart_read_dirname, restore_num
      );
   /*
    * After creating all objects and initializing their state, we
    * print the input database and variable database contents to
    * the log file.
    */
   tbox::plog << "\nCheck input data and variables before simulation:" << endl;
   tbox::plog << "Input database..." << endl;
   input_db->printClassData(tbox::plog);
   tbox::plog << "\nVariable database..." << endl;
   hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

   input_db.reset();
}

void AMPE::run()
{
   d_pfm->Run();

   d_time_man->print(tbox::plog);
}

