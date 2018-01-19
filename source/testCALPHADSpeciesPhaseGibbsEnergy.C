#include "CALPHADSpeciesPhaseGibbsEnergy.h"

#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   cout<<"Test CALPHAD Species Gibbs Energy functions."<<endl;

   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();


   std::string input_filename;
   input_filename = argv[1];
   cout<<"Input file = "<<input_filename<<endl;

   boost::shared_ptr<tbox::MemoryDatabase> input_db =
      tbox::InputManager::getManager()->parseInputFile(input_filename);
   input_db->printClassData(cout);

   if ( input_db->keyExists( "ModelParameters" ) ){
      boost::shared_ptr<tbox::Database> model_db =
         input_db->getDatabase("ModelParameters");

      cout<<"Read T"<<endl;
      boost::shared_ptr<tbox::Database> temperature_db = model_db->getDatabase( "Temperature" );

      double temperature = temperature_db->getDouble( "temperature" );
      cout<<"T="<<temperature<<endl;

      boost::shared_ptr<tbox::Database> conc_db(model_db->getDatabase( "ConcentrationModel" ));
      boost::shared_ptr<tbox::Database> dcalphad_db=conc_db->getDatabase( "Calphad" );
      std::string calphad_filename = dcalphad_db->getString( "filename" );
      boost::shared_ptr<tbox::MemoryDatabase> calphad_db ( new tbox::MemoryDatabase( "calphad_db" ) );
      tbox::InputManager::getManager()->parseInputFile( calphad_filename, calphad_db );

      boost::shared_ptr<tbox::Database> speciesC_db = calphad_db->getDatabase( "SpeciesC" );
      speciesC_db->printClassData(cout);

      string name = speciesC_db->getStringWithDefault( "name", "unknown" );
      
      CALPHADSpeciesPhaseGibbsEnergy gspecies;
      gspecies.initialize(name,speciesC_db->getDatabase( "PhaseL" ) );

      double energy = gspecies.fenergy(temperature);
      cout<<"Energy = "<<energy<<endl;
   }else{
      cout<<"No ModelParameters found in input file"<<endl;
      return 1;
  }

   return(0);
}
