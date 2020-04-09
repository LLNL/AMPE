#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

using namespace SAMRAI;

int main(int argc, char** argv)
{
   // Initialize tbox::MPI and SAMRAI
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      int np = mpi.getSize();
      tbox::pout << "Run MPI with " << np << " tasks\n";
   }

   // Shutdown SAMRAI and tbox::MPI.
   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   // return 0 for SUCCESS
   return 0;
}
