#include "ParabolicEqConcSolverBinary.h"
#include "ParabolicTools.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <boost/optional/optional.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;
using namespace SAMRAI;

int main(int argc, char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   int ret = 0;
   {
      std::string databasename(argv[1]);
      double temperature = atof(argv[2]);
      std::cout << "Temperature = " << temperature << std::endl;
      double c0 = 0.5;
      double c1 = 1.5;
      if (argc > 3) {
         c0 = atof(argv[3]);
         c1 = atof(argv[4]);
      }

      std::shared_ptr<tbox::MemoryDatabase> input_db(
          new tbox::MemoryDatabase("db"));
      std::cout << "Filename = " << databasename << std::endl;
      tbox::InputManager::getManager()->parseInputFile(databasename, input_db);

      double coeffL[3][2];
      readParabolicData(input_db, "Liquid", coeffL);

      double coeffA[3][2];
      readParabolicData(input_db, "PhaseA", coeffA);

      double coeffB[3][2];
      readParabolicData(input_db, "PhaseB", coeffB);

      double Tref = input_db->getDouble("Tref");

      input_db->printClassData(std::cout);

      Thermo4PFM::ParabolicEqConcSolverBinary solver;
      solver.setup(temperature - Tref, coeffL, coeffA);

      const double toln = 1.e-12;
      const int max_iters = 100;
      const double alpha = 1.;

      double sol[2] = {c0, c1};
      std::cout << "Phases L,A..." << std::endl;
      int ret = solver.ComputeConcentration(sol, toln, max_iters, alpha);
      std::cout << ret << " iterations" << std::endl;
      std::cout << "Solution: " << sol[0] << " " << sol[1] << std::endl;

      sol[0] = c0;
      sol[1] = c1;

      std::cout << "Phases L,B..." << std::endl;
      solver.setup(temperature - Tref, coeffL, coeffB);
      ret = solver.ComputeConcentration(sol, toln, max_iters, alpha);
      std::cout << ret << " iterations" << std::endl;
      std::cout << "Solution: " << sol[0] << " " << sol[1] << std::endl;

      sol[0] = c0;
      sol[1] = c1;

      std::cout << "Phases A,B..." << std::endl;
      solver.setup(temperature - Tref, coeffA, coeffB);
      ret = solver.ComputeConcentration(sol, toln, max_iters, alpha);
      std::cout << ret << " iterations" << std::endl;
      std::cout << "Solution: " << sol[0] << " " << sol[1] << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
