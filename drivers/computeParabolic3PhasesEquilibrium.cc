#include "ParabolicEqConcSolverBinary.h"

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
      std::shared_ptr<tbox::Database> liquid_db =
          input_db->getDatabase("Liquid");
      coeffL[0][0] = liquid_db->getDouble("a0");
      coeffL[0][1] = liquid_db->getDouble("a1");
      coeffL[1][0] = liquid_db->getDouble("b0");
      coeffL[1][1] = liquid_db->getDouble("b1");
      coeffL[2][0] = liquid_db->getDouble("c0");
      coeffL[2][1] = liquid_db->getDouble("c1");

      double coeffA[3][2];
      std::shared_ptr<tbox::Database> phasea_db =
          input_db->getDatabase("PhaseA");
      coeffA[0][0] = phasea_db->getDouble("a0");
      coeffA[0][1] = phasea_db->getDouble("a1");
      coeffA[1][0] = phasea_db->getDouble("b0");
      coeffA[1][1] = phasea_db->getDouble("b1");
      coeffA[2][0] = phasea_db->getDouble("c0");
      coeffA[2][1] = phasea_db->getDouble("c1");

      double coeffB[3][2];
      std::shared_ptr<tbox::Database> phaseb_db =
          input_db->getDatabase("PhaseB");
      coeffB[0][0] = phaseb_db->getDouble("a0");
      coeffB[0][1] = phaseb_db->getDouble("a1");
      coeffB[1][0] = phaseb_db->getDouble("b0");
      coeffB[1][1] = phaseb_db->getDouble("b1");
      coeffB[2][0] = phaseb_db->getDouble("c0");
      coeffB[2][1] = phaseb_db->getDouble("c1");

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
