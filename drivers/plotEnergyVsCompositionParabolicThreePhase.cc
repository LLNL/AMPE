#include "ParabolicFreeEnergyFunctionsBinaryThreePhase.h"

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
      double cmin = 0.;
      double cmax = 1.;
      if (argc > 3) {
         cmin = atof(argv[3]);
         cmax = atof(argv[4]);
      }

      Thermo4PFM::EnergyInterpolationType energy_interp_func_type =
          Thermo4PFM::EnergyInterpolationType::PBG;
      Thermo4PFM::ConcInterpolationType conc_interp_func_type =
          Thermo4PFM::ConcInterpolationType::LINEAR;

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

      const int npts = 100;
      std::cout << " Compute energies..." << std::endl;

      Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase qfe(
          Tref, coeffL, coeffA, coeffB, energy_interp_func_type,
          conc_interp_func_type);

      std::stringstream ss;
      ss << std::fixed;
      ss << std::setprecision(2);
      ss << temperature;
      std::ofstream os("ParabolicFvsC" + ss.str() + ".csv", std::ios::out);

      qfe.printEnergyVsComposition(temperature, os, cmin, cmax, npts);
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
