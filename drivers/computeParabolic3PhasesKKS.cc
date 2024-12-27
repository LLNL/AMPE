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
      double phiL = atof(argv[3]);
      double phiA = atof(argv[4]);
      double phiB = atof(argv[5]);
      double c = atof(argv[6]);

      std::cout << "Temperature = " << temperature << std::endl;

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

      Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase
          parabolic_fenergy(Tref, coeffL, coeffA, coeffB,
                            Thermo4PFM::EnergyInterpolationType::PBG,
                            Thermo4PFM::ConcInterpolationType::LINEAR);

      double phi[3] = {phiL, phiA, phiB};
      double conc[3] = {c, c, c};
      parabolic_fenergy.computePhaseConcentrations(temperature, &c, phi, conc);

      std::cout << "phi = (" << phiL << "," << phiA << "," << phiB << ")"
                << ", c=" << c << std::endl;
      std::cout << "KKS solution: cl = " << conc[0] << ", ca=" << conc[1]
                << " and cb = " << conc[2] << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
