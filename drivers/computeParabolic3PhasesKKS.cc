#include "ParabolicFreeEnergyFunctionsBinaryThreePhase.h"
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
      readParabolicData(input_db, "Liquid", coeffL);

      double coeffA[3][2];
      readParabolicData(input_db, "PhaseA", coeffA);

      double coeffB[3][2];
      readParabolicData(input_db, "PhaseB", coeffB);

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
