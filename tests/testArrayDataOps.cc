// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Box.h"

#include "ArrayOperation.h"

#include <vector>

using namespace SAMRAI;

#define NGHOSTS (1)

int main(int argc, char* argv[])
{
   /*
    * Initialize MPI, SAMRAI.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   int ret = 0;
   {
      const tbox::Dimension dim(static_cast<unsigned short>(NDIM));

      /*
       * Start logging.
       */
      const std::string log_file_name = "test.log";
      tbox::PIO::logOnlyNodeZero(log_file_name);


      // create a rectangular box
      hier::Index box_lower(dim, 0);
      hier::Index box_upper(dim);
      for (int d = 0; d < dim.getValue(); ++d) {
         box_upper(d) = (d + 2) * 3;
      }

      hier::Box box(box_lower, box_upper, hier::BlockId(0));
      tbox::pout << "Box :" << box << std::endl;
      std::shared_ptr<pdat::SideData<double>> sdata1(
          new pdat::SideData<double>(box, 1, hier::IntVector(dim, NGHOSTS)));
      sdata1->fill(0., box);

      const double value = 0.2;
      std::shared_ptr<pdat::SideData<double>> sdata2(
          new pdat::SideData<double>(box, 1, hier::IntVector(dim, NGHOSTS)));
      sdata2->fill(value, box);

      const hier::IntVector src_shift(box.getDim(), 0);
      AddOperation<double> addop;
      const unsigned int num_depth = 1;
      const unsigned int src_depth = 0;
      const unsigned int dst_depth = 0;

      // add sdata2 to sdata1 using SAMRAI function "doArrayDataOperationOnBox"
      for (int side_normal = 0; side_normal < dim.getValue(); side_normal++) {

         hier::Box sbox(box);
         hier::Index up(box.upper());
         up[side_normal]++;
         sbox.setUpper(up);

         pdat::ArrayDataOperationUtilities<double, AddOperation<double>>::
             doArrayDataOperationOnBox(sdata1->getArrayData(side_normal),
                                       sdata2->getArrayData(side_normal), sbox,
                                       src_shift, dst_depth, src_depth,
                                       num_depth, addop);
      }

      // verify result
      for (int side_normal = 0; side_normal < dim.getValue(); ++side_normal) {
         tbox::pout << "Verify component " << side_normal << " of SideData..."
                    << std::endl;
         pdat::SideIterator iend(pdat::SideGeometry::end(box, side_normal));
         int count = 0;
         for (pdat::SideIterator si(
                  pdat::SideGeometry::begin(box, side_normal));
              si != iend; ++si) {
            pdat::SideIndex side = *si;
            if (std::abs((*sdata1)(side)-value) > 1.e-6) {
               tbox::pout << "expected value = " << value
                          << ", extracted = " << (*sdata1)(side) << std::endl;
               ret = 1;
            }
            count++;
         }
         tbox::pout << "Verified " << count << " values" << std::endl;
      }
   }

   if (ret == 0) tbox::pout << "\nPASSED" << std::endl;

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
