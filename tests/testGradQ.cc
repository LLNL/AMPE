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
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Box.h"

#include "computeQDiffs.h"


using namespace SAMRAI;


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
      const int qlen = 4;

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
         box_upper(d) = (d + 4) * 3;
      }

      hier::Box box(box_lower, box_upper, hier::BlockId(0));

      std::shared_ptr<pdat::CellData<double>> quat(
          new pdat::CellData<double>(box, qlen, hier::IntVector(dim, 1)));
      std::shared_ptr<pdat::SideData<double>> quat_diffs(
          new pdat::SideData<double>(box, qlen, hier::IntVector(dim, 1)));

      // initialize quat field
      const double alpha[3] = {
          2.,
          3.,
          4,
      };

      hier::Box gbox(quat->getGhostBox());

      int x_lower = gbox.lower(0);
      int y_lower = gbox.lower(1);
#if (NDIM == 3)
      int z_lower = gbox.lower(2);
#endif

      pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
      for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox)); i != iend;
           ++i) {
         pdat::CellIndex cell = *i;
         int ix = cell(0) - x_lower;
         int iy = cell(1) - y_lower;
#if (NDIM == 3)
         int iz = cell(2) - z_lower;
#endif
         for (int q = 0; q < qlen; q++) {
            (*quat)(cell, q) = ix * alpha[0] + iy * alpha[1]
#if (NDIM == 3)
                               + iz * alpha[2]
#endif
                ;
         }
      }

      // compute diffs
      computeQDiffs(quat, quat_diffs, false, nullptr);

      // verify result
      for (int axis = 0; axis < dim.getValue(); ++axis) {
         pdat::SideIterator iend(pdat::SideGeometry::end(box, axis));
         for (pdat::SideIterator si(pdat::SideGeometry::begin(box, axis));
              si != iend; ++si) {
            pdat::SideIndex side = *si;
            for (int q = 0; q < qlen; q++) {
               if (((*quat_diffs)(side, q) - alpha[axis]) > 1.e-6) {
                  tbox::pout << "expected diff = " << alpha[axis]
                             << ", computed = " << (*quat_diffs)(side, q)
                             << std::endl;
                  ret = 1;
               }
            }
         }
      }
   }

   if (ret == 0) tbox::pout << "\nPASSED" << std::endl;


   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
