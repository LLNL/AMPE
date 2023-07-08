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
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Box.h"

#include "computeQDiffs.h"

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
         box_upper(d) = (d + 2) * 3;
      }

      hier::Box box(box_lower, box_upper, hier::BlockId(0));

      std::shared_ptr<pdat::CellData<double>> quat(
          new pdat::CellData<double>(box, qlen, hier::IntVector(dim, NGHOSTS)));
      std::shared_ptr<pdat::SideData<double>> quat_diffs(
          new pdat::SideData<double>(box, qlen, hier::IntVector(dim, NGHOSTS)));

      // initialize quat fields as linear functions of x,y,z
      std::vector<double> alpha(NDIM * qlen);
      for (int axis = 0; axis < dim.getValue(); ++axis)
         for (int q = 0; q < qlen; q++)
            alpha[axis * qlen + q] = (double)(axis + 1) + 0.1 * (double)q;

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
            (*quat)(cell, q) = ix * alpha[0 * qlen + q] +
                               iy * alpha[1 * qlen + q]
#if (NDIM == 3)
                               + iz * alpha[2 * qlen + q]
#endif
                ;
         }
      }

      // compute diffs
      computeQDiffs(quat, quat_diffs, false, nullptr);

      // verify result
      for (int axis = 0; axis < dim.getValue(); ++axis) {
         tbox::pout << "Verify component " << axis << " of diff Q..."
                    << std::endl;
         pdat::SideIterator iend(pdat::SideGeometry::end(box, axis));
         for (pdat::SideIterator si(pdat::SideGeometry::begin(box, axis));
              si != iend; ++si) {
            pdat::SideIndex side = *si;
            for (int q = 0; q < qlen; q++) {
               if (((*quat_diffs)(side, q) - alpha[axis * qlen + q]) > 1.e-6) {
                  tbox::pout << "expected diff = " << alpha[axis * qlen + q]
                             << ", computed = " << (*quat_diffs)(side, q)
                             << std::endl;
                  ret = 1;
               }
            }
         }
      }

      // compute gradients at cell centers
      tbox::pout << "Gradients at cell centers..." << std::endl;
      double dx[3] = {0.1, 0.11, 0.12};
      std::shared_ptr<pdat::CellData<double>> quat_grad(
          new pdat::CellData<double>(box, NDIM * qlen,
                                     hier::IntVector(dim, 0)));
      computeQGrad(quat_diffs, quat_grad, dx, false, nullptr);

      // verify result
      for (int axis = 0; axis < dim.getValue(); ++axis) {
         tbox::pout << "Verify component " << axis << " of grad Q..."
                    << std::endl;
         pdat::CellIterator iend(pdat::CellGeometry::end(box));
         for (pdat::CellIterator ci(pdat::CellGeometry::begin(box)); ci != iend;
              ++ci) {
            pdat::CellIndex cell = *ci;
            for (int q = 0; q < qlen; q++) {
               double val = (*quat_grad)(cell, q + axis * qlen);
               double expected = alpha[axis * qlen + q] / dx[axis];
               if ((val - expected) > 1.e-6) {
                  tbox::pout << "q=" << q << ": expected grad = " << expected
                             << ", computed = " << val << std::endl;
                  ret = 1;
               }
            }
         }
      }

      // compute gradients at cell sides
      tbox::pout << "Gradients at cell sides..." << std::endl;
      std::shared_ptr<pdat::SideData<double>> quat_grad_side(
          new pdat::SideData<double>(box, NDIM * qlen,
                                     hier::IntVector(dim, 0)));
      computeQGradSide(quat_diffs, quat_grad_side, dx, false, false, nullptr);

      // verify result
      for (int axis = 0; axis < dim.getValue(); ++axis) {
         tbox::pout << "Verify component " << axis << " of grad Q..."
                    << std::endl;
         pdat::SideIterator iend(pdat::SideGeometry::end(box, axis));
         for (pdat::SideIterator si(pdat::SideGeometry::begin(box, axis));
              si != iend; ++si) {
            pdat::SideIndex side = *si;
            for (int q = 0; q < qlen; q++) {
               double val = (*quat_grad_side)(side, q + axis * qlen);
               double expected = alpha[axis * qlen + q] / dx[axis];
               if ((val - expected) > 1.e-6) {
                  tbox::pout << "q=" << q << ": expected grad = " << expected
                             << ", computed = " << val << std::endl;
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
