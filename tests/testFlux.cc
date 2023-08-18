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

#include "ConcFort.h"

#include <vector>

using namespace SAMRAI;

#define NGHOSTS (2)

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

      const int depth = 1;

      const double dx[3] = {0.2, 0.15, 0.1};

      std::shared_ptr<pdat::CellData<double>> data(
          new pdat::CellData<double>(box, depth,
                                     hier::IntVector(dim, NGHOSTS)));
      std::shared_ptr<pdat::SideData<double>> diffusion(
          new pdat::SideData<double>(box, depth * depth,
                                     hier::IntVector(dim, NGHOSTS)));

      std::shared_ptr<pdat::SideData<double>> flux(
          new pdat::SideData<double>(box, depth,
                                     hier::IntVector(dim, NGHOSTS)));

      // initialize slopes for data fields
      std::vector<double> alpha(NDIM * depth);
      for (int axis = 0; axis < dim.getValue(); ++axis)
         for (int d = 0; d < depth; d++)
            alpha[axis * depth + d] = (double)(axis + 1) + 0.2 * (double)d;

      hier::Box gbox(data->getGhostBox());

      int x_lower = gbox.lower(0);
      int y_lower = gbox.lower(1);
#if (NDIM == 3)
      int z_lower = gbox.lower(2);
#endif

      // initialize data field as linear functions of x,y,z,
      // including ghost cells
      pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
      for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox)); i != iend;
           ++i) {
         pdat::CellIndex cell = *i;
         int ix = cell(0) - x_lower;
         int iy = cell(1) - y_lower;
#if (NDIM == 3)
         int iz = cell(2) - z_lower;
#endif
         for (int d = 0; d < depth; d++) {
            (*data)(cell, d) = ix * alpha[0 * depth + d] +
                               iy * alpha[1 * depth + d]
#if (NDIM == 3)
                               + iz * alpha[2 * depth + d]
#endif
                ;
         }
      }

      // initialize diffusion
      const double Dcoeff = 111.;
      diffusion->fillAll(Dcoeff);

      // Case 1: touches all physical boundaries
      //      int physbc[6] = {1, 1, 1, 1, 1, 1};
      for (int i = 0; i < 2; i++) {
         // Case i=0: touches no physical boundaries
         // Case i=1: touches all physical boundaries
         int physbc[6] = {i, i, i, i, i, i};

         // compute flux
         flux->fillAll(0.);

         const hier::Index& ifirst(box.lower());
         const hier::Index& ilast(box.upper());

         ADD_FLUX_4TH(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                      ifirst(2), ilast(2),
#endif
                      dx, data->getPointer(), data->getGhostCellWidth()[0],
                      depth, diffusion->getPointer(0), diffusion->getPointer(1),
#if (NDIM == 3)
                      diffusion->getPointer(2),
#endif
                      diffusion->getGhostCellWidth()[0], flux->getPointer(0),
                      flux->getPointer(1),
#if (NDIM == 3)
                      flux->getPointer(2),
#endif
                      flux->getGhostCellWidth()[0], physbc);


         // verify result
         for (int axis = 0; axis < dim.getValue(); ++axis) {
            tbox::pout << "Verify component " << axis << " of flux..."
                       << std::endl;
            pdat::SideIterator iend(pdat::SideGeometry::end(box, axis));
            for (pdat::SideIterator si(pdat::SideGeometry::begin(box, axis));
                 si != iend; ++si) {
               pdat::SideIndex side = *si;
               for (int d = 0; d < depth; d++) {
                  const double eflux =
                      Dcoeff * alpha[axis * depth + d] / dx[axis];
                  if (std::abs((*flux)(side, d) - eflux) > 1.e-6) {
                     tbox::pout << side << ", expected flux = " << eflux
                                << ", computed = " << (*flux)(side, d)
                                << std::endl;
                     ret = 1;
                  }
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
