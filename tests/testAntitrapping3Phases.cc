// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "ConcFort.h"

#include "SAMRAI/SAMRAI_config.h"

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


#include <vector>
#include <iomanip>

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
         box_upper(d) = 3 + d;
      }

      hier::Box box(box_lower, box_upper, hier::BlockId(0));

      std::shared_ptr<pdat::CellData<double>> cl(
          new pdat::CellData<double>(box, 1, hier::IntVector(dim, 1)));
      std::shared_ptr<pdat::CellData<double>> ca(
          new pdat::CellData<double>(box, 1, hier::IntVector(dim, 1)));
      std::shared_ptr<pdat::CellData<double>> cb(
          new pdat::CellData<double>(box, 1, hier::IntVector(dim, 1)));

      std::shared_ptr<pdat::CellData<double>> phi(
          new pdat::CellData<double>(box, 3, hier::IntVector(dim, 1)));
      std::shared_ptr<pdat::CellData<double>> dphidt(
          new pdat::CellData<double>(box, 3, hier::IntVector(dim, 1)));

      // initialize phi fields as linear functions of x,y,z
      hier::Box gbox(phi->getGhostBox());

      const hier::Index& ifirst = box.lower();
      const hier::Index& ilast = box.upper();

      // different dx for 3 directions should lead to same
      // different results for 3 directions since used only
      // for computing normal vector
      double dx[3] = {0.1, 0.11, 0.12};

      // loop over solid phases
      for (int ip = 1; ip < 3; ip++) {
         std::cout << "Solid phase: " << ip << std::endl;
         // test for solid front facing one direction at a time
         for (int d = 0; d < dim.getValue(); ++d) {

            pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox));
                 i != iend; ++i) {
               pdat::CellIndex cell = *i;
               double delta = (double)(cell(d) - box.lower(d)) /
                              (double)(box.upper(d) - box.lower(d));
               (*phi)(cell, 0) = 0.5 * (1. + tanh(delta));
               (*phi)(cell, 1) = 0.;
               (*phi)(cell, 2) = 0.;
               (*phi)(cell, ip) = 1. - (*phi)(cell, 0);

               (*dphidt)(cell, 0) = 0.88;  // arbitrary value
               (*dphidt)(cell, 1) = 0.88;
               (*dphidt)(cell, 2) = 0.88;

               (*cl)(cell) = 0.1;  // arbitrary value
               (*ca)(cell) = 0.2;
               (*cb)(cell) = 0.2;
            }

            std::shared_ptr<pdat::SideData<double>> flux(
                new pdat::SideData<double>(box, 1, hier::IntVector(dim, 0)));
            std::shared_ptr<pdat::SideData<double>> flux3(
                new pdat::SideData<double>(box, 1, hier::IntVector(dim, 0)));

            double alpha = 1.5;  // arbitrary value

            // compute antitrapping for 1 phi model
            ADDCONCENTRATIONFLUXFROMANTITRAPPING(
                ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                dx, phi->getPointer(), phi->getGhostCellWidth()[0],
                cl->getPointer(), ca->getPointer(), cl->getGhostCellWidth()[0],
                1, dphidt->getPointer(), dphidt->getGhostCellWidth()[0], alpha,
                flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                flux->getPointer(2),
#endif
                flux->getGhostCellWidth()[0]);

            // compute antitrapping for 3 phis model
            ADDCONCENTRATIONFLUXFROMANTITRAPPING3PHASES(
                ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                dx, phi->getPointer(), phi->getGhostCellWidth()[0],
                cl->getPointer(), ca->getPointer(), cb->getPointer(),
                cl->getGhostCellWidth()[0], dphidt->getPointer(),
                dphidt->getGhostCellWidth()[0], alpha, flux3->getPointer(0),
                flux3->getPointer(1),
#if (NDIM == 3)
                flux3->getPointer(2),
#endif
                flux3->getGhostCellWidth()[0]);

            // verify result
            std::cout << "Verify result for solid front in direction " << d
                      << std::endl;
            std::cout << std::setprecision(10);

            for (int dd = 0; dd < NDIM; dd++) {
               std::cout << "Side: " << dd << std::endl;
               pdat::SideIterator send(pdat::SideGeometry::end(box, dd));
               for (pdat::SideIterator si(pdat::SideGeometry::begin(box, dd));
                    si != send; ++si) {
                  pdat::SideIndex side = *si;
                  double diff = (*flux)(side) - (*flux3)(side);
                  std::cout << "1 phase: " << (*flux)(side)
                            << ", 3 phases: " << (*flux3)(side) << std::endl;
                  if (std::abs(diff) > 1.e-6) {
                     std::cout << "Computed diff = " << diff << std::endl;
                     ret = 1;
                  }
               }
            }
         }
      }
   }

   if (ret == 0) std::cout << "\nPASSED" << std::endl;

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return ret;
}
