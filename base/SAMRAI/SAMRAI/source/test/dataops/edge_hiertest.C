/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program to test edge-centered patch data ops
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
// using namespace std;

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/EdgeData.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyEdgeDataOpsComplex.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyEdgeDataOpsReal.h"
#include "SAMRAI/pdat/EdgeIndex.h"
#include "SAMRAI/pdat/EdgeIterator.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

#include <boost/shared_ptr.hpp>

using namespace SAMRAI;

/* Helper function declarations */
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   boost::shared_ptr<hier::PatchHierarchy> hierarchy);

#define NVARS 4

int main(
   int argc,
   char* argv[]) {

   int num_failures = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   if (argc < 2) {
      TBOX_ERROR("Usage: " << argv[0] << " [dimension]");
   }

   const unsigned short d = static_cast<unsigned short>(atoi(argv[1]));
   TBOX_ASSERT(d > 0);
   TBOX_ASSERT(d <= tbox::Dimension::MAXIMUM_DIMENSION_VALUE);
   const tbox::Dimension dim(d);

   if (dim != tbox::Dimension(2)) {
      TBOX_ERROR("This test code is completed only for 2D!!!");
   }

   const std::string log_fn = std::string("edge_hiertest.")
      + tbox::Utilities::intToString(dim.getValue(), 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      const tbox::Dimension dim2d(2);

      int ln, iv;

      // Make a dummy hierarchy domain
      double lo[2] = { 0.0, 0.0 };
      double hi[2] = { 1.0, 0.5 };

      hier::Box coarse0(hier::Index(0, 0), hier::Index(9, 2), hier::BlockId(0));
      hier::Box coarse1(hier::Index(0, 3), hier::Index(9, 4), hier::BlockId(0));
      hier::Box fine0(hier::Index(4, 4), hier::Index(7, 7), hier::BlockId(0));
      hier::Box fine1(hier::Index(8, 4), hier::Index(13, 7), hier::BlockId(0));
      hier::IntVector ratio(dim2d, 2);

      coarse0.initialize(coarse0, hier::LocalId(0), 0);
      coarse1.initialize(coarse1, hier::LocalId(1), 0);
      fine0.initialize(fine0, hier::LocalId(0), 0);
      fine1.initialize(fine1, hier::LocalId(1), 0);

      hier::BoxContainer coarse_domain;
      hier::BoxContainer fine_boxes;
      coarse_domain.pushBack(coarse0);
      coarse_domain.pushBack(coarse1);
      fine_boxes.pushBack(fine0);
      fine_boxes.pushBack(fine1);

      boost::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            lo,
            hi,
            coarse_domain));

      boost::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      hierarchy->setMaxNumberOfLevels(2);
      hierarchy->setRatioToCoarserLevel(ratio, 1);

      // Note: For these simple tests we allow at most 2 processors.
      const int nproc = mpi.getSize();
      TBOX_ASSERT(nproc < 3);

      const int n_coarse_boxes = coarse_domain.size();
      const int n_fine_boxes = fine_boxes.size();

      hier::BoxLevel layer0(hier::IntVector(dim, 1), geometry);
      hier::BoxLevel layer1(ratio, geometry);

      hier::BoxContainer::iterator coarse_itr(coarse_domain);
      for (int ib = 0; ib < n_coarse_boxes; ib++, ++coarse_itr) {
         if (nproc > 1) {
            if (ib == layer0.getMPI().getRank()) {
               layer0.addBox(hier::Box(*coarse_itr, hier::LocalId(ib),
                     layer0.getMPI().getRank()));
            }
         } else {
            layer0.addBox(hier::Box(*coarse_itr, hier::LocalId(ib), 0));
         }
      }

      hier::BoxContainer::iterator fine_itr(fine_boxes);
      for (int ib = 0; ib < n_fine_boxes; ib++, ++fine_itr) {
         if (nproc > 1) {
            if (ib == layer1.getMPI().getRank()) {
               layer1.addBox(hier::Box(*fine_itr, hier::LocalId(ib),
                     layer1.getMPI().getRank()));
            }
         } else {
            layer1.addBox(hier::Box(*fine_itr, hier::LocalId(ib), 0));
         }
      }

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);

      // Create instance of hier::Variable database
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      boost::shared_ptr<hier::VariableContext> dummy(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim2d, 0);

      // Make some dummy variables and data on the hierarchy
      boost::shared_ptr<pdat::EdgeVariable<double> > fvar[NVARS];
      int svindx[NVARS];
      fvar[0].reset(new pdat::EdgeVariable<double>(dim2d, "fvar0", 1));
      svindx[0] = variable_db->registerVariableAndContext(
            fvar[0], dummy, no_ghosts);
      fvar[1].reset(new pdat::EdgeVariable<double>(dim2d, "fvar1", 1));
      svindx[1] = variable_db->registerVariableAndContext(
            fvar[1], dummy, no_ghosts);
      fvar[2].reset(new pdat::EdgeVariable<double>(dim2d, "fvar2", 1));
      svindx[2] = variable_db->registerVariableAndContext(
            fvar[2], dummy, no_ghosts);
      fvar[3].reset(new pdat::EdgeVariable<double>(dim2d, "fvar3", 1));
      svindx[3] = variable_db->registerVariableAndContext(
            fvar[3], dummy, no_ghosts);

      boost::shared_ptr<pdat::EdgeVariable<double> > swgt(
         new pdat::EdgeVariable<double>(dim2d, "swgt", 1));
      int swgt_id = variable_db->registerVariableAndContext(
            swgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(swgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->allocatePatchData(svindx[iv]);
         }
      }

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > edge_ops(
         new math::HierarchyEdgeDataOpsReal<double>(
            hierarchy,
            0,
            1));
      TBOX_ASSERT(edge_ops);

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > swgt_ops(
         new math::HierarchyEdgeDataOpsReal<double>(
            hierarchy,
            0,
            1));

      boost::shared_ptr<hier::Patch> patch;

      // Initialize control volume data for edge-centered components
      hier::Box coarse_fine = fine0 + fine1;
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = *ip;
            boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
               patch->getPatchGeometry(),
               boost::detail::dynamic_cast_tag());
            const double* dx = pgeom->getDx();
            const double edge_vol = dx[0] * dx[1];
            boost::shared_ptr<pdat::EdgeData<double> > data(
               patch->getPatchData(swgt_id),
               boost::detail::dynamic_cast_tag());
            data->fillAll(edge_vol);
            pdat::EdgeIndex fi(dim);
            int plo0 = patch->getBox().lower(0);
            int phi0 = patch->getBox().upper(0);
            int plo1 = patch->getBox().lower(1);
            int phi1 = patch->getBox().upper(1);
            int ic;

            if (ln == 0) {
               data->fillAll(0.0, (coarse_fine * patch->getBox()));

               if (patch->getLocalId() == 0) {
                  //bottom edge boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(ic,
                              plo1), pdat::EdgeIndex::X, pdat::EdgeIndex::Lower);
                     (*data)(fi) *= 0.5;
                  }
                  //left and right edge boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(plo0,
                              ic), pdat::EdgeIndex::Y, pdat::EdgeIndex::Lower);
                     (*data)(fi) *= 0.5;
                     fi = pdat::EdgeIndex(hier::Index(phi0,
                              ic), pdat::EdgeIndex::Y, pdat::EdgeIndex::Upper);
                     (*data)(fi) *= 0.5;
                  }
               } else {
                  //top and bottom edge boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(ic,
                              plo1), pdat::EdgeIndex::X, pdat::EdgeIndex::Lower);
                     (*data)(fi) = 0.0;
                     fi = pdat::EdgeIndex(hier::Index(ic,
                              phi1), pdat::EdgeIndex::X, pdat::EdgeIndex::Upper);
                     (*data)(fi) *= 0.5;
                  }
                  //left and right edge boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(plo0,
                              ic), pdat::EdgeIndex::Y, pdat::EdgeIndex::Lower);
                     (*data)(fi) *= 0.5;
                     fi = pdat::EdgeIndex(hier::Index(phi0,
                              ic), pdat::EdgeIndex::Y, pdat::EdgeIndex::Upper);
                     (*data)(fi) *= 0.5;
                  }
               }
            } else {
               if (patch->getLocalId() == 0) {
                  // top and bottom coarse-fine edge boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(ic,
                              plo1), pdat::EdgeIndex::X, pdat::EdgeIndex::Lower);
                     (*data)(fi) *= 1.5;
                     fi = pdat::EdgeIndex(hier::Index(ic,
                              phi1), pdat::EdgeIndex::X, pdat::EdgeIndex::Upper);
                     (*data)(fi) *= 1.5;
                  }
                  //left coarse-fine edge boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(plo0,
                              ic), pdat::EdgeIndex::Y, pdat::EdgeIndex::Lower);
                     (*data)(fi) *= 1.5;
                  }
               } else {
                  // top and bottom coarse-fine edge boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(ic,
                              plo1), pdat::EdgeIndex::X, pdat::EdgeIndex::Lower);
                     (*data)(fi) *= 1.5;
                     fi = pdat::EdgeIndex(hier::Index(ic,
                              phi1), pdat::EdgeIndex::X, pdat::EdgeIndex::Upper);
                     (*data)(fi) *= 1.5;
                  }
                  //left and right coarse-fine edge boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::EdgeIndex(hier::Index(plo0,
                              ic), pdat::EdgeIndex::Y, pdat::EdgeIndex::Lower);
                     (*data)(fi) = 0.0;
                     fi = pdat::EdgeIndex(hier::Index(phi0,
                              ic), pdat::EdgeIndex::Y, pdat::EdgeIndex::Upper);
                     (*data)(fi) *= 1.5;
                  }
               }
            }
         }
      }

      // Test #1: Print out control volume data and compute its integral

      // Test #1a: Check control volume data set properly
      // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
      // 0.0025 on fine level
/*   bool vol_test_passed = true;
 *   for (ln = 0; ln < 2; ln++) {
 *   for (hier::PatchLevel::iterator ip(hierarchy->getPatchLevel(ln)->begin()); ip != hierarchy->getPatchLevel(ln)->end(); ++ip) {
 *   patch = hierarchy->getPatchLevel(ln)->getPatch(ip());
 *   boost::shared_ptr< pdat::EdgeData<double> > cvdata = patch->getPatchData(cwgt_id);
 *
 *   pdat::EdgeIterator cend(cvdata->getBox(), 1, false);
 *   for (pdat::EdgeIterator c(cvdata->getBox(), 1, true); c != cend && vol_test_passed; ++c) {
 *   pdat::EdgeIndex edge_index = *c;
 *
 *   if (ln == 0) {
 *   if ((coarse_fine * patch->getBox()).contains(edge_index)) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(edge_index),0.0) ) {
 *   vol_test_passed = false;
 *   }
 *   } else {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(edge_index),0.01) ) {
 *   vol_test_passed = false;
 *   }
 *   }
 *   }
 *
 *   if (ln == 1) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(edge_index),0.0025) ) {
 *   vol_test_passed = false;
 *   }
 *   }
 *   }
 *   }
 *   }
 *   if (!vol_test_passed) {
 *   num_failures++;
 *   tbox::perr << "FAILED: - Test #1a: Check control volume data set properly" << std::endl;
 *   cwgt_ops->printData(cwgt_id, tbox::plog);
 *   }
 */
      // Print out control volume data and compute its integral
/*   tbox::plog << "edge control volume data" << std::endl;
 *   swgt_ops->printData(swgt_id, tbox::plog);
 */

      // Test #1b: math::HierarchyEdgeDataOpsReal::sumControlVolumes()
      // Expected: norm = 1.0
      double norm =
         edge_ops->sumControlVolumes(svindx[0], swgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm, 1.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #1b: math::HierarchyEdgeDataOpsReal::sumControlVolumes()\n"
         << "Expected value = 1.0 , Computed value = "
         << norm << std::endl;
      }

      // Test #2: math::HierarchyEdgeDataOpsReal::numberOfEntries()
      // Expected: num_data_points = 209
      int num_data_points = edge_ops->numberOfEntries(svindx[0]);
      if (num_data_points != 209) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #2: math::HierarchyEdgeDataOpsReal::numberOfEntries()\n"
         << "Expected value = 209 , Computed value = "
         << num_data_points << std::endl;
      }

      // Test #3a: math::HierarchyEdgeDataOpsReal::setToScalar()
      // Expected: v0 = 2.0
      double val0 = double(2.0);
      edge_ops->setToScalar(svindx[0], val0);
      if (!doubleDataSameAsValue(svindx[0], val0, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #3a: math::HierarchyEdgeDataOpsReal::setToScalar()\n"
         << "Expected: v0 = " << val0 << std::endl;
         edge_ops->printData(svindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchyEdgeDataOpsReal::setToScalar()
      // Expected: v1 = (4.0)
      edge_ops->setToScalar(svindx[1], 4.0);
      double val1 = 4.0;
      if (!doubleDataSameAsValue(svindx[1], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #3b: math::HierarchyEdgeDataOpsReal::setToScalar()\n"
         << "Expected: v1 = " << val1 << std::endl;
         edge_ops->printData(svindx[1], tbox::plog);
      }

      // Test #4: math::HierarchyEdgeDataOpsReal::copyData()
      // Expected:  v2 = v1 = (4.0)
      edge_ops->copyData(svindx[2], svindx[1]);
      if (!doubleDataSameAsValue(svindx[2], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #4: math::HierarchyEdgeDataOpsReal::copyData()\n"
         << "Expected: v2 = " << val1 << std::endl;
         edge_ops->printData(svindx[2], tbox::plog);
      }

      // Test #5: math::HierarchyEdgeDataOpsReal::swapData()
      // Expected:  v0 = (4.0), v1 = (2.0)
      edge_ops->swapData(svindx[0], svindx[1]);
      if (!doubleDataSameAsValue(svindx[0], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #5a: math::HierarchyEdgeDataOpsReal::swapData()\n"
         << "Expected: v0 = " << val1 << std::endl;
         edge_ops->printData(svindx[0], tbox::plog);
      }
      if (!doubleDataSameAsValue(svindx[1], val0, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #5b: math::HierarchyEdgeDataOpsReal::swapData()\n"
         << "Expected: v1 = " << val0 << std::endl;
         edge_ops->printData(svindx[1], tbox::plog);
      }

      // Test #6: math::HierarchyEdgeDataOpsReal::scale()
      // Expected:  v2 = 0.25 * v2 = (1.0)
      edge_ops->scale(svindx[2], 0.25, svindx[2]);
      double val_scale = 1.0;
      if (!doubleDataSameAsValue(svindx[2], val_scale, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #6: math::HierarchyEdgeDataOpsReal::scale()\n"
         << "Expected: v2 = " << val_scale << std::endl;
         edge_ops->printData(svindx[2], tbox::plog);
      }

      // Test #7: math::HierarchyEdgeDataOpsReal::add()
      // Expected: v3 = v0 + v1 = (6.0)
      edge_ops->add(svindx[3], svindx[0], svindx[1]);
      double val_add = 6.0;
      if (!doubleDataSameAsValue(svindx[3], val_add, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #7: math::HierarchyEdgeDataOpsReal::add()\n"
         << "Expected: v3 = " << val_add << std::endl;
         edge_ops->printData(svindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0)
      edge_ops->setToScalar(svindx[0], 0.0);

      // Test #8: math::HierarchyEdgeDataOpsReal::subtract()
      // Expected: v1 = v3 - v0 = (6.0)
      edge_ops->subtract(svindx[1], svindx[3], svindx[0]);
      double val_sub = 6.0;
      if (!doubleDataSameAsValue(svindx[1], val_sub, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #8: math::HierarchyEdgeDataOpsReal::subtract()\n"
         << "Expected: v1 = " << val_sub << std::endl;
         edge_ops->printData(svindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchyEdgeDataOpsReal::addScalar()
      // Expected:  v1 = v1 + (0.0) = (6.0)
      edge_ops->addScalar(svindx[1], svindx[1], 0.0);
      double val_addScalar = 6.0;
      if (!doubleDataSameAsValue(svindx[1], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9a: math::HierarchyEdgeDataOpsReal::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << std::endl;
         edge_ops->printData(svindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchyEdgeDataOpsReal::addScalar()
      // Expected:  v2 = v2 + (0.0) = (1.0)
      edge_ops->addScalar(svindx[2], svindx[2], 0.0);
      val_addScalar = 1.0;
      if (!doubleDataSameAsValue(svindx[2], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9b: math::HierarchyEdgeDataOpsReal::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         edge_ops->printData(svindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchyEdgeDataOpsReal::addScalar()
      // Expected:  v2 = v2 + (3.0) = (4.0)
      edge_ops->addScalar(svindx[2], svindx[2], 3.0);
      val_addScalar = 4.0;
      if (!doubleDataSameAsValue(svindx[2], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9c: math::HierarchyEdgeDataOpsReal::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         edge_ops->printData(svindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5)
      edge_ops->setToScalar(svindx[3], 0.5);

      // Test #10: math::HierarchyEdgeDataOpsReal::multiply()
      // Expected: v1 = v3 * v1 = (3.0)
      edge_ops->multiply(svindx[1], svindx[3], svindx[1]);
      double val_mult = 3.0;
      if (!doubleDataSameAsValue(svindx[1], val_mult, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #10 math::HierarchyEdgeDataOpsReal::multiply()\n"
         << "Expected: v1 = " << val_mult << std::endl;
         edge_ops->printData(svindx[1], tbox::plog);
      }

      // Test #11: math::HierarchyEdgeDataOpsReal::divide()
      // Expected: v0 = v2 / v1 =  1.33333333333333
      edge_ops->divide(svindx[0], svindx[2], svindx[1]);
      double val_div = 1.333333333333;
      if (!doubleDataSameAsValue(svindx[0], val_div, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #11 math::HierarchyEdgeDataOpsReal::divide()\n"
         << "Expected: v0 = " << val_div << std::endl;
         edge_ops->printData(svindx[0], tbox::plog);
      }

      // Test #12: math::HierarchyEdgeDataOpsReal::reciprocal()
      // Expected: v1 = 1 / v1 = (0.333333333)
      edge_ops->reciprocal(svindx[1], svindx[1]);
      double val_rec = 0.33333333333333;
      if (!doubleDataSameAsValue(svindx[1], val_rec, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #12 math::HierarchyEdgeDataOpsReal::reciprocal()\n"
         << "Expected: v1 = " << val_rec << std::endl;
         edge_ops->printData(svindx[1], tbox::plog);
      }

      // Test #13: math::HierarchyEdgeDataOpsReal::abs()
      // Expected: v3 = abs(v2) = 4.0
      edge_ops->abs(svindx[3], svindx[2]);
      double val_abs = 4.0;
      if (!doubleDataSameAsValue(svindx[3], val_abs, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #13 math::HierarchyEdgeDataOpsReal::abs()\n"
         << "Expected: v3 = " << val_abs << std::endl;
         edge_ops->printData(svindx[3], tbox::plog);
      }

      // Test #14: Place some bogus values on coarse level
      boost::shared_ptr<pdat::EdgeData<double> > cdata;

      // set values
      boost::shared_ptr<hier::PatchLevel> level_zero(
        hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level_zero->begin());
           ip != level_zero->end(); ++ip) {
         patch = *ip;
         cdata = boost::dynamic_pointer_cast<pdat::EdgeData<double>,
                                             hier::PatchData>(patch->getPatchData(svindx[2]));
         hier::Index index0(2, 2);
         hier::Index index1(5, 3);
         if (patch->getBox().contains(index0)) {
            (*cdata)(pdat::EdgeIndex(index0, pdat::EdgeIndex::Y,
                        pdat::EdgeIndex::Lower), 0) = 100.0;
         }
         if (patch->getBox().contains(index1)) {
            (*cdata)(pdat::EdgeIndex(index1, pdat::EdgeIndex::Y,
                        pdat::EdgeIndex::Upper), 0) = -1000.0;
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level_zero->begin());
           ipp != level_zero->end(); ++ipp) {
         patch = *ipp;
         cdata = boost::dynamic_pointer_cast<pdat::EdgeData<double>,
                                             hier::PatchData>(patch->getPatchData(svindx[2]));
         pdat::EdgeIndex index0(hier::Index(2,
                                   2), pdat::EdgeIndex::Y,
                                pdat::EdgeIndex::Lower);
         pdat::EdgeIndex index1(hier::Index(5,
                                   3), pdat::EdgeIndex::Y,
                                pdat::EdgeIndex::Upper);

         // check X axis data
         pdat::EdgeIterator cend(cdata->getBox(), pdat::EdgeIndex::X, false);
         for (pdat::EdgeIterator c(cdata->getBox(), pdat::EdgeIndex::X, true);
              c != cend && bogus_value_test_passed;
              ++c) {
            pdat::EdgeIndex edge_index = *c;

            if (!tbox::MathUtilities<double>::equalEps((*cdata)(edge_index),
                   4.0)) {
               bogus_value_test_passed = false;
            }
         }

         // check Y axis data
         pdat::EdgeIterator ccend(cdata->getBox(), pdat::EdgeIndex::Y, false);
         for (pdat::EdgeIterator cc(cdata->getBox(), pdat::EdgeIndex::Y, true);
              cc != ccend && bogus_value_test_passed;
              ++cc) {
            pdat::EdgeIndex edge_index = *cc;

            if (edge_index == index0) {
               if (!tbox::MathUtilities<double>::equalEps((*cdata)(edge_index),
                      100.0)) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (edge_index == index1) {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            edge_index), -1000.0)) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            edge_index), 4.0)) {
                     bogus_value_test_passed = false;
                  }
               }
            }
         }
      }
      if (!bogus_value_test_passed) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #14:  Place some bogus values on coarse level"
         << std::endl;
         edge_ops->printData(svindx[2], tbox::plog);
      }

      // Test #15: math::HierarchyEdgeDataOpsReal::L1Norm() - w/o control weights
      // Expected: bogus_l1_norm = 1984.00
      double bogus_l1_norm = edge_ops->L1Norm(svindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, 1984.00)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #15: math::HierarchyEdgeDataOpsReal::L1Norm()"
         << " - w/o control weights\n"
         << "Expected value = 1984.00, Computed value = "
         << std::setprecision(12) << bogus_l1_norm << std::endl;
      }

      // Test #16: math::HierarchyEdgeDataOpsReal::L1Norm() - w/control weights
      // Expected: correct_l1_norm = 4.0
      double correct_l1_norm = edge_ops->L1Norm(svindx[2], swgt_id);
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 4.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #16: math::HierarchyEdgeDataOpsReal::L1Norm()"
         << " - w/control weights\n"
         << "Expected value = 4.0, Computed value = "
         << correct_l1_norm << std::endl;
      }

      // Test #17: math::HierarchyEdgeDataOpsReal::L2Norm()
      // Expected: l2_norm = 4.0
      double l2_norm = edge_ops->L2Norm(svindx[2], swgt_id);
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 4.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #17: math::HierarchyEdgeDataOpsReal::L2Norm()\n"
         << "Expected value = 4.0, Computed value = "
         << l2_norm << std::endl;
      }

      // Test #18: math::HierarchyEdgeDataOpsReal::maxNorm() - w/o control weights
      // Expected: bogus_max_norm = 1000.0
      double bogus_max_norm = edge_ops->maxNorm(svindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #18: math::HierarchyEdgeDataOpsReal::maxNorm()"
         << " - w/o control weights\n"
         << "Expected value = 1000.0, Computed value = "
         << bogus_max_norm << std::endl;
      }

      // Test #19: math::HierarchyEdgeDataOpsReal::maxNorm() - w/control weights
      // Expected: max_norm = 4.0
      double max_norm = edge_ops->maxNorm(svindx[2], swgt_id);
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 4.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #19: math::HierarchyEdgeDataOpsReal::maxNorm()"
         << " - w/control weights\n"
         << "Expected value = 4.0, Computed value = "
         << max_norm << std::endl;
      }

      // Reset data and test sums, axpy's
      edge_ops->setToScalar(svindx[0], 1.0);
      edge_ops->setToScalar(svindx[1], 2.5);
      edge_ops->setToScalar(svindx[2], 7.0);

      // Test #20: math::HierarchyEdgeDataOpsReal::linearSum()
      // Expected: v3 = 5.0
      edge_ops->linearSum(svindx[3], 2.0, svindx[1], 0.0, svindx[0]);
      double val_linearSum = 5.0;
      if (!doubleDataSameAsValue(svindx[3], val_linearSum, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #20: math::HierarchyEdgeDataOpsReal::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << std::endl;
         edge_ops->printData(svindx[3], tbox::plog);
      }

      // Test #21: math::HierarchyEdgeDataOpsReal::axmy()
      // Expected: v3 = 6.5
      edge_ops->axmy(svindx[3], 3.0, svindx[1], svindx[0]);
      double val_axmy = 6.5;
      if (!doubleDataSameAsValue(svindx[3], val_axmy, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #21: math::HierarchyEdgeDataOpsReal::axmy()\n"
         << "Expected: v3 = " << val_axmy << std::endl;
         edge_ops->printData(svindx[3], tbox::plog);
      }

      // Test #22a: math::HierarchyEdgeDataOpsReal::dot() - (ind2) * (ind1)
      // Expected: cdot = 17.5
      double cdot = edge_ops->dot(svindx[2], svindx[1], swgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 17.5)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #22a: math::HierarchyEdgeDataOpsReal::dot() - (ind2) * (ind1)\n"
         << "Expected Value = 17.5, Computed Value = "
         << cdot << std::endl;
      }

      // Test #22b: math::HierarchyEdgeDataOpsReal::dot() - (ind2) * (ind1)
      // Expected: cdot = 17.5
      cdot = edge_ops->dot(svindx[1], svindx[2], swgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 17.5)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #22b: math::HierarchyEdgeDataOpsReal::dot() - (ind2) * (ind1)\n"
         << "Expected Value = 17.5, Computed Value = "
         << cdot << std::endl;
      }

      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(swgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(svindx[iv]);
         }
      }

      for (iv = 0; iv < NVARS; iv++) {
         fvar[iv].reset();
      }
      swgt.reset();

      geometry.reset();
      hierarchy.reset();
      edge_ops.reset();
      swgt_ops.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  edge hiertest" << std::endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}

/*
 * Returns true if all the data in the hierarchy is equal to the specified
 * value.  Returns false otherwise.
 */
static bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   boost::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   bool test_passed = true;

   int ln;
   boost::shared_ptr<hier::Patch> patch;
   for (ln = 0; ln < 2; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         patch = *ip;
         boost::shared_ptr<pdat::EdgeData<double> > cvdata(
               patch->getPatchData(desc_id),
               boost::detail::dynamic_cast_tag());

         pdat::EdgeIterator cend(cvdata->getBox(), 1, false);
         for (pdat::EdgeIterator c(cvdata->getBox(), 1, true);
              c != cend && test_passed; ++c) {
            pdat::EdgeIndex edge_index = *c;
            if (!tbox::MathUtilities<double>::equalEps((*cvdata)(edge_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
}
