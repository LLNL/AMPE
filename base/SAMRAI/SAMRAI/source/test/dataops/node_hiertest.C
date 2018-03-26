/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program to test node-centered patch data ops
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
using namespace std;

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/SAMRAIManager.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyNodeDataOpsComplex.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyNodeDataOpsReal.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeVariable.h"
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

   const std::string log_fn = std::string("node_hiertest.")
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

      int coarse0_array[2][2] = { { 0, 0 }, { 9, 2 } };
      hier::Box coarse0((hier::Index(dim2d, coarse0_array[0])),
                        (hier::Index(dim2d, coarse0_array[1])),
                        hier::BlockId(0));

      int coarse1_array[2][2] = { { 0, 3 }, { 9, 4 } };
      hier::Box coarse1((hier::Index(dim2d, coarse1_array[0])),
                        (hier::Index(dim2d, coarse1_array[1])),
                        hier::BlockId(0));

      int fine0_array[2][2] = { { 4, 4 }, { 7, 7 } };
      hier::Box fine0((hier::Index(dim2d, fine0_array[0])),
                      (hier::Index(dim2d, fine0_array[1])),
                      hier::BlockId(0));

      int fine1_array[2][2] = { { 8, 4 }, { 13, 7 } };
      hier::Box fine1((hier::Index(dim2d, fine1_array[0])),
                      (hier::Index(dim2d, fine1_array[1])),
                      hier::BlockId(0));

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
      boost::shared_ptr<pdat::NodeVariable<double> > nvar[NVARS];
      int nvindx[NVARS];
      nvar[0].reset(new pdat::NodeVariable<double>(dim2d, "nvar0", 1));
      nvindx[0] = variable_db->registerVariableAndContext(
            nvar[0], dummy, no_ghosts);
      nvar[1].reset(new pdat::NodeVariable<double>(dim2d, "nvar1", 1));
      nvindx[1] = variable_db->registerVariableAndContext(
            nvar[1], dummy, no_ghosts);
      nvar[2].reset(new pdat::NodeVariable<double>(dim2d, "nvar2", 1));
      nvindx[2] = variable_db->registerVariableAndContext(
            nvar[2], dummy, no_ghosts);
      nvar[3].reset(new pdat::NodeVariable<double>(dim2d, "nvar3", 1));
      nvindx[3] = variable_db->registerVariableAndContext(
            nvar[3], dummy, no_ghosts);

      boost::shared_ptr<pdat::NodeVariable<double> > nwgt(
         new pdat::NodeVariable<double>(dim2d, "nwgt", 1));
      int nwgt_id = variable_db->registerVariableAndContext(
            nwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         level->allocatePatchData(nwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            level->allocatePatchData(nvindx[iv]);
         }
      }

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > node_ops(
         new math::HierarchyNodeDataOpsReal<double>(
            hierarchy,
            0,
            1));
      TBOX_ASSERT(node_ops);

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > nwgt_ops(
         new math::HierarchyNodeDataOpsReal<double>(
            hierarchy,
            0,
            1));

      boost::shared_ptr<hier::Patch> patch;

      // Initialize control volume data for node-centered components
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
            const double node_vol = dx[0] * dx[1];
            boost::shared_ptr<pdat::NodeData<double> > data(
               patch->getPatchData(nwgt_id),
               boost::detail::dynamic_cast_tag());
            data->fillAll(node_vol);
            pdat::NodeIndex ni(dim);
            int plo0 = patch->getBox().lower(0);
            int phi0 = patch->getBox().upper(0);
            int plo1 = patch->getBox().lower(1);
            int phi1 = patch->getBox().upper(1);
            int ic;

            if (ln == 0) {
               data->fillAll(0.0, (coarse_fine * patch->getBox()));

               if (patch->getLocalId() == 0) {
                  //bottom boundaries
                  for (ic = plo0; ic < phi0; ic++) {
                     ni = pdat::NodeIndex(hier::Index(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) *= 0.5;
                  }
                  //left and right boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index(plo0,
                              ic), pdat::NodeIndex::UpperLeft);
                     (*data)(ni) *= 0.5;
                     ni = pdat::NodeIndex(hier::Index(phi0,
                              ic), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 0.5;
                  }
                  // corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index(plo0,
                                plo1), pdat::NodeIndex::LowerLeft)) *= 0.25;
                  (*data)(pdat::NodeIndex(hier::Index(phi0,
                                plo1), pdat::NodeIndex::LowerRight)) *= 0.25;
               } else {
                  //top and bottom boundaries
                  for (ic = plo0; ic < phi0; ic++) {
                     ni = pdat::NodeIndex(hier::Index(ic,
                              phi1), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 0.5;
                     ni = pdat::NodeIndex(hier::Index(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) = 0.0;
                  }
                  //left and right boundaries
                  for (ic = plo1; ic < phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index(plo0,
                              ic), pdat::NodeIndex::UpperLeft);
                     (*data)(ni) *= 0.5;
                     ni = pdat::NodeIndex(hier::Index(phi0,
                              ic), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 0.5;
                  }
                  // corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index(plo0,
                                plo1), pdat::NodeIndex::LowerLeft)) = 0.0;
                  (*data)(pdat::NodeIndex(hier::Index(plo0,
                                phi1), pdat::NodeIndex::UpperLeft)) *= 0.25;
                  (*data)(pdat::NodeIndex(hier::Index(phi0,
                                plo1), pdat::NodeIndex::LowerRight)) = 0.0;
                  (*data)(pdat::NodeIndex(hier::Index(phi0,
                                phi1), pdat::NodeIndex::UpperRight)) *= 0.25;
               }
            } else {
               if (patch->getLocalId() == 0) {
                  // top and bottom coarse-fine boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     ni = pdat::NodeIndex(hier::Index(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) *= 1.5;
                     ni = pdat::NodeIndex(hier::Index(ic,
                              phi1), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 1.5;
                  }
                  //left coarse-fine boundaries
                  for (ic = plo1; ic < phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index(plo0,
                              ic), pdat::NodeIndex::UpperLeft);
                     (*data)(ni) *= 1.5;
                  }
                  // coarse-fine corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index(plo0,
                                plo1), pdat::NodeIndex::LowerLeft)) *= 2.25;
                  (*data)(pdat::NodeIndex(hier::Index(plo0,
                                phi1), pdat::NodeIndex::UpperLeft)) *= 2.25;
               } else {
                  // top and bottom coarse-fine boundaries
                  for (ic = plo0; ic < phi0; ic++) {
                     ni = pdat::NodeIndex(hier::Index(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) *= 1.5;
                     ni = pdat::NodeIndex(hier::Index(ic,
                              phi1), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 1.5;
                  }
                  //right coarse-fine boundaries
                  for (ic = plo1; ic < phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index(phi0,
                              ic), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 1.5;
                  }
                  // coarse-fine corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index(phi0,
                                plo1), pdat::NodeIndex::LowerRight)) *= 2.25;
                  (*data)(pdat::NodeIndex(hier::Index(phi0,
                                phi1), pdat::NodeIndex::UpperRight)) *= 2.25;
                  //shared left boundaries
                  for (ic = plo1; ic <= phi1 + 1; ic++) {
                     ni = pdat::NodeIndex(hier::Index(plo0,
                              ic), pdat::NodeIndex::LowerLeft);
                     (*data)(ni) = 0;
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
 *   boost::shared_ptr< pdat::NodeData<double> > cvdata = patch->getPatchData(cwgt_id);
 *
 *   pdat::NodeIterator cend(cvdata->getBox(), false);
 *   for (pdat::NodeIterator c(cvdata->getBox(), true); c != cend && vol_test_passed; ++c) {
 *   pdat::NodeIndex node_index = *c;
 *
 *   if (ln == 0) {
 *   if ((coarse_fine * patch->getBox()).contains(node_index)) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(node_index),0.0) ) {
 *   vol_test_passed = false;
 *   }
 *   } else {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(node_index),0.01) ) {
 *   vol_test_passed = false;
 *   }
 *   }
 *   }
 *
 *   if (ln == 1) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*cvdata)(node_index),0.0025) ) {
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
/*
 * tbox::plog << "node control volume data" << std::endl;
 * nwgt_ops->printData(nwgt_id, tbox::plog);
 */

      // Test #1b: math::HierarchyNodeDataOpsReal::sumControlVolumes()
      // Expected: norm = 0.5
      double norm = node_ops->sumControlVolumes(nvindx[0], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm, 0.5)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #1b: math::HierarchyNodeDataOpsReal::sumControlVolumes()\n"
         << "Expected value = 0.5 , Computed value = "
         << norm << std::endl;
      }

      // Test #2: math::HierarchyNodeDataOpsReal::numberOfEntries()
      // Expected: num_data_points = 121
      int num_data_points = node_ops->numberOfEntries(nvindx[0]);
      if (num_data_points != 121) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #2: math::HierarchyNodeDataOpsReal::numberOfEntries()\n"
         << "Expected value = 121 , Computed value = "
         << num_data_points << std::endl;
      }

      // Test #3a: math::HierarchyNodeDataOpsReal::setToScalar()
      // Expected: v0 = 2.0
      double val0 = double(2.0);
      node_ops->setToScalar(nvindx[0], val0);
      if (!doubleDataSameAsValue(nvindx[0], val0, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #3a: math::HierarchyNodeDataOpsReal::setToScalar()\n"
         << "Expected: v0 = " << val0 << std::endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchyNodeDataOpsReal::setToScalar()
      // Expected: v1 = (4.0)
      node_ops->setToScalar(nvindx[1], 4.0);
      double val1 = 4.0;
      if (!doubleDataSameAsValue(nvindx[1], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #3b: math::HierarchyNodeDataOpsReal::setToScalar()\n"
         << "Expected: v1 = " << val1 << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #4: math::HierarchyNodeDataOpsReal::copyData()
      // Expected: v2 = v1 = (4.0)
      node_ops->copyData(nvindx[2], nvindx[1]);
      if (!doubleDataSameAsValue(nvindx[2], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #4: math::HierarchyNodeDataOpsReal::setToScalar()\n"
         << "Expected: v2 = " << val1 << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #5: math::HierarchyNodeDataOpsReal::swapData()
      // Expected: v0 = (4.0), v1 = (2.0)
      node_ops->swapData(nvindx[0], nvindx[1]);
      if (!doubleDataSameAsValue(nvindx[0], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #5a: math::HierarchyNodeDataOpsReal::setToScalar()\n"
         << "Expected: v0 = " << val1 << std::endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }
      if (!doubleDataSameAsValue(nvindx[1], val0, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #5b: math::HierarchyNodeDataOpsReal::setToScalar()\n"
         << "Expected: v1 = " << val0 << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #6: math::HierarchyNodeDataOpsReal::scale()
      // Expected: v2 = 0.25 * v2 = (1.0)
      node_ops->scale(nvindx[2], 0.25, nvindx[2]);
      double val_scale = 1.0;
      if (!doubleDataSameAsValue(nvindx[2], val_scale, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #6: math::HierarchyNodeDataOpsReal::scale()\n"
         << "Expected: v2 = " << val_scale << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #7: math::HierarchyNodeDataOpsReal::add()
      // Expected: v3 = v0 + v1 = (6.0)
      node_ops->add(nvindx[3], nvindx[0], nvindx[1]);
      double val_add = 6.0;
      if (!doubleDataSameAsValue(nvindx[3], val_add, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #7: math::HierarchyNodeDataOpsReal::add()\n"
         << "Expected: v3 = " << val_add << std::endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0)
      node_ops->setToScalar(nvindx[0], 0.0);

      // Test #8: math::HierarchyNodeDataOpsReal::subtract()
      // Expected: v1 = v3 - v0 = (6.0)
      node_ops->subtract(nvindx[1], nvindx[3], nvindx[0]);
      double val_sub = 6.0;
      if (!doubleDataSameAsValue(nvindx[1], val_sub, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #8: math::HierarchyNodeDataOpsReal::subtract()\n"
         << "Expected: v1 = " << val_sub << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchyNodeDataOpsReal::addScalar()
      // Expected: v1 = v1 + (0.0) = (6.0)
      node_ops->addScalar(nvindx[1], nvindx[1], 0.0);
      double val_addScalar = 6.0;
      if (!doubleDataSameAsValue(nvindx[1], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9a: math::HierarchyNodeDataOpsReal::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchyNodeDataOpsReal::addScalar()
      // Expected: v2 = v2 + (0.0) = (1.0)
      node_ops->addScalar(nvindx[2], nvindx[2], 0.0);
      val_addScalar = 1.0;
      if (!doubleDataSameAsValue(nvindx[2], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9b: math::HierarchyNodeDataOpsReal::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchyNodeDataOpsReal::addScalar()
      // Expected:  v2 = v2 + (3.0) = (4.0)
      node_ops->addScalar(nvindx[2], nvindx[2], 3.0);
      val_addScalar = 4.0;
      if (!doubleDataSameAsValue(nvindx[2], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9c: math::HierarchyNodeDataOpsReal::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5)
      node_ops->setToScalar(nvindx[3], 0.5);

      // Test #10: math::HierarchyNodeDataOpsReal::multiply()
      // Expected:  v1 = v3 * v1 = (3.0)
      node_ops->multiply(nvindx[1], nvindx[3], nvindx[1]);
      double val_mult = 3.0;
      if (!doubleDataSameAsValue(nvindx[1], val_mult, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #10: math::HierarchyNodeDataOpsReal::multiply()\n"
         << "Expected: v1 = " << val_mult << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #11: math::HierarchyNodeDataOpsReal::divide()
      // Expected:  v0 = v2 / v1 = 1.3333333333333
      node_ops->divide(nvindx[0], nvindx[2], nvindx[1]);
      double val_div = 1.333333333333333;
      if (!doubleDataSameAsValue(nvindx[0], val_div, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #11: math::HierarchyNodeDataOpsReal::divide()\n"
         << "Expected: v0 = " << val_div << std::endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #12: math::HierarchyNodeDataOpsReal::reciprocal()
      // Expected:  v1 = 1 / v1 = (0.333333333)
      node_ops->reciprocal(nvindx[1], nvindx[1]);
      double val_rec = 0.33333333333333;
      if (!doubleDataSameAsValue(nvindx[1], val_rec, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #12: math::HierarchyNodeDataOpsReal::reciprocal()\n"
         << "Expected: v1 = " << val_rec << std::endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #13: math::HierarchyNodeDataOpsReal::abs()
      // Expected:  v3 = abs(v2) = 4.0
      node_ops->abs(nvindx[3], nvindx[2]);
      double val_abs = 4.0;
      if (!doubleDataSameAsValue(nvindx[3], val_abs, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #13: math::HierarchyNodeDataOpsReal::abs()\n"
         << "Expected: v3 = " << val_abs << std::endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #14: Place some bogus values on coarse level
      boost::shared_ptr<pdat::NodeData<double> > ndata;

      // set values
      boost::shared_ptr<hier::PatchLevel> level_zero(
         hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level_zero->begin());
           ip != level_zero->end(); ++ip) {
         patch = *ip;
         ndata = boost::dynamic_pointer_cast<pdat::NodeData<double>,
                                             hier::PatchData>(patch->getPatchData(nvindx[2]));
         hier::Index index0(2, 2);
         hier::Index index1(5, 3);
         if (patch->getBox().contains(index0)) {
            (*ndata)(pdat::NodeIndex(index0,
                        pdat::NodeIndex::LowerLeft), 0) = 100.0;
         }
         if (patch->getBox().contains(index1)) {
            (*ndata)(pdat::NodeIndex(index1,
                        pdat::NodeIndex::UpperRight), 0) = -1000.0;
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level_zero->begin());
           ipp != level_zero->end(); ++ipp) {
         patch = *ipp;
         ndata = boost::dynamic_pointer_cast<pdat::NodeData<double>,
                                             hier::PatchData>(patch->getPatchData(nvindx[2]));
         pdat::NodeIndex index0(hier::Index(2, 2), pdat::NodeIndex::LowerLeft);
         pdat::NodeIndex index1(hier::Index(5, 3), pdat::NodeIndex::UpperRight);

         pdat::NodeIterator cend(ndata->getBox(), false);
         for (pdat::NodeIterator c(ndata->getBox(), true);
              c != cend && bogus_value_test_passed;
              ++c) {
            pdat::NodeIndex node_index = *c;

            if (node_index == index0) {
               if (!tbox::MathUtilities<double>::equalEps((*ndata)(node_index),
                      100.0)) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (node_index == index1) {
                  if (!tbox::MathUtilities<double>::equalEps((*ndata)(
                            node_index), -1000.0)) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<double>::equalEps((*ndata)(
                            node_index), 4.0)) {
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
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #15: math::HierarchyNodeDataOpsReal::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1640.00
      double bogus_l1_norm = node_ops->L1Norm(nvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, 1640.00)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #15: math::HierarchyNodeDataOpsReal::L1Norm()"
         << " - w/o control weight\n"
         << "Expected value = 1640.00, Computed value = "
         << std::setprecision(12) << bogus_l1_norm << std::endl;
      }

      // Test #16: math::HierarchyNodeDataOpsReal::L1Norm() - w/control weight
      // Expected:  correct_l1_norm = 2.0
      double correct_l1_norm = node_ops->L1Norm(nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 2.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #16: math::HierarchyNodeDataOpsReal::L1Norm()"
         << " - w/control weight\n"
         << "Expected value = 2.0, Computed value = "
         << correct_l1_norm << std::endl;
      }

      // Test #17: math::HierarchyNodeDataOpsReal::L2Norm()
      // Expected:  l2_norm = 2.8284271
      double l2_norm = node_ops->L2Norm(nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 2.82842712475)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #17: math::HierarchyNodeDataOpsReal::L2Norm()\n"
         << "Expected value = 2.82842712475, Computed value = "
         << l2_norm << std::endl;
      }

      // Test #18: math::HierarchyNodeDataOpsReal::maxNorm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.0
      double bogus_max_norm = node_ops->maxNorm(nvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #18: math::HierarchyNodeDataOpsReal::maxNorm()"
         << " - w/o control weight\n"
         << "Expected value = 1000.0, Computed value = "
         << bogus_max_norm << std::endl;
      }

      // Test #19: math::HierarchyNodeDataOpsReal::maxNorm() - w/control weight
      // Expected:  max_norm = 4.0
      double max_norm = node_ops->maxNorm(nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 4.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #19: math::HierarchyNodeDataOpsReal::maxNorm()"
         << " - w/control weight\n"
         << "Expected value = 4.0, Computed value = "
         << max_norm << std::endl;
      }

      // Reset data and test sums, axpy's
      node_ops->setToScalar(nvindx[0], 1.0);
      node_ops->setToScalar(nvindx[1], 2.5);
      node_ops->setToScalar(nvindx[2], 7.0);

      // Test #20: math::HierarchyNodeDataOpsReal::linearSum()
      // Expected:  v3 = 5.0
      double val_linearSum = 5.0;
      node_ops->linearSum(nvindx[3], 2.0, nvindx[1], 0.00, nvindx[0]);
      if (!doubleDataSameAsValue(nvindx[3], val_linearSum, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #20: math::HierarchyNodeDataOpsReal::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << std::endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #21: math::HierarchyNodeDataOpsReal::axmy()
      // Expected:  v3 = 6.5
      node_ops->axmy(nvindx[3], 3.0, nvindx[1], nvindx[0]);
      double val_axmy = 6.5;
      if (!doubleDataSameAsValue(nvindx[3], val_axmy, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #21: math::HierarchyNodeDataOpsReal::axmy()\n"
         << "Expected: v3 = " << val_axmy << std::endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #22a: math::HierarchyNodeDataOpsReal::dot() - (ind2) * (ind1)
      // Expected:  cdot = 8.75
      double cdot = node_ops->dot(nvindx[2], nvindx[1], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #22a: math::HierarchyNodeDataOpsReal::dot() - (ind2) * (ind1)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << std::endl;
      }

      // Test #22a: math::HierarchyNodeDataOpsReal::dot() - (ind1) * (ind2)
      // Expected:  cdot = 8.75
      cdot = node_ops->dot(nvindx[1], nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #22a: math::HierarchyNodeDataOpsReal::dot() - (ind1) * (ind2)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << std::endl;
      }

      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(nwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(nvindx[iv]);
         }
      }

      for (iv = 0; iv < NVARS; iv++) {
         nvar[iv].reset();
      }
      nwgt.reset();

      geometry.reset();
      hierarchy.reset();
      node_ops.reset();
      nwgt_ops.reset();

      tbox::pout << "\nPASSED:  node hiertest" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
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
         boost::shared_ptr<pdat::NodeData<double> > nvdata(
               patch->getPatchData(desc_id),
               boost::detail::dynamic_cast_tag());

         pdat::NodeIterator cend(nvdata->getBox(), false);
         for (pdat::NodeIterator c(nvdata->getBox(), true);
              c != cend && test_passed; ++c) {
            pdat::NodeIndex node_index = *c;
            if (!tbox::MathUtilities<double>::equalEps((*nvdata)(node_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
}
