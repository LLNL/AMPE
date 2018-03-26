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
#include "SAMRAI/hier/BoxArray.h"
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
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

#include <boost/shared_ptr.hpp>

using namespace SAMRAI;

/* Helper function declarations */
bool
doubleDataSameAsValue(
   int desc_id,
   double value,
   boost::shared_ptr<hier::PatchHierarchy> hierarchy);

#define NVARS 4

int main(
   int argc,
   char* argv[]) {

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
// tbox::PIO::logOnlyNodeZero("node_hierops.log");
   tbox::PIO::logAllNodes("node_hierops.log");

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      int ln, iv;

      /*
       * Make a simple 2-level hierarchy.
       */
      double lo[2] = { 0.0, 0.0 };
      double hi[2] = { 1.0, 0.5 };

      hier::Box<2> coarse0(hier::Index<2>(0, 0), hier::Index<2>(9, 2));
      hier::Box<2> coarse1(hier::Index<2>(0, 3), hier::Index<2>(9, 4));
      hier::Box<2> fine0(hier::Index<2>(4, 4), hier::Index<2>(7, 7));
      hier::Box<2> fine1(hier::Index<2>(8, 4), hier::Index<2>(13, 7));
      hier::IntVector<2> ratio(2);

      hier::BoxArray<2> coarse_domain(2);
      hier::BoxArray<2> fine_domain(2);
      coarse_domain(0) = coarse0;
      coarse_domain(1) = coarse1;
      fine_domain(0) = fine0;
      fine_domain(1) = fine1;

      boost::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            lo,
            hi,
            coarse_domain));

      boost::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      // Note: For these simple tests we allow at most 2 processors.
      tbox::SAMRAI_MPI mpi(SAMRAIManager::getSAMRAICommWorld());
      const int nproc = mpi.getSize();
      TBOX_ASSERT(nproc < 3);

      const int n_coarse_boxes = coarse_domain.getNumberOfBoxes();
      const int n_fine_boxes = fine_domain.getNumberOfBoxes();
      hier::ProcessorMapping mapping0(n_coarse_boxes);
      hier::ProcessorMapping mapping1(n_fine_boxes);

      int ib;
      for (ib = 0; ib < n_coarse_boxes; ib++) {
         if (nproc > 1) {
            mapping0.setProcessorAssignment(ib, ib);
         } else {
            mapping0.setProcessorAssignment(ib, 0);
         }
      }

      for (ib = 0; ib < n_fine_boxes; ib++) {
         if (nproc > 1) {
            mapping1.setProcessorAssignment(ib, ib);
         } else {
            mapping1.setProcessorAssignment(ib, 0);
         }
      }

      hierarchy->makeNewPatchLevel(0, hier::IntVector<2>(
            1), coarse_domain, mapping0);
      hierarchy->makeNewPatchLevel(1, ratio, fine_domain, mapping1);

      /*
       * Create some variables, a context, and register them with
       * the variable database.
       */
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      boost::shared_ptr<hier::VariableContext> dummy(
         variable_db->getContext("dummy"));
      const hier::IntVector<2> no_ghosts(0);

      boost::shared_ptr<pdat::NodeVariable<double> > nvar[NVARS];
      int nvindx[NVARS];
      nvar[0].reset(new pdat::NodeVariable<double>("nvar0", 1));
      nvindx[0] = variable_db->registerVariableAndContext(
            nvar[0], dummy, no_ghosts);
      nvar[1].reset(new pdat::NodeVariable<double>("nvar1", 1));
      nvindx[1] = variable_db->registerVariableAndContext(
            nvar[1], dummy, no_ghosts);
      nvar[2].reset(new pdat::NodeVariable<double>("nvar2", 1));
      nvindx[2] = variable_db->registerVariableAndContext(
            nvar[2], dummy, no_ghosts);
      nvar[3].reset(new pdat::NodeVariable<double>("nvar3", 1));
      nvindx[3] = variable_db->registerVariableAndContext(
            nvar[3], dummy, no_ghosts);

      boost::shared_ptr<pdat::NodeVariable<double> > nwgt(
         new pdat::NodeVariable<double>("nwgt", 1));
      int nwgt_id = variable_db->registerVariableAndContext(
            nwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(nwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->allocatePatchData(nvindx[iv]);
         }
      }

      /*
       * Create instances of hierarchy operations to apply certain
       * mathematical operators.  e.g. scale(), axpy(), min(), etc.
       */
      int coarsest = 0;
      int finest = 1;
      boost::shared_ptr<math::HierarchyDataOpsReal<double> > node_ops(
         new math::HierarchyNodeDataOpsReal<double>(
            hierarchy,
            coarsest,
            finest));
      TBOX_ASSERT(node_ops);

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > nwgt_ops(
         new math::HierarchyNodeDataOpsReal<double>(
            hierarchy,
            coarsest,
            finest));

      boost::shared_ptr<hier::Patch> patch;
      boost::shared_ptr<geom::CartesianPatchGeometry> pgeom;

      // Initialize control volume data for node-centered components
      hier::Box<2> coarse_fine = fine0 + fine1;
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            boost::shared_ptr<pdat::NodeData<double> > data;
            patch = level->getPatch(ip());
            pgeom = patch->getPatchGeometry();
            const double* dx = pgeom->getDx();
            const double node_vol = dx[0] * dx[1];
            data = patch->getPatchData(nwgt_id);
            data->fillAll(node_vol);
            pdat::NodeIndex ni;
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
                     ni = pdat::NodeIndex(hier::Index<2>(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) *= 0.5;
                  }
                  //left and right boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(plo0,
                              ic), pdat::NodeIndex::UpperLeft);
                     (*data)(ni) *= 0.5;
                     ni = pdat::NodeIndex(hier::Index<2>(phi0,
                              ic), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 0.5;
                  }
                  // corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index<2>(plo0,
                                plo1), pdat::NodeIndex::LowerLeft)) *= 0.25;
                  (*data)(pdat::NodeIndex(hier::Index<2>(phi0,
                                plo1), pdat::NodeIndex::LowerRight)) *= 0.25;
               } else {
                  //top and bottom boundaries
                  for (ic = plo0; ic < phi0; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(ic,
                              phi1), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 0.5;
                     ni = pdat::NodeIndex(hier::Index<2>(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) = 0.0;
                  }
                  //left and right boundaries
                  for (ic = plo1; ic < phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(plo0,
                              ic), pdat::NodeIndex::UpperLeft);
                     (*data)(ni) *= 0.5;
                     ni = pdat::NodeIndex(hier::Index<2>(phi0,
                              ic), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 0.5;
                  }
                  // corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index<2>(plo0,
                                plo1), pdat::NodeIndex::LowerLeft)) = 0.0;
                  (*data)(pdat::NodeIndex(hier::Index<2>(plo0,
                                phi1), pdat::NodeIndex::UpperLeft)) *= 0.25;
                  (*data)(pdat::NodeIndex(hier::Index<2>(phi0,
                                plo1), pdat::NodeIndex::LowerRight)) = 0.0;
                  (*data)(pdat::NodeIndex(hier::Index<2>(phi0,
                                phi1), pdat::NodeIndex::UpperRight)) *= 0.25;
               }
            } else {
               if (patch->getLocalId() == 0) {
                  // top and bottom coarse-fine boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) *= 1.5;
                     ni = pdat::NodeIndex(hier::Index<2>(ic,
                              phi1), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 1.5;
                  }
                  //left coarse-fine boundaries
                  for (ic = plo1; ic < phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(plo0,
                              ic), pdat::NodeIndex::UpperLeft);
                     (*data)(ni) *= 1.5;
                  }
                  // coarse-fine corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index<2>(plo0,
                                plo1), pdat::NodeIndex::LowerLeft)) *= 2.25;
                  (*data)(pdat::NodeIndex(hier::Index<2>(plo0,
                                phi1), pdat::NodeIndex::UpperLeft)) *= 2.25;
               } else {
                  // top and bottom coarse-fine boundaries
                  for (ic = plo0; ic < phi0; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(ic,
                              plo1), pdat::NodeIndex::LowerRight);
                     (*data)(ni) *= 1.5;
                     ni = pdat::NodeIndex(hier::Index<2>(ic,
                              phi1), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 1.5;
                  }
                  //right coarse-fine boundaries
                  for (ic = plo1; ic < phi1; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(phi0,
                              ic), pdat::NodeIndex::UpperRight);
                     (*data)(ni) *= 1.5;
                  }
                  // coarse-fine corner boundaries
                  (*data)(pdat::NodeIndex(hier::Index<2>(phi0,
                                plo1), pdat::NodeIndex::LowerRight)) *= 2.25;
                  (*data)(pdat::NodeIndex(hier::Index<2>(phi0,
                                phi1), pdat::NodeIndex::UpperRight)) *= 2.25;
                  //shared left boundaries
                  for (ic = plo1; ic <= phi1 + 1; ic++) {
                     ni = pdat::NodeIndex(hier::Index<2>(plo0,
                              ic), pdat::NodeIndex::LowerLeft);
                     (*data)(ni) = 0;
                  }
               }
            }
         }
      }

      // Test #1b: HierarchyNodeDataOpsReal2::sumControlVolumes()
      // Expected: norm = 0.5
      double norm = node_ops->sumControlVolumes(nvindx[0], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm, 0.5)) {
         tbox::perr
         << "FAILED: - Test #1b: HierarchyNodeDataOpsReal2::sumControlVolumes()\n"
         << "Expected value = 0.5 , Computed value = "
         << norm << endl;
      }

      // Test #2: HierarchyNodeDataOpsReal2::numberOfEntries()
      // Expected: num_data_points = 121
      int num_data_points = node_ops->numberOfEntries(nvindx[0]);
      if (num_data_points != 121) {
         tbox::perr
         << "FAILED: - Test #2: HierarchyNodeDataOpsReal2::numberOfEntries()\n"
         << "Expected value = 121 , Computed value = "
         << num_data_points << endl;
      }

      // Test #3a: HierarchyNodeDataOpsReal2::setToScalar()
      // Expected: v0 = 2.0
      double val0 = double(2.0);
      node_ops->setToScalar(nvindx[0], val0);
      if (!doubleDataSameAsValue(nvindx[0], val0, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #3a: HierarchyNodeDataOpsReal2::setToScalar()\n"
         << "Expected: v0 = " << val0 << endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #3b: HierarchyNodeDataOpsReal2::setToScalar()
      // Expected: v1 = (4.0)
      node_ops->setToScalar(nvindx[1], 4.0);
      double val1 = 4.0;
      if (!doubleDataSameAsValue(nvindx[1], val1, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #3b: HierarchyNodeDataOpsReal2::setToScalar()\n"
         << "Expected: v1 = " << val1 << endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #4: HierarchyNodeDataOpsReal2::copyData()
      // Expected: v2 = v1 = (4.0)
      node_ops->copyData(nvindx[2], nvindx[1]);
      if (!doubleDataSameAsValue(nvindx[2], val1, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #4: HierarchyNodeDataOpsReal2::setToScalar()\n"
         << "Expected: v2 = " << val1 << endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #5: HierarchyNodeDataOpsReal2::swapData()
      // Expected: v0 = (4.0), v1 = (2.0)
      node_ops->swapData(nvindx[0], nvindx[1]);
      if (!doubleDataSameAsValue(nvindx[0], val1, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #5a: HierarchyNodeDataOpsReal2::setToScalar()\n"
         << "Expected: v0 = " << val1 << endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }
      if (!doubleDataSameAsValue(nvindx[1], val0, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #5b: HierarchyNodeDataOpsReal2::setToScalar()\n"
         << "Expected: v1 = " << val0 << endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #6: HierarchyNodeDataOpsReal2::scale()
      // Expected: v2 = 0.25 * v2 = (1.0)
      node_ops->scale(nvindx[2], 0.25, nvindx[2]);
      double val_scale = 1.0;
      if (!doubleDataSameAsValue(nvindx[2], val_scale, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #6: HierarchyNodeDataOpsReal2::scale()\n"
         << "Expected: v2 = " << val_scale << endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #7: HierarchyNodeDataOpsReal2::add()
      // Expected: v3 = v0 + v1 = (6.0)
      node_ops->add(nvindx[3], nvindx[0], nvindx[1]);
      double val_add = 6.0;
      if (!doubleDataSameAsValue(nvindx[3], val_add, hierarchy)) {
         tbox::perr << "FAILED: - Test #7: HierarchyNodeDataOpsReal2::add()\n"
                    << "Expected: v3 = " << val_add << endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0)
      node_ops->setToScalar(nvindx[0], 0.0);

      // Test #8: HierarchyNodeDataOpsReal2::subtract()
      // Expected: v1 = v3 - v0 = (6.0)
      node_ops->subtract(nvindx[1], nvindx[3], nvindx[0]);
      double val_sub = 6.0;
      if (!doubleDataSameAsValue(nvindx[1], val_sub, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #8: HierarchyNodeDataOpsReal2::subtract()\n"
         << "Expected: v1 = " << val_sub << endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9a: HierarchyNodeDataOpsReal2::addScalar()
      // Expected: v1 = v1 + (0.0) = (6.0)
      node_ops->addScalar(nvindx[1], nvindx[1], 0.0);
      double val_addScalar = 6.0;
      if (!doubleDataSameAsValue(nvindx[1], val_addScalar, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #9a: HierarchyNodeDataOpsReal2::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #9b: HierarchyNodeDataOpsReal2::addScalar()
      // Expected: v2 = v2 + (0.0) = (1.0)
      node_ops->addScalar(nvindx[2], nvindx[2], 0.0);
      val_addScalar = 1.0;
      if (!doubleDataSameAsValue(nvindx[2], val_addScalar, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #9b: HierarchyNodeDataOpsReal2::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #9c: HierarchyNodeDataOpsReal2::addScalar()
      // Expected:  v2 = v2 + (3.0) = (4.0)
      node_ops->addScalar(nvindx[2], nvindx[2], 3.0);
      val_addScalar = 4.0;
      if (!doubleDataSameAsValue(nvindx[2], val_addScalar, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #9c: HierarchyNodeDataOpsReal2::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5)
      node_ops->setToScalar(nvindx[3], 0.5);

      // Test #10: HierarchyNodeDataOpsReal2::multiply()
      // Expected:  v1 = v3 * v1 = (3.0)
      node_ops->multiply(nvindx[1], nvindx[3], nvindx[1]);
      double val_mult = 3.0;
      if (!doubleDataSameAsValue(nvindx[1], val_mult, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #10: HierarchyNodeDataOpsReal2::multiply()\n"
         << "Expected: v1 = " << val_mult << endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #11: HierarchyNodeDataOpsReal2::divide()
      // Expected:  v0 = v2 / v1 = 1.3333333333333
      node_ops->divide(nvindx[0], nvindx[2], nvindx[1]);
      double val_div = 1.333333333333333;
      if (!doubleDataSameAsValue(nvindx[0], val_div, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #11: HierarchyNodeDataOpsReal2::divide()\n"
         << "Expected: v0 = " << val_div << endl;
         node_ops->printData(nvindx[0], tbox::plog);
      }

      // Test #12: HierarchyNodeDataOpsReal2::reciprocal()
      // Expected:  v1 = 1 / v1 = (0.333333333)
      node_ops->reciprocal(nvindx[1], nvindx[1]);
      double val_rec = 0.33333333333333;
      if (!doubleDataSameAsValue(nvindx[1], val_rec, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #12: HierarchyNodeDataOpsReal2::reciprocal()\n"
         << "Expected: v1 = " << val_rec << endl;
         node_ops->printData(nvindx[1], tbox::plog);
      }

      // Test #13: HierarchyNodeDataOpsReal2::abs()
      // Expected:  v3 = abs(v2) = 4.0
      node_ops->abs(nvindx[3], nvindx[2]);
      double val_abs = 4.0;
      if (!doubleDataSameAsValue(nvindx[3], val_abs, hierarchy)) {
         tbox::perr << "FAILED: - Test #13: HierarchyNodeDataOpsReal2::abs()\n"
                    << "Expected: v3 = " << val_abs << endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #14: Place some bogus values on coarse level
      boost::shared_ptr<pdat::NodeData<double> > ndata;

      // set values
      boost::shared_ptr<hier::PatchLevel> level_zero(
         hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level_zero->begin());
           ip != level_zero->end(); ++ip) {
         patch = level_zero->getPatch(ip());
         ndata = patch->getPatchData(nvindx[2]);
         hier::Index<2> index0(2, 2);
         hier::Index<2> index1(5, 3);
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
         patch = level_zero->getPatch(ipp());
         ndata = patch->getPatchData(nvindx[2]);
         pdat::NodeIndex index0(hier::Index<2>(2,
                                               2), pdat::NodeIndex::LowerLeft);
         pdat::NodeIndex index1(hier::Index<2>(5,
                                               3), pdat::NodeIndex::UpperRight);

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
         tbox::perr
         << "FAILED: - Test #14:  Place some bogus values on coarse level"
         << endl;
         node_ops->printData(nvindx[2], tbox::plog);
      }

      // Test #15: HierarchyNodeDataOpsReal2::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1640.00
      double bogus_l1_norm = node_ops->L1Norm(nvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, 1640.00)) {
         tbox::perr
         << "FAILED: - Test #15: HierarchyNodeDataOpsReal2::L1Norm()"
         << " - w/o control weight\n"
         << "Expected value = 1640.00, Computed value = "
         << setprecision(12) << bogus_l1_norm << endl;
      }

      // Test #16: HierarchyNodeDataOpsReal2::L1Norm() - w/control weight
      // Expected:  correct_l1_norm = 2.0
      double correct_l1_norm = node_ops->L1Norm(nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 2.0)) {
         tbox::perr
         << "FAILED: - Test #16: HierarchyNodeDataOpsReal2::L1Norm()"
         << " - w/control weight\n"
         << "Expected value = 2.0, Computed value = "
         << correct_l1_norm << endl;
      }

      // Test #17: HierarchyNodeDataOpsReal2::L2Norm()
      // Expected:  l2_norm = 2.8284271
      double l2_norm = node_ops->L2Norm(nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 2.82842712475)) {
         tbox::perr
         << "FAILED: - Test #17: HierarchyNodeDataOpsReal2::L2Norm()\n"
         << "Expected value = 2.82842712475, Computed value = "
         << l2_norm << endl;
      }

      // Test #18: HierarchyNodeDataOpsReal2::maxNorm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.0
      double bogus_max_norm = node_ops->maxNorm(nvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.0)) {
         tbox::perr
         << "FAILED: - Test #18: HierarchyNodeDataOpsReal2::maxNorm()"
         << " - w/o control weight\n"
         << "Expected value = 1000.0, Computed value = "
         << bogus_max_norm << endl;
      }

      // Test #19: HierarchyNodeDataOpsReal2::maxNorm() - w/control weight
      // Expected:  max_norm = 4.0
      double max_norm = node_ops->maxNorm(nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 4.0)) {
         tbox::perr
         << "FAILED: - Test #19: HierarchyNodeDataOpsReal2::maxNorm()"
         << " - w/control weight\n"
         << "Expected value = 4.0, Computed value = "
         << max_norm << endl;
      }

      // Reset data and test sums, axpy's
      node_ops->setToScalar(nvindx[0], 1.0);
      node_ops->setToScalar(nvindx[1], 2.5);
      node_ops->setToScalar(nvindx[2], 7.0);

      // Test #20: HierarchyNodeDataOpsReal2::linearSum()
      // Expected:  v3 = 5.0
      double val_linearSum = 5.0;
      node_ops->linearSum(nvindx[3], 2.0, nvindx[1], 0.00, nvindx[0]);
      if (!doubleDataSameAsValue(nvindx[3], val_linearSum, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #20: HierarchyNodeDataOpsReal2::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #21: HierarchyNodeDataOpsReal2::axmy()
      // Expected:  v3 = 6.5
      node_ops->axmy(nvindx[3], 3.0, nvindx[1], nvindx[0]);
      double val_axmy = 6.5;
      if (!doubleDataSameAsValue(nvindx[3], val_axmy, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #21: HierarchyNodeDataOpsReal2::axmy()\n"
         << "Expected: v3 = " << val_axmy << endl;
         node_ops->printData(nvindx[3], tbox::plog);
      }

      // Test #22a: HierarchyNodeDataOpsReal2::dot() - (ind2) * (ind1)
      // Expected:  cdot = 8.75
      double cdot = node_ops->dot(nvindx[2], nvindx[1], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         tbox::perr
         << "FAILED: - Test #22a: HierarchyNodeDataOpsReal2::dot() - (ind2) * (ind1)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << endl;
      }

      // Test #22a: HierarchyNodeDataOpsReal2::dot() - (ind1) * (ind2)
      // Expected:  cdot = 8.75
      cdot = node_ops->dot(nvindx[1], nvindx[2], nwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         tbox::perr
         << "FAILED: - Test #22a: HierarchyNodeDataOpsReal2::dot() - (ind1) * (ind2)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << endl;
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
bool
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
         patch = level->getPatch(ip());
         boost::shared_ptr<pdat::NodeData<double> > nvdata(
            patch->getPatchData(desc_id));

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
