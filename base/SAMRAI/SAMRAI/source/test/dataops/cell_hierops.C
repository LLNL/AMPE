/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program to test cell-centered patch data ops
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
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyCellDataOpsComplex.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
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

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

// tbox::PIO::logOnlyNodeZero("cell_hierops.log");
      tbox::PIO::logAllNodes("cell_hierops.log");

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
      coarse_domain[0] = coarse0;
      coarse_domain[1] = coarse1;
      fine_domain[0] = fine0;
      fine_domain[1] = fine1;

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

      boost::shared_ptr<pdat::CellVariable<double> > cvar[NVARS];
      int cvindx[NVARS];
      cvar[0].reset(new pdat::CellVariable<double>("cvar0", 1));
      cvindx[0] = variable_db->registerVariableAndContext(
            cvar[0], dummy, no_ghosts);
      cvar[1].reset(new pdat::CellVariable<double>("cvar1", 1));
      cvindx[1] = variable_db->registerVariableAndContext(
            cvar[1], dummy, no_ghosts);
      cvar[2].reset(new pdat::CellVariable<double>("cvar2", 1));
      cvindx[2] = variable_db->registerVariableAndContext(
            cvar[2], dummy, no_ghosts);
      cvar[3].reset(new pdat::CellVariable<double>("cvar3", 1));
      cvindx[3] = variable_db->registerVariableAndContext(
            cvar[3], dummy, no_ghosts);

      boost::shared_ptr<pdat::CellVariable<double> > cwgt(
         new pdat::CellVariable<double>("cwgt", 1);
      int cwgt_id = variable_db->registerVariableAndContext(
            cwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(cwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->allocatePatchData(cvindx[iv]);
         }
      }

      /*
       * Create instances of hierarchy operations to apply certain
       * mathematical operators.  e.g. scale(), axpy(), min(), etc.
       */
      int coarsest = 0;
      int finest = 1;
      boost::shared_ptr<math::HierarchyDataOpsReal<double> > cell_ops(
         new math::HierarchyCellDataOpsReal<double>(
            hierarchy,
            coarsest,
            finest));
      TBOX_ASSERT(cell_ops);

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > cwgt_ops(
         new math::HierarchyCellDataOpsReal<double>(
            hierarchy,
            coarsest,
            finest));

      boost::shared_ptr<hier::Patch> patch;
      boost::shared_ptr<geom::CartesianPatchGeometry> pgeom;

      // Initialize control volume data for cell-centered components
      hier::Box<2> coarse_fine = fine0 + fine1;
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = level->getPatch(ip());
            pgeom = patch->getPatchGeometry();
            const double* dx = pgeom->getDx();
            const double cell_vol = dx[0] * dx[1];
            boost::shared_ptr<pdat::CellData<double> > cvdata(
               patch->getPatchData(cwgt_id));
            cvdata->fillAll(cell_vol);
            if (ln == 0) cvdata->fillAll(0.0, (coarse_fine * patch->getBox()));
         }
      }

      /*
       * Apply various operations to the hierarchy data. Test the
       * result to assure its accuracy.
       */
      // Test #1a: Check control volume data set properly
      // Expected: cwgt = 0.01 on coarse (except where finer patch exists) and
      // 0.0025 on fine level
      bool vol_test_passed = true;
      for (ln = 0; ln < 2; ln++) {

         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));

         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            patch = level->getPatch(ip());
            boost::shared_ptr<pdat::CellData<double> > cvdata(
               patch->getPatchData(cwgt_id));

            pdat::CellIterator cend(cvdata->getBox(), false);
            for (pdat::CellIterator c(cvdata->getBox(), true);
                 c != cend && vol_test_passed;
                 ++c) {
               pdat::CellIndex cell_index = *c;

               if (ln == 0) {
                  if ((coarse_fine * patch->getBox()).contains(cell_index)) {
                     if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                               cell_index), 0.0)) {
                        vol_test_passed = false;
                     }
                  } else {
                     if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                               cell_index), 0.01)) {
                        vol_test_passed = false;
                     }
                  }
               }

               if (ln == 1) {
                  if (!tbox::MathUtilities<double>::equalEps((*cvdata)(
                            cell_index), 0.0025)) {
                     vol_test_passed = false;
                  }
               }
            }
         }
      }
      if (!vol_test_passed) {
         tbox::perr
         << "FAILED: - Test #1a: Check control volume data set properly"
         << endl;
         cwgt_ops->printData(cwgt_id, tbox::pout);
      }

      // Test #1b: HierarchyCellDataOpsReal2::sumControlVolumes()
      // Expected: norm = 0.5
      double norm = cell_ops->sumControlVolumes(cvindx[0], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm, 0.5)) {
         tbox::perr
         << "FAILED: - Test #1b: HierarchyCellDataOpsReal2::sumControlVolumes()\n"
         << "Expected value = 0.5 , Computed value = "
         << norm << endl;
      }

      // Test #2: HierarchyCellDataOpsReal2::numberOfEntries()
      // Expected: num_data_points = 90
      int num_data_points = cell_ops->numberOfEntries(cvindx[0]);
      if (num_data_points != 90) {
         tbox::perr
         << "FAILED: - Test #2: HierarchyCellDataOpsReal2::numberOfEntries()\n"
         << "Expected value = 90 , Computed value = "
         << num_data_points << endl;
      }

      // Test #3a: HierarchyCellDataOpsReal2::setToScalar()
      // Expected: v0 = 2.0
      double val0 = 2.0;
      cell_ops->setToScalar(cvindx[0], val0);
      if (!doubleDataSameAsValue(cvindx[0], val0, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #3a: HierarchyCellDataOpsReal2::setToScalar()\n"
         << "Expected: v0 = " << val0 << endl;
         cell_ops->printData(cvindx[0], tbox::pout);
      }

      // Test #3b: HierarchyCellDataOpsReal2::setToScalar()
      // Expected: v1 = (4.0)
      cell_ops->setToScalar(cvindx[1], 4.0);
      double val1 = 4.0;
      if (!doubleDataSameAsValue(cvindx[1], val1, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #3b: HierarchyCellDataOpsReal2::setToScalar()\n"
         << "Expected: v1 = " << val1 << endl;
         cell_ops->printData(cvindx[1], tbox::pout);
      }

      // Test #4: HierarchyCellDataOpsReal2::copyData()
      // Expected: v2 = v1 = (4.0)
      cell_ops->copyData(cvindx[2], cvindx[1]);
      if (!doubleDataSameAsValue(cvindx[2], val1, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #4: HierarchyCellDataOpsReal2::copyData()\n"
         << "Expected: v2 = " << val1 << endl;
         cell_ops->printData(cvindx[2], tbox::pout);
      }

      // Test #5: HierarchyCellDataOpsReal2::swapData()
      // Expected: v0 = (4.0), v1 = (2.0)
      cell_ops->swapData(cvindx[0], cvindx[1]);
      if (!doubleDataSameAsValue(cvindx[0], val1, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #5a: HierarchyCellDataOpsReal2::swapData()\n"
         << "Expected: v0 = " << val1 << endl;
         cell_ops->printData(cvindx[0], tbox::pout);
      }
      if (!doubleDataSameAsValue(cvindx[1], val0, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #5b: HierarchyCellDataOpsReal2::swapData()\n"
         << "Expected: v1 = " << val0 << endl;
         cell_ops->printData(cvindx[1], tbox::pout);
      }

      // Test #6: HierarchyCellDataOpsReal2::scale()
      // Expected: v2 = 0.25 * v2 = (1.0)
      cell_ops->scale(cvindx[2], 0.25, cvindx[2]);
      double val_scale = 1.0;
      if (!doubleDataSameAsValue(cvindx[2], val_scale, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #6: HierarchyCellDataOpsReal2::scale()\n"
         << "Expected: v2 = " << val_scale << endl;
         cell_ops->printData(cvindx[2], tbox::pout);
      }

      // Test #7: HierarchyCellDataOpsReal2::add()
      // Expected: v3 = v0 + v1 = (6.0)
      cell_ops->add(cvindx[3], cvindx[0], cvindx[1]);
      double val_add = 6.0;
      if (!doubleDataSameAsValue(cvindx[3], val_add, hierarchy)) {
         tbox::perr << "FAILED: - Test #7: HierarchyCellDataOpsReal2::add()\n"
                    << "Expected: v3 = " << val_add << endl;
         cell_ops->printData(cvindx[3], tbox::pout);
      }

      // Reset v0: v0 = (0.0)
      cell_ops->setToScalar(cvindx[0], 0.0);

      // Test #8: HierarchyCellDataOpsReal2::subtract()
      // Expected: v1 = v3 - v0 = (6.0)
      cell_ops->subtract(cvindx[1], cvindx[3], cvindx[0]);
      double val_sub = 6.0;
      if (!doubleDataSameAsValue(cvindx[1], val_sub, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #8: HierarchyCellDataOpsReal2::subtract()\n"
         << "Expected: v1 = " << val_sub << endl;
         cell_ops->printData(cvindx[1], tbox::pout);
      }

      // Test #9a: HierarchyCellDataOpsReal2::addScalar()
      // Expected: v1 = v1 + (0.0) = (6.0)
      cell_ops->addScalar(cvindx[1], cvindx[1], 0.0);
      double val_addScalar = 6.0;
      if (!doubleDataSameAsValue(cvindx[1], val_addScalar, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #9a: HierarchyCellDataOpsReal2::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << endl;
         cell_ops->printData(cvindx[1], tbox::pout);
      }

      // Test #9b: HierarchyCellDataOpsReal2::addScalar()
      // Expected: v2 = v2 + (0.0) = (1.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], 0.0);
      val_addScalar = 1.0;
      if (!doubleDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #9b: HierarchyCellDataOpsReal2::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << endl;
         cell_ops->printData(cvindx[2], tbox::pout);
      }

      // Test #9c: HierarchyCellDataOpsReal2::addScalar()
      // Expected: v2 = v2 + (3.0) = (4.0)
      cell_ops->addScalar(cvindx[2], cvindx[2], 3.0);
      val_addScalar = 4.0;
      if (!doubleDataSameAsValue(cvindx[2], val_addScalar, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #9c: HierarchyCellDataOpsReal2::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << endl;
         cell_ops->printData(cvindx[2], tbox::pout);
      }

      // Reset v3:  v3 = (0.5)
      cell_ops->setToScalar(cvindx[3], 0.5);

      // Test #10: HierarchyCellDataOpsReal2::multiply()
      // Expected: v1 = v3 * v1 = (3.0)
      cell_ops->multiply(cvindx[1], cvindx[3], cvindx[1]);
      double val_mult = 3.0;
      if (!doubleDataSameAsValue(cvindx[1], val_mult, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #10: HierarchyCellDataOpsReal2::multiply()\n"
         << "Expected: v1 = " << val_mult << endl;
         cell_ops->printData(cvindx[1], tbox::pout);
      }

      // Test #11: HierarchyCellDataOpsReal2::divide()
      // Expected: v0 = v2 / v1 = 1.3333333333
      cell_ops->divide(cvindx[0], cvindx[2], cvindx[1]);
      double val_div = 1.33333333333;
      if (!doubleDataSameAsValue(cvindx[0], val_div, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #11: HierarchyCellDataOpsReal2::divide()\n"
         << "Expected: v0 = " << val_div << endl;
         cell_ops->printData(cvindx[0], tbox::pout);
      }

      // Test #12: HierarchyCellDataOpsReal2::reciprocal()
      // Expected:  v1 = 1 / v1 = (0.333333333)
      cell_ops->reciprocal(cvindx[1], cvindx[1]);
      double val_rec = 0.33333333333;
      if (!doubleDataSameAsValue(cvindx[1], val_rec, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #12: HierarchyCellDataOpsReal2::reciprocal()\n"
         << "Expected: v1 = " << val_rec << endl;
         cell_ops->printData(cvindx[1], tbox::pout);
      }

      // Test #13: HierarchyCellDataOpsReal2::abs()
      // Expected:  v3 = abs(v2) = 4.0
      cell_ops->abs(cvindx[3], cvindx[2]);
      double val_abs = 4.0;
      if (!doubleDataSameAsValue(cvindx[3], val_abs, hierarchy)) {
         tbox::perr << "FAILED: - Test #13: HierarchyCellDataOpsReal2::abs()\n"
                    << "Expected: v3 = " << val_abs << endl;
         cell_ops->printData(cvindx[3], tbox::pout);
      }

      // Test #14: Place some bogus values on coarse level
      boost::shared_ptr<pdat::CellData<double> > cdata;

      // set values
      boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         patch = level->getPatch(ip());
         cdata = patch->getPatchData(cvindx[2]);
         hier::Index<2> index0(2, 2);
         hier::Index<2> index1(5, 3);
         if (patch->getBox().contains(index0)) {
            (*cdata)(pdat::CellIndex(index0), 0) = 100.0;
         }
         if (patch->getBox().contains(index1)) {
            (*cdata)(pdat::CellIndex(index1), 0) = -1000.0;
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level->begin());
           ipp != level->end(); ++ipp) {
         patch = level->getPatch(ipp());
         cdata = patch->getPatchData(cvindx[2]);
         hier::Index<2> index0(2, 2);
         hier::Index<2> index1(5, 3);

         pdat::CellIterator cend(cdata->getBox(), false);
         for (pdat::CellIterator c(cdata->getBox(), true);
              c != cend && bogus_value_test_passed;
              ++c) {
            pdat::CellIndex cell_index = *c;

            if (cell_index == index0) {
               if (!tbox::MathUtilities<double>::equalEps((*cdata)(cell_index),
                      100.0)) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (cell_index == index1) {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            cell_index), -1000.0)) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<double>::equalEps((*cdata)(
                            cell_index), 4.0)) {
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
         cell_ops->printData(cvindx[2], tbox::pout);
      }

      // Test #15: HierarchyCellDataOpsReal2::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 1452
      double bogus_l1_norm = cell_ops->L1Norm(cvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, 1452)) {
         tbox::perr
         << "FAILED: - Test #15: HierarchyCellDataOpsReal2::L1Norm()"
         << " - w/o control weight\n"
         << "Expected value = 1452, Computed value = "
         << setprecision(12) << bogus_l1_norm << endl;
      }

      // Test #16: HierarchyCellDataOpsReal2::L1Norm() - w/control weight
      // Expected:  correct_l1_norm = 2.0
      double correct_l1_norm = cell_ops->L1Norm(cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 2.0)) {
         tbox::perr
         << "FAILED: - Test #16: HierarchyCellDataOpsReal2::L1Norm()"
         << " - w/control weight\n"
         << "Expected value = 2.0, Computed value = "
         << correct_l1_norm << endl;
      }

      // Test #17: HierarchyCellDataOpsReal2::L2Norm()
      // Expected:  l2_norm = 2.82842712475
      double l2_norm = cell_ops->L2Norm(cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 2.82842712475)) {
         tbox::perr
         << "FAILED: - Test #17: HierarchyCellDataOpsReal2::L2Norm()\n"
         << "Expected value = 2.82842712475, Computed value = "
         << l2_norm << endl;
      }

      // Test #18: HierarchyCellDataOpsReal2::L2Norm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.0
      double bogus_max_norm = cell_ops->maxNorm(cvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.0)) {
         tbox::perr
         << "FAILED: - Test #18: HierarchyCellDataOpsReal2::L2Norm()"
         << " - w/o control weight\n"
         << "Expected value = 1000.0, Computed value = "
         << bogus_max_norm << endl;
      }

      // Test #19: HierarchyCellDataOpsReal2::L2Norm() - w/control weight
      // Expected:  max_norm = 4.0
      double max_norm = cell_ops->maxNorm(cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 4.0)) {
         tbox::perr
         << "FAILED: - Test #19: HierarchyCellDataOpsReal2::L2Norm()"
         << " - w/control weight\n"
         << "Expected value = 4.0, Computed value = "
         << max_norm << endl;
      }

      // Reset data and test sums, axpy's
      cell_ops->setToScalar(cvindx[0], 1.00);
      cell_ops->setToScalar(cvindx[1], 2.5);
      cell_ops->setToScalar(cvindx[2], 7.0);

      // Test #20: HierarchyCellDataOpsReal2::linearSum()
      // Expected:  v3 = 5.0
      cell_ops->linearSum(cvindx[3], 2.0, cvindx[1], 0.00, cvindx[0]);
      double val_linearSum = 5.0;
      if (!doubleDataSameAsValue(cvindx[3], val_linearSum, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #20: HierarchyCellDataOpsReal2::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << endl;
         cell_ops->printData(cvindx[3], tbox::pout);
      }

      // Test #21: HierarchyCellDataOpsReal2::axmy()
      // Expected:  v3 = 6.5
      cell_ops->axmy(cvindx[3], 3.0, cvindx[1], cvindx[0]);
      double val_axmy = 6.5;
      if (!doubleDataSameAsValue(cvindx[3], val_axmy, hierarchy)) {
         tbox::perr
         << "FAILED: - Test #21: HierarchyCellDataOpsReal2::axmy()\n"
         << "Expected: v3 = " << val_axmy << endl;
         cell_ops->printData(cvindx[3], tbox::pout);
      }

      // Test #22a: HierarchyCellDataOpsReal2::dot() - (ind2) * (ind1)
      // Expected:  cdot = 8.75
      double cdot = cell_ops->dot(cvindx[2], cvindx[1], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         tbox::perr
         << "FAILED: - Test #22a: HierarchyCellDataOpsReal2::dot() - (ind2) * (ind1)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << endl;
      }

      // Test #22b: HierarchyCellDataOpsReal2::dot() - (ind1) * (ind2)
      // Expected:  cdot = 8.75
      cdot = cell_ops->dot(cvindx[1], cvindx[2], cwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(cdot, 8.75)) {
         tbox::perr
         << "FAILED: - Test #22b: HierarchyCellDataOpsReal2::dot() - (ind1) * (ind2)\n"
         << "Expected Value = 8.75, Computed Value = "
         << cdot << endl;
      }

      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(cwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(cvindx[iv]);
         }
      }

      for (iv = 0; iv < NVARS; iv++) {
         cvar[iv].reset();
      }
      cwgt.reset();

      geometry.reset();
      hierarchy.reset();
      cell_ops.reset();
      cwgt_ops.reset();

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
         boost::shared_ptr<pdat::CellData<double> > cvdata(
            patch->getPatchData(desc_id));

         pdat::CellIterator cend(cvdata->getBox(), false);
         for (pdat::CellIterator c(cvdata->getBox(), true);
              c != cend && test_passed; ++c) {
            pdat::CellIndex cell_index = *c;
            if (!tbox::MathUtilities<double>::equalEps((*cvdata)(cell_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
}
