/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program to test face-centered complex patch data ops
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
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/math/HierarchyDataOpsComplex.h"
#include "SAMRAI/math/HierarchyFaceDataOpsComplex.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/math/HierarchyFaceDataOpsReal.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/pdat/FaceVariable.h"
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
complexDataSameAsValue(
   int desc_id,
   dcomplex value,
   boost::shared_ptr<hier::PatchHierarchy> hierarchy);
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

   if (argc < 2) {
      TBOX_ERROR("Usage: " << argv[0] << " [dimension]");
   }

   const int d = atoi(argv[1]);
   TBOX_ASSERT(d > 0);
   TBOX_ASSERT(d <= SAMRAI_MAXIMUM_DIMENSION_VALUE);
   const tbox::Dimension dim(d);

   if (dim != 2) {
      TBOX_ERROR("This test code is completed only for 2D!!!");
   }

   const std::string log_fn = std::string("face_cplxtest.")
      + tbox::Utilities::intToString(dim, 1) + "d.log";
   tbox::PIO::logAllNodes(log_fn);

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      const tbox::Dimension dim2(2);

      int ln, iv;

      // Make a dummy hierarchy domain
      double lo[2] = { 0.0, 0.0 };
      double hi[2] = { 1.0, 0.5 };

      hier::Box coarse0(hier::Index(0, 0), hier::Index(9, 2));
      hier::Box coarse1(hier::Index(0, 3), hier::Index(9, 4));
      hier::Box fine0(hier::Index(4, 4), hier::Index(7, 7));
      hier::Box fine1(hier::Index(8, 4), hier::Index(13, 7));
      hier::IntVector ratio(dim2, 2);

      hier::BoxContainer coarse_domain(dim);
      hier::BoxContainer fine_boxes(dim);
      coarse_domain.appendItem(coarse0);
      coarse_domain.appendItem(coarse1);
      fine_boxes.appendItem(fine0);
      fine_boxes.appendItem(fine1);

      boost::shared_ptr<geom::CartesianGridGeometry> geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            lo,
            hi,
            coarse_domain));

      boost::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy("PatchHierarchy", geometry));

      // Note: For these simple tests we allow at most 2 processors.
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      const int nproc = mpi.getSize();
      TBOX_ASSERT(nproc < 3);

      const int n_coarse_boxes = coarse_domain.getNumberOfBoxes();
      const int n_fine_boxes = fine_boxes.getNumberOfBoxes();

      hier::BoxLevel layer0(hier::IntVector(dim, 1), geometry);
      hier::BoxLevel layer1(ratio, geometry);

      hier::BoxContainer::iterator coarse_itr(coarse_domain);
      for (int ib = 0; ib < n_coarse_boxes; ib++, ++coarse_itr) {
         if (nproc > 1) {
            if (ib == layer0.getRank()) {
               layer0.addBox(hier::Box(*coarse_itr, ib,
                     layer0.getRank()));
            }
         } else {
            layer0.addBox(hier::Box(*coarse_itr, ib, 0));
         }
      }

      hier::BoxContainer::iterator fine_itr(fine_boxes);
      for (int ib = 0; ib < n_fine_boxes; ib++, ++fine_itr) {
         if (nproc > 1) {
            if (ib == layer1.getRank()) {
               layer1.addBox(hier::Box(*fine_itr, ib,
                     layer1.getRank()));
            }
         } else {
            layer1.addBox(hier::Box(*fine_itr, ib, 0));
         }
      }

      hierarchy->makeNewPatchLevel(0, layer0);
      hierarchy->makeNewPatchLevel(1, layer1);

      // Create instance of hier::Variable database
      hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
      boost::shared_ptr<hier::VariableContext> dummy(
         variable_db->getContext("dummy"));
      const hier::IntVector no_ghosts(dim2, 0);

      // Make some dummy variables and data on the hierarchy
      boost::shared_ptr<pdat::FaceVariable<dcomplex> > fvar[NVARS];
      int fvindx[NVARS];
      fvar[0].reset(new pdat::FaceVariable<dcomplex>(dim, "fvar0", 1));
      fvindx[0] = variable_db->registerVariableAndContext(
            fvar[0], dummy, no_ghosts);
      fvar[1].reset(new pdat::FaceVariable<dcomplex>(dim, "fvar1", 1));
      fvindx[1] = variable_db->registerVariableAndContext(
            fvar[1], dummy, no_ghosts);
      fvar[2].reset(new pdat::FaceVariable<dcomplex>(dim, "fvar2", 1));
      fvindx[2] = variable_db->registerVariableAndContext(
            fvar[2], dummy, no_ghosts);
      fvar[3].reset(new pdat::FaceVariable<dcomplex>(dim, "fvar3", 1));
      fvindx[3] = variable_db->registerVariableAndContext(
            fvar[3], dummy, no_ghosts);

      boost::shared_ptr<pdat::FaceVariable<double> > fwgt(
         new pdat::FaceVariable<double>(dim, "fwgt", 1));
      int fwgt_id = variable_db->registerVariableAndContext(
            fwgt, dummy, no_ghosts);

      // allocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->allocatePatchData(fwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->allocatePatchData(fvindx[iv]);
         }
      }

      boost::shared_ptr<math::HierarchyDataOpsComplex> face_ops(
         new math::HierarchyFaceDataOpsComplex(
            hierarchy,
            0,
            1));
      TBOX_ASSERT(face_ops);

      boost::shared_ptr<math::HierarchyDataOpsReal<double> > fwgt_ops(
         new math::HierarchyFaceDataOpsReal<double>(
            hierarchy,
            0,
            1));

      boost::shared_ptr<hier::Patch> patch;
      boost::shared_ptr<geom::CartesianPatchGeometry> pgeom;

      // Initialize control volume data for face-centered components
      hier::Box coarse_fine = fine0 + fine1;
      coarse_fine.coarsen(ratio);
      for (ln = 0; ln < 2; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            boost::shared_ptr<pdat::FaceData<double> > data;
            patch = level->getPatch(ip());
            pgeom = patch->getPatchGeometry();
            const double* dx = pgeom->getDx();
            const double face_vol = dx[0] * dx[1];
            data = patch->getPatchData(fwgt_id);
            data->fillAll(face_vol);
            pdat::FaceIndex fi(dim);
            int plo0 = patch->getBox().lower(0);
            int phi0 = patch->getBox().upper(0);
            int plo1 = patch->getBox().lower(1);
            int phi1 = patch->getBox().upper(1);
            int ic;

            if (ln == 0) {
               data->fillAll(0.0, (coarse_fine * patch->getBox()));

               if (patch->getLocalId() == 0) {
                  //bottom face boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::FaceIndex(hier::Index(ic,
                              plo1), pdat::FaceIndex::Y, pdat::FaceIndex::Lower);
                     (*data)(fi) *= 0.5;
                  }
                  //left and right face boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::FaceIndex(hier::Index(plo0,
                              ic), pdat::FaceIndex::X, pdat::FaceIndex::Lower);
                     (*data)(fi) *= 0.5;
                     fi = pdat::FaceIndex(hier::Index(phi0,
                              ic), pdat::FaceIndex::X, pdat::FaceIndex::Upper);
                     (*data)(fi) *= 0.5;
                  }
               } else {
                  //top and bottom face boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::FaceIndex(hier::Index(ic,
                              plo1), pdat::FaceIndex::Y, pdat::FaceIndex::Lower);
                     (*data)(fi) = 0.0;
                     fi = pdat::FaceIndex(hier::Index(ic,
                              phi1), pdat::FaceIndex::Y, pdat::FaceIndex::Upper);
                     (*data)(fi) *= 0.5;
                  }
                  //left and right face boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::FaceIndex(hier::Index(plo0,
                              ic), pdat::FaceIndex::X, pdat::FaceIndex::Lower);
                     (*data)(fi) *= 0.5;
                     fi = pdat::FaceIndex(hier::Index(phi0,
                              ic), pdat::FaceIndex::X, pdat::FaceIndex::Upper);
                     (*data)(fi) *= 0.5;
                  }
               }
            } else {
               if (patch->getLocalId() == 0) {
                  // top and bottom coarse-fine face boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::FaceIndex(hier::Index(ic,
                              plo1), pdat::FaceIndex::Y, pdat::FaceIndex::Lower);
                     (*data)(fi) *= 1.5;
                     fi = pdat::FaceIndex(hier::Index(ic,
                              phi1), pdat::FaceIndex::Y, pdat::FaceIndex::Upper);
                     (*data)(fi) *= 1.5;
                  }
                  //left coarse-fine face boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::FaceIndex(hier::Index(plo0,
                              ic), pdat::FaceIndex::X, pdat::FaceIndex::Lower);
                     (*data)(fi) *= 1.5;
                  }
               } else {
                  // top and bottom coarse-fine face boundaries
                  for (ic = plo0; ic <= phi0; ic++) {
                     fi = pdat::FaceIndex(hier::Index(ic,
                              plo1), pdat::FaceIndex::Y, pdat::FaceIndex::Lower);
                     (*data)(fi) *= 1.5;
                     fi = pdat::FaceIndex(hier::Index(ic,
                              phi1), pdat::FaceIndex::Y, pdat::FaceIndex::Upper);
                     (*data)(fi) *= 1.5;
                  }
                  //left and right coarse-fine face boundaries
                  for (ic = plo1; ic <= phi1; ic++) {
                     fi = pdat::FaceIndex(hier::Index(plo0,
                              ic), pdat::FaceIndex::X, pdat::FaceIndex::Lower);
                     (*data)(fi) = 0.0;
                     fi = pdat::FaceIndex(hier::Index(phi0,
                              ic), pdat::FaceIndex::X, pdat::FaceIndex::Upper);
                     (*data)(fi) *= 1.5;
                  }
               }
            }
         }
      }

      // Test #1: Print out control volume data and compute its integral

      // Test #1a: Check control volume data set properly
      // Expected: cwgt =
/*   bool vol_test_passed = true;
 *   for (ln = 0; ln < 2; ln++) {
 *   for (hier::PatchLevel::iterator ip(hierarchy->getPatchLevel(ln)->begin()); ip != hierarchy->getPatchLevel(ln)->end(); ++ip) {
 *   patch = hierarchy->getPatchLevel(ln)->getPatch(ip());
 *   boost::shared_ptr< pdat::FaceData<double> > fvdata = patch->getPatchData(fwgt_id);
 *
 *   pdat::FaceIterator cend(fvdata->getBox(),1, false);
 *   for (pdat::FaceIterator c(fvdata->getBox(), 1, true); c != cend && vol_test_passed; ++c) {
 *   pdat::FaceIndex face_index = *c;
 *
 *   if (ln == 0) {
 *   if ((coarse_fine * patch->getBox()).contains(face_index)) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*fvdata)(face_index),0.0) ) {
 *   vol_test_passed = false;
 *   }
 *   } else {
 *   if ( !tbox::MathUtilities<double>::equalEps((*fvdata)(face_index),0.01) ) {
 *   vol_test_passed = false;
 *   }
 *   }
 *   }
 *
 *   if (ln == 1) {
 *   if ( !tbox::MathUtilities<double>::equalEps((*fvdata)(face_index),0.0025) ) {
 *   vol_test_passed = false;
 *   }
 *   }
 *   }
 *   }
 *   }
 *   if (!vol_test_passed) {
 *   num_failures++;
 *   tbox::perr << "FAILED: - Test #1a: Check control volume data set properly" << std::endl;
 *   fwgt_ops->printData(fwgt_id, tbox::plog);
 *   }
 */
      // Test #1b: math::HierarchyFaceDataOpsComplex::sumControlVolumes()
      // Expected: norm = 1.0
      double norm =
         face_ops->sumControlVolumes(fvindx[0], fwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(norm, 1.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #1b: math::HierarchyFaceDataOpsComplex::sumControlVolumes()\n"
         << "Expected value = 1.0 , Computed value = "
         << norm << std::endl;
      }

      // Test #2: math::HierarchyFaceDataOpsComplex::numberOfEntries()
      // Expected: num_data_points = 209
      int num_data_points = face_ops->numberOfEntries(fvindx[0]);
      if (num_data_points != 209) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #2: math::HierarchyFaceDataOpsComplex::numberOfEntries()\n"
         << "Expected value = 209 , Computed value = "
         << num_data_points << std::endl;
      }

      // Test #3a: math::HierarchyFaceDataOpsComplex::setToScalar()
      // Expected: v0 = (2.0,1.5)
      dcomplex val0 = dcomplex(2.0, 1.5);
      face_ops->setToScalar(fvindx[0], val0);
      if (!complexDataSameAsValue(fvindx[0], val0, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #3a: math::HierarchyFaceDataOpsComplex::setToScalar()\n"
         << "Expected: v0 = " << val0 << std::endl;
         face_ops->printData(fvindx[0], tbox::plog);
      }

      // Test #3b: math::HierarchyFaceDataOpsComplex::setToScalar()
      // Expected: v1 = (4.0,3.0)
      dcomplex val1(4.0, 3.0);
      face_ops->setToScalar(fvindx[1], val1);
      if (!complexDataSameAsValue(fvindx[1], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #3b: math::HierarchyFaceDataOpsComplex::setToScalar()\n"
         << "Expected: v1 = " << val1 << std::endl;
         face_ops->printData(fvindx[1], tbox::plog);
      }

      // Test #4: math::HierarchyFaceDataOpsComplex::copyData()
      // Expected: v2 = v1 = (4.0,3.0)
      face_ops->copyData(fvindx[2], fvindx[1]);
      if (!complexDataSameAsValue(fvindx[2], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #3b: math::HierarchyFaceDataOpsComplex::copyData()\n"
         << "Expected: v2 = v1 = " << val1 << std::endl;
         face_ops->printData(fvindx[2], tbox::plog);
      }

      // Test #5: math::HierarchyFaceDataOpsComplex::swapData()
      // Expected:  v0 = (4.0, 3.0),  v1 = (2.0,1.5)
      face_ops->swapData(fvindx[0], fvindx[1]);
      if (!complexDataSameAsValue(fvindx[0], val1, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #5a: math::HierarchyFaceDataOpsComplex::swapData()\n"
         << "Expected: v0 = " << val1 << std::endl;
         face_ops->printData(fvindx[0], tbox::plog);
      }
      if (!complexDataSameAsValue(fvindx[1], val0, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #5b: math::HierarchyFaceDataOpsComplex::swapData()\n"
         << "Expected: v1 = " << val0 << std::endl;
         face_ops->printData(fvindx[1], tbox::plog);
      }

      // Test #6: math::HierarchyFaceDataOpsComplex::scale()
      // Expected:  v2 = 0.25 * v2 = (1.0,0.75)
      face_ops->scale(fvindx[2], 0.25, fvindx[2]);
      dcomplex val_scale(1.0, 0.75);
      if (!complexDataSameAsValue(fvindx[2], val_scale, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #6: math::HierarchyFaceDataOpsComplex::scale()\n"
         << "Expected: v2 = " << val_scale << std::endl;
         face_ops->printData(fvindx[2], tbox::plog);
      }

      // Test #7: math::HierarchyFaceDataOpsComplex::add()
      // Expected:  v3 = v0 + v1 = (6.0, 4.5)
      face_ops->add(fvindx[3], fvindx[0], fvindx[1]);
      dcomplex val_add(6.0, 4.5);
      if (!complexDataSameAsValue(fvindx[3], val_add, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #7: math::HierarchyFaceDataOpsComplex::add()\n"
         << "Expected: v3 = " << val_add << std::endl;
         face_ops->printData(fvindx[3], tbox::plog);
      }

      // Reset v0: v0 = (0.0, 4.5)
      face_ops->setToScalar(fvindx[0], dcomplex(0.0, 4.5));

      // Test #8: math::HierarchyFaceDataOpsComplex::subtract()
      // Expected:  v1 = v3 - v0 = (6.0,0.0)
      face_ops->subtract(fvindx[1], fvindx[3], fvindx[0]);
      dcomplex val_sub(6.0, 0.0);
      if (!complexDataSameAsValue(fvindx[1], val_sub, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #8: math::HierarchyFaceDataOpsComplex::subtract()\n"
         << "Expected: v1 = " << val_sub << std::endl;
         face_ops->printData(fvindx[1], tbox::plog);
      }

      // Test #9a: math::HierarchyFaceDataOpsComplex::addScalar()
      // Expected:  v1 = v1 + (0.0,-4.0) = (6.0,-4.0)
      face_ops->addScalar(fvindx[1], fvindx[1], dcomplex(0.0, -4.0));
      dcomplex val_addScalar(6.0, -4.0);
      if (!complexDataSameAsValue(fvindx[1], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9a: math::HierarchyFaceDataOpsComplex::addScalar()\n"
         << "Expected: v1 = " << val_addScalar << std::endl;
         face_ops->printData(fvindx[1], tbox::plog);
      }

      // Test #9b: math::HierarchyFaceDataOpsComplex::addScalar()
      // Expected:  v2 = v2 + (0.0,0.25) = (1.0,1.0)
      face_ops->addScalar(fvindx[2], fvindx[2], dcomplex(0.0, 0.25));
      val_addScalar = dcomplex(1.0, 1.0);
      if (!complexDataSameAsValue(fvindx[2], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9b: math::HierarchyFaceDataOpsComplex::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         face_ops->printData(fvindx[2], tbox::plog);
      }

      // Test #9c: math::HierarchyFaceDataOpsComplex::addScalar()
      // Expected:  v2 = v2 + (3.0,-4.0) = (4.0,-3.0)
      face_ops->addScalar(fvindx[2], fvindx[2], dcomplex(3.0, -4.0));
      val_addScalar = dcomplex(4.0, -3.0);
      if (!complexDataSameAsValue(fvindx[2], val_addScalar, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #9c: math::HierarchyFaceDataOpsComplex::addScalar()\n"
         << "Expected: v2 = " << val_addScalar << std::endl;
         face_ops->printData(fvindx[2], tbox::plog);
      }

      // Reset v3: v3 = (0.5, 0.0)
      face_ops->setToScalar(fvindx[3], dcomplex(0.5, 0.0));

      // Test #10: math::HierarchyFaceDataOpsComplex::multiply()
      // Expected:  v1 = v3 * v1 = (3.0,-2.0)
      face_ops->multiply(fvindx[1], fvindx[3], fvindx[1]);
      dcomplex val_mult(3.0, -2.0);
      if (!complexDataSameAsValue(fvindx[1], val_mult, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #10: math::HierarchyFaceDataOpsComplex::multiply()\n"
         << "Expected: v1 = " << val_mult << std::endl;
         face_ops->printData(fvindx[1], tbox::plog);
      }

      // Test #11: math::HierarchyFaceDataOpsComplex::divide()
      // Expected:  v0 = v2 / v1 = (1.3846153846154,-0.076923076923077)
      face_ops->divide(fvindx[0], fvindx[2], fvindx[1]);
      dcomplex val_div(1.3846153846154, -0.076923076923077);
      if (!complexDataSameAsValue(fvindx[0], val_div, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #11: math::HierarchyFaceDataOpsComplex::divide()\n"
         << "Expected: v0 = " << val_div << std::endl;
         face_ops->printData(fvindx[0], tbox::plog);
      }

      // Test #12: math::HierarchyFaceDataOpsComplex::reciprocal()
      // Expected:  v1 = 1 / v1 = (0.23076923076923, 0.15384615384615)
      face_ops->reciprocal(fvindx[1], fvindx[1]);
      dcomplex val_rec(0.23076923076923, 0.15384615384615);
      if (!complexDataSameAsValue(fvindx[1], val_rec, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #12: math::HierarchyFaceDataOpsComplex::reciprocal()\n"
         << "Expected: v1 = " << val_rec << std::endl;
         face_ops->printData(fvindx[1], tbox::plog);
      }

      // Test #13: Place some bogus values on coarse level and test norms
      boost::shared_ptr<pdat::FaceData<dcomplex> > cdata;

      // set values
      boost::shared_ptr<hier::PatchLevel> level_zero(
         hierarchy->getPatchLevel(0));
      for (hier::PatchLevel::iterator ip(level_zero->begin());
           ip != level_zero->end(); ++ip) {
         patch = level_zero->getPatch(ip());
         cdata = patch->getPatchData(fvindx[2]);
         hier::Index index0(2, 2);
         hier::Index index1(5, 3);
         if (patch->getBox().contains(index0)) {
            (*cdata)(pdat::FaceIndex(index0, pdat::FaceIndex::Y,
                        pdat::FaceIndex::Lower), 0) = dcomplex(100.0, -50.0);
         }
         if (patch->getBox().contains(index1)) {
            (*cdata)(pdat::FaceIndex(index1, pdat::FaceIndex::Y,
                        pdat::FaceIndex::Upper), 0) = dcomplex(-1000.0, 20.0);
         }
      }

      // check values
      bool bogus_value_test_passed = true;
      for (hier::PatchLevel::iterator ipp(level_zero->begin());
           ipp != level_zero->end(); ++ipp) {
         patch = level_zero->getPatch(ipp());
         cdata = patch->getPatchData(fvindx[2]);
         pdat::FaceIndex index0(hier::Index(2,
                                   2), pdat::FaceIndex::Y,
                                pdat::FaceIndex::Lower);
         pdat::FaceIndex index1(hier::Index(5,
                                   3), pdat::FaceIndex::Y,
                                pdat::FaceIndex::Upper);

         // check X axis data
         pdat::FaceIterator cend(cdata->getBox(), pdat::FaceIndex::X, false);
         for (pdat::FaceIterator c(cdata->getBox(), pdat::FaceIndex::X, true);
              c != cend && bogus_value_test_passed;
              ++c) {
            pdat::FaceIndex face_index = *c;

            if (!tbox::MathUtilities<dcomplex>::equalEps((*cdata)(face_index),
                   dcomplex(4.0, -3.0))) {
               bogus_value_test_passed = false;
            }
         }

         // check Y axis data
         pdat::FaceIterator ccend(cdata->getBox(), pdat::FaceIndex::Y, false);
         for (pdat::FaceIterator cc(cdata->getBox(), pdat::FaceIndex::Y, true);
              cc != ccend && bogus_value_test_passed;
              ++cc) {
            pdat::FaceIndex face_index = *cc;

            if (face_index == index0) {
               if (!tbox::MathUtilities<dcomplex>::equalEps((*cdata)(face_index),
                      dcomplex(100.0, -50.0))) {
                  bogus_value_test_passed = false;
               }
            } else {
               if (face_index == (pdat::FaceIndex)index1) {
                  if (!tbox::MathUtilities<dcomplex>::equalEps((*cdata)(
                            face_index),
                         dcomplex(-1000.0, 20.0))) {
                     bogus_value_test_passed = false;
                  }
               } else {
                  if (!tbox::MathUtilities<dcomplex>::equalEps((*cdata)(
                            face_index),
                         dcomplex(4.0, -3.0))) {
                     bogus_value_test_passed = false;
                  }
               }
            }
         }
      }
      if (!bogus_value_test_passed) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #13:  Place some bogus values on coarse level"
         << std::endl;
         face_ops->printData(fvindx[2], tbox::plog);
      }

      // Test norms on patch data with fvindex[2] on hierarchy with bogus values

      // Test #14: math::HierarchyFaceDataOpsComplex::L1Norm() - w/o control weight
      // Expected:  bogus_l1_norm = 2217.003379
      double bogus_l1_norm = face_ops->L1Norm(fvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_l1_norm, 2217.003379)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #14: math::HierarchyFaceDataOpsComplex::L1Norm()"
         << " - w/o control weight\n"
         << "Expected value = 2217.003379, Computed value = "
         << std::setprecision(12) << bogus_l1_norm << std::endl;
      }

      // Test #15: math::HierarchyFaceDataOpsComplex::L1Norm() - w/ control weight
      // Expected:  correct_l1_norm = 5.0
      double correct_l1_norm = face_ops->L1Norm(fvindx[2], fwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(correct_l1_norm, 5.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #15: math::HierarchyFaceDataOpsComplex::L1Norm()"
         << " - w/ control weight\n"
         << "Expected value = 5.0, Computed value = "
         << std::setprecision(12) << correct_l1_norm << std::endl;
      }

      // Test #16: math::HierarchyFaceDataOpsComplex::L2Norm()
      // Expected:  l2_norm = 5.0
      double l2_norm = face_ops->L2Norm(fvindx[2], fwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(l2_norm, 5.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #16: math::HierarchyFaceDataOpsComplex::L2Norm()\n"
         << "Expected value = 5.0, Computed value = "
         << std::setprecision(12) << l2_norm << std::endl;
      }

      // Test #17: math::HierarchyFaceDataOpsComplex::maxNorm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.19998
      double bogus_max_norm = face_ops->maxNorm(fvindx[2]);
      if (!tbox::MathUtilities<double>::equalEps(bogus_max_norm, 1000.19998)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #17: math::HierarchyFaceDataOpsComplex::maxNorm()"
         << " - w/o control weight\n"
         << "Expected value = 1000.19998, Computed value = "
         << std::setprecision(12) << bogus_max_norm << std::endl;
      }

      // Test #18: math::HierarchyFaceDataOpsComplex::maxNorm() - w/o control weight
      // Expected:  bogus_max_norm = 1000.19998
      double max_norm = face_ops->maxNorm(fvindx[2], fwgt_id);
      if (!tbox::MathUtilities<double>::equalEps(max_norm, 5.0)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #18: math::HierarchyFaceDataOpsComplex::maxNorm()"
         << " - w/ control weight\n"
         << "Expected value = 5.0, Computed value = "
         << std::setprecision(12) << max_norm << std::endl;
      }

      // Reset data and test sums, axpy's
      face_ops->setToScalar(fvindx[0], dcomplex(1.0, -3.0));
      face_ops->setToScalar(fvindx[1], dcomplex(2.5, 3.0));
      face_ops->setToScalar(fvindx[2], dcomplex(7.0, 0.0));

      // Test #19: math::HierarchyFaceDataOpsComplex::linearSum()
      // Expected:  v3 = (2.0,5.0)
      face_ops->linearSum(fvindx[3],
         dcomplex(2.0, 0.0), fvindx[1], dcomplex(0.0, -1.0), fvindx[0]);
      dcomplex val_linearSum(2.0, 5.0);
      if (!complexDataSameAsValue(fvindx[3], val_linearSum, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #19: math::HierarchyFaceDataOpsComplex::linearSum()\n"
         << "Expected: v3 = " << val_linearSum << std::endl;
         face_ops->printData(fvindx[3], tbox::plog);
      }

      // Test #20: math::HierarchyFaceDataOpsComplex::axmy()
      // Expected:  v3 = (6.5,12.0)
      face_ops->axmy(fvindx[3], 3.0, fvindx[1], fvindx[0]);
      dcomplex val_axmy(6.5, 12.0);
      if (!complexDataSameAsValue(fvindx[3], val_axmy, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #20: math::HierarchyFaceDataOpsComplex::axmy()\n"
         << "Expected: v3 = " << val_axmy << std::endl;
         face_ops->printData(fvindx[3], tbox::plog);
      }

      // Test #21a: math::HierarchyFaceDataOpsComplex::dot()
      // Expected:  cdot = (17.5,-21.0)
      dcomplex cdot = face_ops->dot(fvindx[2], fvindx[1], fwgt_id);
      dcomplex ans_2_dot_1(17.5, -21.0);
      if (!tbox::MathUtilities<dcomplex>::equalEps(cdot, ans_2_dot_1)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #21a: math::HierarchyFaceDataOpsComplex::dot()\n"
         << "Expected value = (17.5,-21.0), Computed value = "
         << cdot << std::endl;
      }

      // Test #21b: math::HierarchyFaceDataOpsComplex::dot()
      // Expected:  cdot = (17.5,21.0)
      dcomplex cdot2 = face_ops->dot(fvindx[1], fvindx[2], fwgt_id);
      dcomplex ans_1_dot_2(17.5, 21.0);
      if (!tbox::MathUtilities<dcomplex>::equalEps(cdot2, ans_1_dot_2)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #21b: math::HierarchyFaceDataOpsComplex::dot()\n"
         << "Expected value = (17.5,21.0), Computed value = "
         << cdot2 << std::endl;
      }

      // Test #22: math::HierarchyFaceDataOpsComplex::abs()
      // Expected: abs(v0) = 5.0
      face_ops->setToScalar(fvindx[0], dcomplex(4.0, -3.0));
      face_ops->abs(fwgt_id, fvindx[0]);
      if (!doubleDataSameAsValue(fwgt_id, 5.0, hierarchy)) {
         num_failures++;
         tbox::perr
         << "FAILED: - Test #22: math::HierarchyFaceDataOpsComplex::abs()\n"
         << "Expected: abs(v0) = 5.0" << std::endl;
         fwgt_ops->printData(fwgt_id, tbox::plog);
      }

      // deallocate data on hierarchy
      for (ln = 0; ln < 2; ln++) {
         hierarchy->getPatchLevel(ln)->deallocatePatchData(fwgt_id);
         for (iv = 0; iv < NVARS; iv++) {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(fvindx[iv]);
         }
      }

      for (iv = 0; iv < NVARS; iv++) {
         fvar[iv].reset();
      }
      fwgt.reset();

      geometry.reset();
      hierarchy.reset();
      face_ops.reset();
      fwgt_ops.reset();

      if (num_failures == 0) {
         tbox::pout << "\nPASSED:  face cplxtest" << std::endl;
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
complexDataSameAsValue(
   int desc_id,
   dcomplex value,
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
         boost::shared_ptr<pdat::FaceData<dcomplex> > fvdata(
            patch->getPatchData(desc_id));
         pdat::FaceIterator fend(fvdata->getBox(), 1, false);
         for (pdat::FaceIterator f(fvdata->getBox(), 1, true);
              f != fend && test_passed; ++f) {
            pdat::FaceIndex face_index = *f;
            if (!tbox::MathUtilities<dcomplex>::equalEps((*fvdata)(face_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
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
         patch = level->getPatch(ip());
         boost::shared_ptr<pdat::FaceData<double> > fvdata(
            patch->getPatchData(desc_id));

         pdat::FaceIterator fend(fvdata->getBox(), 1, false);
         for (pdat::FaceIterator f(fvdata->getBox(), 1, true);
              f != fend && test_passed; ++f) {
            pdat::FaceIndex face_index = *f;
            if (!tbox::MathUtilities<double>::equalEps((*fvdata)(face_index),
                   value)) {
               test_passed = false;
            }
         }
      }
   }

   return test_passed;
}
