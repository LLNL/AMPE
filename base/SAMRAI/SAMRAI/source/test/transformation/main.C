/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Test program to demonstrate/test the Transformation class
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

// Headers for basic SAMRAI objects used in this code.
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/SideGeometry.h"

using namespace std;

using namespace SAMRAI;

int
testGeometryTransformations(
   const hier::Transformation& transformation,
   const hier::Box& box);

int main(
   int argc,
   char* argv[])
{
   int fail_count = 0;

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      tbox::Dimension dim2(2);
      tbox::Dimension dim3(3);

      hier::BlockId block_zero(0);
      hier::BlockId block_one(1);

      /*
       * For each test, create an arbitrary Box, do a transform and an
       * inverse transform on the Box, and check to see that it has
       * returned to its original state.
       */

      /*
       * 2D trivial transformation.
       */
      {
         hier::Transformation zero_trans(hier::IntVector::getZero(dim2));

         hier::Index index_lo(3, 5);
         hier::Index index_hi(21, 17);
         hier::Box ref_box(index_lo, index_hi, block_zero);
         hier::Box trans_box(ref_box);

         zero_trans.transform(trans_box);
         zero_trans.inverseTransform(trans_box);

         if (!trans_box.isSpatiallyEqual(ref_box)) {
            fail_count++;
            tbox::perr << "FAILED: - 2D trivial tranformation" << endl;
         }

         fail_count += testGeometryTransformations(zero_trans, ref_box);
      }

      /*
       * 3D trivial transformation.
       */
      {
         hier::Transformation zero_trans(hier::IntVector::getZero(dim3));

         hier::Index index_lo(4, 6, 1);
         hier::Index index_hi(15, 38, 12);
         hier::Box ref_box(index_lo, index_hi, block_zero);
         hier::Box trans_box(ref_box);

         zero_trans.transform(trans_box);
         zero_trans.inverseTransform(trans_box);

         if (!trans_box.isSpatiallyEqual(ref_box)) {
            fail_count++;
            tbox::perr << "FAILED: - 3D trivial tranformation" << endl;
         }

         fail_count += testGeometryTransformations(zero_trans, ref_box);
      }

      /*
       * 2D shift transformation with no rotation.
       */
      {
         hier::IntVector offset(dim2);
         offset[0] = 9;
         offset[1] = -7;
         hier::Transformation shift_trans(hier::Transformation::NO_ROTATE,
                                          offset,
                                          block_zero,
                                          block_one);

         hier::Index index_lo(12, 28);
         hier::Index index_hi(33, 44);
         hier::Box ref_box(index_lo, index_hi, block_zero);
         hier::Box trans_box(ref_box);

         shift_trans.transform(trans_box);
         shift_trans.inverseTransform(trans_box);

         if (!trans_box.isSpatiallyEqual(ref_box)) {
            fail_count++;
            tbox::perr << "FAILED: - 2D shift tranformation" << endl;
         }

         fail_count += testGeometryTransformations(shift_trans, ref_box);
      }

      /*
       * 3D shift transformation with no rotation.
       */
      {
         hier::IntVector offset(dim3);
         offset[0] = 13;
         offset[1] = 28;
         offset[2] = 11;
         hier::Transformation shift_trans(hier::Transformation::NO_ROTATE,
                                          offset,
                                          block_zero,
                                          block_one);

         hier::Index index_lo(25, 22, 9);
         hier::Index index_hi(35, 37, 34);
         hier::Box ref_box(index_lo, index_hi, block_zero);
         hier::Box trans_box(ref_box);

         shift_trans.transform(trans_box);
         shift_trans.inverseTransform(trans_box);

         if (!trans_box.isSpatiallyEqual(ref_box)) {
            fail_count++;
            tbox::perr << "FAILED: - 3D shift tranformation" << endl;
         }

         fail_count += testGeometryTransformations(shift_trans, ref_box);
      }

      /*
       * 2D rotation and shift transformations.
       */
      {
         // loop over all 2D rotations
         for (int i = 0; i < 4; i++) {
            hier::Transformation::RotationIdentifier rotate =
               static_cast<hier::Transformation::RotationIdentifier>(i);

            hier::IntVector offset(dim2);
            offset[0] = 14 + i;
            offset[1] = 42 - i;

            hier::Transformation trans(rotate, offset, block_zero, block_one);

            int i_sqr = i * i;
            hier::Index index_lo(11 - i_sqr, 27 - i_sqr);
            hier::Index index_hi(21 + i_sqr, 38 + i_sqr);
            hier::Box ref_box(index_lo, index_hi, block_zero);
            hier::Box trans_box(ref_box);

            trans.transform(trans_box);
            trans.inverseTransform(trans_box);

            if (!trans_box.isSpatiallyEqual(ref_box)) {
               fail_count++;
               tbox::perr << "FAILED: - 2D rotate/shift transformation" << endl;
            }

            fail_count += testGeometryTransformations(trans, ref_box);
         }
      }

      /*
       * 3D rotation and shift transformations.
       */
      {
         // loop over all 3D rotations
         for (int i = 0; i < 24; i++) {
            hier::Transformation::RotationIdentifier rotate =
               static_cast<hier::Transformation::RotationIdentifier>(i);

            hier::IntVector offset(dim3);
            offset[0] = 22 - i;
            offset[1] = 23 + i;
            offset[2] = i - 5;

            hier::Transformation trans(rotate, offset, block_zero, block_one);

            hier::Index index_lo(24 - 2 * i, -8 - 2 * i, 20 - 2 * i);
            hier::Index index_hi(33 + i, 13 + i, 36 + i);
            hier::Box ref_box(index_lo, index_hi, block_zero);
            hier::Box trans_box(ref_box);

            trans.transform(trans_box);
            trans.inverseTransform(trans_box);

            if (!trans_box.isSpatiallyEqual(ref_box)) {
               fail_count++;
               tbox::perr << "FAILED: - 3D rotate/shift transformation" << endl;
            }

            fail_count += testGeometryTransformations(trans, ref_box);
         }
      }

   }

   if (fail_count == 0) {
      tbox::pout << "\nPASSED:  transformation" << std::endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
   return fail_count;
}

int testGeometryTransformations(const hier::Transformation& transformation,
                                const hier::Box& box)
{
   const tbox::Dimension& dim = box.getDim();
   int fail_count = 0;
   hier::BlockId block_zero(0);
   hier::BlockId block_one(1);

   hier::Transformation::RotationIdentifier reverse_rotation =
      hier::Transformation::getReverseRotationIdentifier(
         transformation.getRotation(),
         dim);

   hier::IntVector reverse_shift(dim);
   hier::Transformation::calculateReverseShift(reverse_shift,
      transformation.getOffset(),
      transformation.getRotation());

   hier::Transformation reverse_trans(reverse_rotation, reverse_shift,
                                      block_one, block_zero);

   /*
    * Cell test
    */
   pdat::CellIndex ref_cindex(box.upper());
   pdat::CellIndex trans_cindex(ref_cindex);
   pdat::CellGeometry::transform(trans_cindex, transformation);
   pdat::CellGeometry::transform(trans_cindex, reverse_trans);

   if (trans_cindex != ref_cindex) {
      fail_count++;
      tbox::perr << "FAILED: - CellIndex transformation" << endl;
   }

   /*
    * Node test
    */
   hier::Box ref_node_box(pdat::NodeGeometry::toNodeBox(box));
   hier::Box trans_node_box(ref_node_box);
   pdat::NodeGeometry::transform(trans_node_box, transformation);
   pdat::NodeGeometry::transform(trans_node_box, reverse_trans);
   if (!trans_node_box.isSpatiallyEqual(ref_node_box)) {
      fail_count++;
      tbox::perr << "FAILED: - Node box transformation" << endl;
   }

   pdat::NodeIndex ref_nindex(box.upper(), hier::IntVector::getOne(dim));
   pdat::NodeIndex trans_nindex(ref_nindex);
   pdat::NodeGeometry::transform(trans_nindex, transformation);
   pdat::NodeGeometry::transform(trans_nindex, reverse_trans);

   if (trans_nindex != ref_nindex) {
      fail_count++;
      tbox::perr << "FAILED: - NodeIndex transformation" << endl;
   }

   /*
    * Side test
    */
   for (int d = 0; d < dim.getValue(); d++) {
      hier::Box ref_side_box(pdat::SideGeometry::toSideBox(box, d));
      hier::Box trans_side_box(ref_side_box);
      int direction = d;
      pdat::SideGeometry::transform(trans_side_box, direction, transformation);
      pdat::SideGeometry::transform(trans_side_box, direction, reverse_trans);
      if (!trans_side_box.isSpatiallyEqual(ref_side_box)) {
         fail_count++;
         tbox::perr << "FAILED: - Side box transformation" << endl;
      }

      pdat::SideIndex ref_sindex(box.upper(), d, 0);
      pdat::SideIndex trans_sindex(ref_sindex);
      pdat::SideGeometry::transform(trans_sindex, transformation);
      pdat::SideGeometry::transform(trans_sindex, reverse_trans);

      if (trans_sindex != ref_sindex) {
         fail_count++;
         tbox::perr << "FAILED: - SideIndex transformation" << endl;
      }
   }

   /*
    * Edge test
    */
   for (int d = 0; d < dim.getValue(); d++) {
      hier::Box ref_edge_box(pdat::EdgeGeometry::toEdgeBox(box, d));
      hier::Box trans_edge_box(ref_edge_box);
      int direction = d;
      pdat::EdgeGeometry::transform(trans_edge_box, direction, transformation);
      pdat::EdgeGeometry::transform(trans_edge_box, direction, reverse_trans);
      if (!trans_edge_box.isSpatiallyEqual(ref_edge_box)) {
         fail_count++;
         tbox::perr << "FAILED: - Edge box transformation" << endl;
      }

      pdat::EdgeIndex ref_eindex(box.upper(), d, 0);
      pdat::EdgeIndex trans_eindex(ref_eindex);
      pdat::EdgeGeometry::transform(trans_eindex, transformation);
      pdat::EdgeGeometry::transform(trans_eindex, reverse_trans);

      if (trans_eindex != ref_eindex) {
         fail_count++;
         tbox::perr << "FAILED: - EdgeIndex transformation" << endl;
      }
   }

   /*
    * Face test
    */
   for (int d = 0; d < dim.getValue(); d++) {
      hier::Box ref_face_box(pdat::FaceGeometry::toFaceBox(box, d));
      hier::Box trans_face_box(ref_face_box);
      int direction = d;
      pdat::FaceGeometry::transform(trans_face_box, direction, transformation);
      pdat::FaceGeometry::transform(trans_face_box, direction, reverse_trans);
      if (!trans_face_box.isSpatiallyEqual(ref_face_box)) {
         fail_count++;
         tbox::perr << "FAILED: - Face box transformation" << endl;
      }

      pdat::FaceIndex ref_findex(box.upper(), d, 0);
      pdat::FaceIndex trans_findex(ref_findex);
      pdat::FaceGeometry::transform(trans_findex, transformation);
      pdat::FaceGeometry::transform(trans_findex, reverse_trans);

      if (trans_findex != ref_findex) {
         fail_count++;
         tbox::perr << "FAILED: - FaceIndex transformation" << endl;
      }
   }

   return fail_count;
}
