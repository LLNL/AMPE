/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated outerface centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceData_C
#define included_pdat_OuterfaceData_C

#include "SAMRAI/pdat/OuterfaceData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int OuterfaceData<TYPE>::PDAT_OUTERFACEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for outerface data objects.  The
 * constructor simply initializes data variables and sets up the
 * array data.
 *
 *************************************************************************
 */

template<class TYPE>
OuterfaceData<TYPE>::OuterfaceData(
   const hier::Box& box,
   int depth):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{

   TBOX_ASSERT(depth > 0);

   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::Box& ghosts = getGhostBox();
      const hier::Box facebox = FaceGeometry::toFaceBox(ghosts, d);
      hier::Box outerfacebox = facebox;
      outerfacebox.upper(0) = facebox.lower(0);
      d_data[d][0].initializeArray(outerfacebox, depth);
      outerfacebox.lower(0) = facebox.upper(0);
      outerfacebox.upper(0) = facebox.upper(0);
      d_data[d][1].initializeArray(outerfacebox, depth);
   }
}

template<class TYPE>
OuterfaceData<TYPE>::~OuterfaceData()
{
}

/*
 *************************************************************************
 *
 * The following are private and cannot be used, but they are defined
 * here for compilers that require that every template declaration have
 * a definition (a stupid requirement, if you ask me).
 *
 *************************************************************************
 */

template<class TYPE>
OuterfaceData<TYPE>::OuterfaceData(
   const OuterfaceData<TYPE>& foo):
   hier::PatchData(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<class TYPE>
void
OuterfaceData<TYPE>::operator = (
   const OuterfaceData<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
int
OuterfaceData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE*
OuterfaceData<TYPE>::getPointer(
   int face_normal,
   int side,
   int d)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_data[face_normal][side].getPointer(d);
}

template<class TYPE>
const TYPE*
OuterfaceData<TYPE>::getPointer(
   int face_normal,
   int side,
   int d) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_data[face_normal][side].getPointer(d);
}

template<class TYPE>
ArrayData<TYPE>&
OuterfaceData<TYPE>::getArrayData(
   int face_normal,
   int side)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[face_normal][side];
}

template<class TYPE>
const ArrayData<TYPE>&
OuterfaceData<TYPE>::getArrayData(
   int face_normal,
   int side) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[face_normal][side];
}

template<class TYPE>
TYPE&
OuterfaceData<TYPE>::operator () (
   const FaceIndex& i,
   int side,
   int depth)
{
   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][side](i, depth);
}

template<class TYPE>
const TYPE&
OuterfaceData<TYPE>::operator () (
   const FaceIndex& i,
   int side,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis][side](i, depth);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between an outerface patch data type (source) and
 * a face patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const FaceData<TYPE> * const t_src =
      dynamic_cast<const FaceData<TYPE> *>(&src);

   TBOX_ASSERT(t_src != NULL);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const ArrayData<TYPE>& face_array = t_src->getArrayData(axis);
      for (int loc = 0; loc < 2; loc++) {
         ArrayData<TYPE>& oface_array = d_data[axis][loc];
         oface_array.copy(face_array, oface_array.getBox());
      }
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   FaceData<TYPE>* t_dst =
      dynamic_cast<FaceData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != NULL);

   for (int d = 0; d < getDim().getValue(); d++) {
      t_dst->getArrayData(d).copy(d_data[d][0], d_data[d][0].getBox());
      t_dst->getArrayData(d).copy(d_data[d][1], d_data[d][1].getBox());
   }
}

/*
 *************************************************************************
 *
 * Copy data from the source into the destination according to the
 * overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{

   NULL_USE(src);
   NULL_USE(overlap);

   TBOX_ERROR("Copy with outerface as destination is not defined yet...");
}

template<class TYPE>
void
OuterfaceData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   FaceData<TYPE>* t_dst =
      dynamic_cast<FaceData<TYPE> *>(&dst);
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);

   const hier::Transformation& transformation = t_overlap->getTransformation();
   TBOX_ASSERT(transformation.getRotation() == hier::Transformation::NO_ROTATE);

   const hier::IntVector& src_offset = transformation.getOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      hier::IntVector face_offset(src_offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); i++) {
            face_offset(i) = src_offset((d + i) % getDim().getValue());
         }
      }
      hier::Transformation face_transform(hier::Transformation::NO_ROTATE,
                                          face_offset,
                                          getBox().getBlockId(),
                                          t_dst->getBox().getBlockId());
 
      const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
      t_dst->getArrayData(d).copy(d_data[d][0], box_list, face_transform);
      t_dst->getArrayData(d).copy(d_data[d][1], box_list, face_transform);
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy from a face data object to this outerface data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::copyDepth(
   int dst_depth,
   const FaceData<TYPE>& src,
   int src_depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const ArrayData<TYPE>& src_face_array = src.getArrayData(axis);
      for (int loc = 0; loc < 2; loc++) {
         ArrayData<TYPE>& dst_oface_array = d_data[axis][loc];
         dst_oface_array.copyDepth(dst_depth,
            src_face_array,
            src_depth,
            dst_oface_array.getBox());
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy to a face data object from this outerface data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::copyDepth2(
   int dst_depth,
   FaceData<TYPE>& dst,
   int src_depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      ArrayData<TYPE>& dst_face_array = dst.getArrayData(axis);
      for (int loc = 0; loc < 2; loc++) {
         const ArrayData<TYPE>& src_oface_array = d_data[axis][loc];
         dst_face_array.copyDepth(dst_depth,
            src_oface_array,
            src_depth,
            src_oface_array.getBox());
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the buffer space needed to pack/unpack messages on the box
 * region using the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
bool
OuterfaceData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
int
OuterfaceData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();

   int size = 0;
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& boxlist = t_overlap->getDestinationBoxContainer(d);
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); i++) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }
      size += d_data[d][0].getDataStreamSize(boxlist, face_offset);
      size += d_data[d][1].getDataStreamSize(boxlist, face_offset);
   }
   return size;
}

/*
 *************************************************************************
 *
 * Pack/unpack data into/out of the message streams using the index
 * space in the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);

      if (boxes.size() > 0) {
         hier::IntVector face_offset(offset);
         if (d > 0) {
            for (int i = 0; i < getDim().getValue(); i++) {
               face_offset(i) = offset((d + i) % getDim().getValue());
            }
         }

         hier::Transformation face_transform(hier::Transformation::NO_ROTATE,
                                             face_offset,
                                             getBox().getBlockId(),
                                             boxes.begin()->getBlockId());

         for (hier::BoxContainer::const_iterator b(boxes);
              b != boxes.end(); ++b) {
            hier::Box src_box(*b);
            face_transform.inverseTransform(src_box);
            for (int f = 0; f < 2; f++) {
               hier::Box intersect(src_box * d_data[d][f].getBox());
               if (!intersect.empty()) {
                  face_transform.transform(intersect); 
                  d_data[d][f].packStream(stream,
                     intersect,
                     face_transform);
               }
            }
         }
      }

   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      hier::IntVector face_offset(offset);
      if (d > 0) {
         for (int i = 0; i < getDim().getValue(); i++) {
            face_offset(i) = offset((d + i) % getDim().getValue());
         }
      }

      for (hier::BoxContainer::const_iterator b(boxes);
           b != boxes.end(); ++b) {
         for (int f = 0; f < 2; f++) {
            const hier::Box intersect = (*b) * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].unpackStream(stream, intersect, face_offset);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  outerface centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
OuterfaceData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth)
{
   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   for (int d = 0; d < box.getDim().getValue(); d++) {
      hier::Box lower = FaceGeometry::toFaceBox(box, d);
      hier::Box upper = FaceGeometry::toFaceBox(box, d);
      lower.upper(d) = box.lower(d);
      upper.lower(d) = box.upper(d);
      size += ArrayData<TYPE>::getSizeOfData(lower, depth);
      size += ArrayData<TYPE>::getSizeOfData(upper, depth);
   }
   return size;
}

/*
 *************************************************************************
 *
 * Fill the outerface centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i][0].fill(t, d);
      d_data[i][1].fill(t, d);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i][0].fill(t, FaceGeometry::toFaceBox(box, i), d);
      d_data[i][1].fill(t, FaceGeometry::toFaceBox(box, i), d);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i][0].fillAll(t);
      d_data[i][1].fillAll(t);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i][0].fillAll(t, FaceGeometry::toFaceBox(box, i));
      d_data[i][1].fillAll(t, FaceGeometry::toFaceBox(box, i));
   }
}

/*
 *************************************************************************
 *
 * Print routines for outerface centered arrays.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int d = 0; d < d_depth; d++) {
      print(box, d, os, prec);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int face_normal = 0; face_normal < getDim().getValue(); face_normal++) {
      os << "Array face normal = " << face_normal << std::endl;
      for (int side = 0; side < 2; side++) {
         os << "side = " << ((side == 0) ? "lower" : "upper") << std::endl;
         printAxisFace(face_normal, side, box, depth, os, prec);
      }
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::printAxisFace(
   int face_normal,
   int side,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxisFace(face_normal, side, box, d, os, prec);
   }
}

template<class TYPE>
void
OuterfaceData<TYPE>::printAxisFace(
   int face_normal,
   int side,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   const hier::Box facebox =
      FaceGeometry::toFaceBox(box, face_normal);
   const hier::Box region =
      facebox * d_data[face_normal][side].getBox();
   os.precision(prec);
   hier::Box::iterator iend(region, false);
   for (hier::Box::iterator i(region, true); i != iend; ++i) {
      os << "array" << *i << " = "
         << d_data[face_normal][side](*i, depth) << std::endl;
      os << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.  If so,
 * reads in d_depth from the database.  Then has each item in d_data
 * read in its data from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::getSpecializedFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("PDAT_OUTERFACEDATA_VERSION");
   if (ver != PDAT_OUTERFACEDATA_VERSION) {
      TBOX_ERROR(
         "OuterfaceData<getDim()>::getSpecializedFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i)
         + "_1";
      array_database = database->getDatabase(array_name);
      (d_data[i][0]).getFromDatabase(array_database);

      array_name = "d_data" + tbox::Utilities::intToString(i) + "_2";
      array_database = database->getDatabase(array_name);
      (d_data[i][1]).getFromDatabase(array_database);
   }
}

/*
 *************************************************************************
 *
 * Writes out class version number, d_depth to the database.
 * Then has each item in d_data write out its data to the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuterfaceData<TYPE>::putSpecializedToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{

   TBOX_ASSERT(database);

   database->putInteger("PDAT_OUTERFACEDATA_VERSION",
      PDAT_OUTERFACEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i)
         + "_1";
      array_database = database->putDatabase(array_name);
      (d_data[i][0]).putUnregisteredToDatabase(array_database);

      array_name = "d_data" + tbox::Utilities::intToString(i) + "_2";
      array_database = database->putDatabase(array_name);
      (d_data[i][1]).putUnregisteredToDatabase(array_database);
   }
}

}
}

#endif
