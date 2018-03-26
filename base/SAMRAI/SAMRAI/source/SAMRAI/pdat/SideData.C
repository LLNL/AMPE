/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated side centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_SideData_C
#define included_pdat_SideData_C

#include "SAMRAI/pdat/SideData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int SideData<TYPE>::PDAT_SIDEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for side data objects.  The constructor
 * simply initializes data variables and sets up the array data.
 *
 *************************************************************************
 */

template<class TYPE>
SideData<TYPE>::SideData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts,
   const hier::IntVector& directions):
   hier::PatchData(box, ghosts),
   d_depth(depth),
   d_directions(directions)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(box, ghosts, directions);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
   TBOX_ASSERT(directions.min() >= 0);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < getDim().getValue(); d++) {
      if (d_directions(d)) {
         const hier::Box side =
            SideGeometry::toSideBox(getGhostBox(), d);
         d_data[d].initializeArray(side, depth);
      } else {
         d_data[d].invalidateArray(dim);
      }
   }
}

template<class TYPE>
SideData<TYPE>::SideData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):

   hier::PatchData(box, ghosts),
   d_depth(depth),
   d_directions(hier::IntVector::getOne(box.getDim()))
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);
   TBOX_ASSERT(d_directions.min() >= 0);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < getDim().getValue(); d++) {
      if (d_directions(d)) {
         const hier::Box side =
            SideGeometry::toSideBox(getGhostBox(), d);
         d_data[d].initializeArray(side, depth);
      } else {
         d_data[d].invalidateArray(dim);
      }
   }
}

template<class TYPE>
SideData<TYPE>::~SideData()
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
SideData<TYPE>::SideData(
   const SideData<TYPE>& foo):
   hier::PatchData(foo.getBox(),
                   foo.getGhostCellWidth()),
   d_directions(tbox::Dimension(getDim()))
{
   NULL_USE(foo);
}

template<class TYPE>
void
SideData<TYPE>::operator = (
   const SideData<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
const hier::IntVector&
SideData<TYPE>::getDirectionVector() const
{
   return d_directions;
}

template<class TYPE>
int
SideData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE*
SideData<TYPE>::getPointer(
   int side_normal,
   int depth)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal].getPointer(depth);
}

template<class TYPE>
const TYPE*
SideData<TYPE>::getPointer(
   int side_normal,
   int depth) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal].getPointer(depth);
}

template<class TYPE>
TYPE&
SideData<TYPE>::operator () (
   const SideIndex& i,
   int depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT(d_directions(axis));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis](i, depth);
}

template<class TYPE>
const TYPE&
SideData<TYPE>::operator () (
   const SideIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   const int axis = i.getAxis();

   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));
   TBOX_ASSERT(d_directions(axis));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[axis](i, depth);
}

template<class TYPE>
ArrayData<TYPE>&
SideData<TYPE>::getArrayData(
   int side_normal)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));

   return d_data[side_normal];
}

template<class TYPE>
const ArrayData<TYPE>&
SideData<TYPE>::getArrayData(
   int side_normal) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT(d_directions(side_normal));

   return d_data[side_normal];
}

/*
 *************************************************************************
 *
 * Perform a fast copy between two side centered arrays where their
 * index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
SideData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, src);

   const SideData<TYPE>* t_src =
      dynamic_cast<const SideData<TYPE> *>(&src);

   if (t_src == NULL) {
      src.copy2(*this);
   } else {

      TBOX_ASSERT(t_src->getDirectionVector() == d_directions);

      for (int d = 0; d < getDim().getValue(); d++) {
         if (d_directions(d)) {
            const hier::Box box =
               d_data[d].getBox() * t_src->d_data[d].getBox();
            if (!box.empty()) {
               d_data[d].copy(t_src->d_data[d], box);
            }
         }
      }
   }
}

template<class TYPE>
void
SideData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, dst);

   SideData<TYPE>* t_dst =
      dynamic_cast<SideData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_dst->getDirectionVector() == d_directions);

   for (int d = 0; d < getDim().getValue(); d++) {
      if (d_directions(d)) {
         const hier::Box box = d_data[d].getBox() * t_dst->d_data[d].getBox();
         if (!box.empty()) {
            t_dst->d_data[d].copy(d_data[d], box);
         }
      }
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
SideData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, src);

   const SideData<TYPE>* t_src =
      dynamic_cast<const SideData<TYPE> *>(&src);
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   if ((t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {

      TBOX_ASSERT(t_src->getDirectionVector() == d_directions);

      const hier::Transformation& transformation =
         t_overlap->getTransformation();         
      if (transformation.getRotation() == hier::Transformation::NO_ROTATE) {

         for (int d = 0; d < getDim().getValue(); d++) {
            if (d_directions(d)) {
               const hier::BoxContainer& box_list =
                  t_overlap->getDestinationBoxContainer(d);
               d_data[d].copy(t_src->d_data[d], box_list, transformation);
            }
         }
      } else {
         copyWithRotation(*t_src, *t_overlap);
      }
   }
}

template<class TYPE>
void
SideData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, dst);

   SideData<TYPE>* t_dst =
      dynamic_cast<SideData<TYPE> *>(&dst);
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);
   TBOX_ASSERT(t_dst->getDirectionVector() == d_directions);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {

      const hier::IntVector& src_offset = t_overlap->getSourceOffset();
      for (int d = 0; d < getDim().getValue(); d++) {
         if (d_directions(d)) {
            const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
            t_dst->d_data[d].copy(d_data[d], box_list, src_offset);
         }
      }
   } else {
      t_dst->copyWithRotation(*this, *t_overlap);
   }
}

template<class TYPE>
void
SideData<TYPE>::copyOnBox(
   const SideData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const hier::Box side_box = SideGeometry::toSideBox(box, axis);
      d_data[axis].copy(src.getArrayData(axis), side_box);
   }

}

template<class TYPE>
void
SideData<TYPE>::copyWithRotation(
   const SideData<TYPE>& src,
   const SideOverlap& overlap)
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(src.getDim());
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Box rotatebox(src.getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   rotatebox.getBlockId(),
                                   getBox().getBlockId());

   for (int i = 0; i < dim.getValue(); i++) {
      if (d_directions(i)) {
         const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer(i);

         hier::Box side_rotatebox(SideGeometry::toSideBox(rotatebox, i));

         for (hier::BoxContainer::const_iterator bi(overlap_boxes);
              bi != overlap_boxes.end(); ++bi) {
            const hier::Box& overlap_box = *bi;

            const hier::Box copybox(side_rotatebox * overlap_box);

            if (!copybox.empty()) {
               const int depth = ((getDepth() < src.getDepth()) ?
                                  getDepth() : src.getDepth());

               hier::Box::iterator ciend(copybox, false);
               for (hier::Box::iterator ci(copybox, true); ci != ciend; ++ci) {

                  SideIndex dst_index(*ci, 0, 0);
                  dst_index.setAxis(i);
                  SideIndex src_index(dst_index);
                  SideGeometry::transform(src_index, back_trans);

                  for (int d = 0; d < depth; d++) {
                     (*this)(dst_index, d) = src(src_index, d);
                  }
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy between two arrays at the
 * specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
SideData<TYPE>::copyDepth(
   int dst_depth,
   const SideData<TYPE>& src,
   int src_depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, src);
   TBOX_ASSERT(src.d_directions == d_directions);

   for (int d = 0; d < getDim().getValue(); d++) {
      if (d_directions(d)) {
         const hier::Box box = d_data[d].getBox() * src.d_data[d].getBox();
         if (!box.empty()) {
            d_data[d].copyDepth(dst_depth, src.d_data[d], src_depth, box);
         }
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
SideData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
int
SideData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();

   int size = 0;
   for (int d = 0; d < getDim().getValue(); d++) {
      if (d_directions(d)) {
         size +=
            d_data[d].getDataStreamSize(t_overlap->getDestinationBoxContainer(d),
               offset);
      }
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
SideData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {

      const hier::IntVector& offset = t_overlap->getSourceOffset();
      for (int d = 0; d < getDim().getValue(); d++) {
         if (d_directions(d)) {
            const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
            if (boxes.size() > 0) {
               d_data[d].packStream(stream, boxes, offset);
            }
         }
      }
   } else {
      packWithRotation(stream, *t_overlap);
   }
}

template<class TYPE>
void
SideData<TYPE>::packWithRotation(
   tbox::MessageStream& stream,
   const SideOverlap& overlap) const
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(getDim());
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Box rotatebox(getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   rotatebox.getBlockId(),
                                   getBox().getBlockId());

   const int depth = getDepth();

   for (int i = 0; i < dim.getValue(); i++) {
      if (d_directions(i)) {
         const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer(i);

         const int size = depth * overlap_boxes.getTotalSizeOfBoxes();
         tbox::Array<TYPE> buffer(size);

         hier::Box side_rotatebox(SideGeometry::toSideBox(rotatebox, i));

         int buf_count = 0;
         for (hier::BoxContainer::const_iterator bi(overlap_boxes);
              bi != overlap_boxes.end(); ++bi) {
            const hier::Box& overlap_box = *bi;

            const hier::Box copybox(side_rotatebox * overlap_box);

            if (!copybox.empty()) {

               for (int d = 0; d < depth; d++) {

                  hier::Box::iterator ciend(copybox, false);
                  for (hier::Box::iterator ci(copybox, true); ci != ciend; ++ci) {

                     SideIndex src_index(*ci, 0, 0);
                     src_index.setAxis(i);
                     SideGeometry::transform(src_index, back_trans);

                     buffer[buf_count] = (*this)(src_index, d);
                     buf_count++;
                  }
               }
            }
         }
         stream.pack(buffer.getPointer(), size);
      }
   }
}

template<class TYPE>
void
SideData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      if (d_directions(d)) {
         const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
         if (boxes.size() > 0) {
            d_data[d].unpackStream(stream, boxes, offset);
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  side centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
SideData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts,
   const hier::IntVector& directions)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(box, ghosts, directions);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(directions.min() >= 0);

   size_t size = 0;
   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   for (int d = 0; d < box.getDim().getValue(); d++) {
      if (directions(d)) {
         const hier::Box side_box = SideGeometry::toSideBox(ghost_box, d);
         size += ArrayData<TYPE>::getSizeOfData(side_box, depth);
      }
   }
   return size;
}

/*
 *************************************************************************
 *
 * Fill the side centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
SideData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_directions(i)) {
         d_data[i].fill(t, d);
      }
   }
}

template<class TYPE>
void
SideData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_directions(i)) {
         d_data[i].fill(t, SideGeometry::toSideBox(box, i), d);
      }
   }
}

template<class TYPE>
void
SideData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_directions(i)) {
         d_data[i].fillAll(t);
      }
   }
}

template<class TYPE>
void
SideData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, box);

   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_directions(i)) {
         d_data[i].fillAll(t, SideGeometry::toSideBox(box, i));
      }
   }
}

/*
 *************************************************************************
 *
 * Print side centered data.  Note:  makes call to specialized printAxis
 * routine in SideDataSpecialized.C
 *
 *************************************************************************
 */

template<class TYPE>
void
SideData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, box);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      os << "Array side normal = " << axis << std::endl;
      printAxis(axis, box, os, prec);
   }
}

template<class TYPE>
void
SideData<TYPE>::print(
   const hier::Box& box,
   int d,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      os << "Array side normal = " << axis << std::endl;
      printAxis(axis, box, d, os, prec);
   }
}

template<class TYPE>
void
SideData<TYPE>::printAxis(
   int axis,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, box);
   TBOX_ASSERT((axis >= 0) && (axis < getDim().getValue()));

   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxis(axis, box, d, os, prec);
   }
}

template<class TYPE>
void
SideData<TYPE>::printAxis(
   int side_normal,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(d_directions, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));

   os.precision(prec);
   if (d_directions(side_normal)) {
      SideIterator iend(box, side_normal, false);
      for (SideIterator i(box, side_normal, true); i != iend; ++i) {
         os << "array" << *i << " = "
            << d_data[side_normal](*i, depth) << std::endl << std::flush;
      }
   } else {
      os << "No side data in " << side_normal << " side normal direction"
         << std::endl << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.  If so,
 * reads in the d_depth data member to the database.  Then tells
 * d_data to read itself in from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
SideData<TYPE>::getSpecializedFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("PDAT_SIDEDATA_VERSION");
   if (ver != PDAT_SIDEDATA_VERSION) {
      TBOX_ERROR("SideData<DIM>::getSpecializedFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_directions(i)) {
         std::string array_name = "d_data" + tbox::Utilities::intToString(i);
         array_database = database->getDatabase(array_name);
         (d_data[i]).getFromDatabase(array_database);
      }
   }
}

/*
 *************************************************************************
 *
 * Write out the class version number, d_depth data member to the
 * database.  Then tells d_data to write itself to the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
SideData<TYPE>::putSpecializedToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   TBOX_ASSERT(database);

   database->putInteger("PDAT_SIDEDATA_VERSION", PDAT_SIDEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_directions(i)) {
         std::string array_name = "d_data" + tbox::Utilities::intToString(i);
         array_database = database->putDatabase(array_name);
         (d_data[i]).putUnregisteredToDatabase(array_database);
      }
   }
}

}
}

#endif
