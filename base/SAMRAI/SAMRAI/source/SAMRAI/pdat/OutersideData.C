/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated outerside centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OutersideData_C
#define included_pdat_OutersideData_C

#include "SAMRAI/pdat/OutersideData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/tbox/Utilities.h"
#include <stdio.h>

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int OutersideData<TYPE>::PDAT_OUTERSIDEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for outerside data objects.  The
 * constructor simply initializes data variables and sets up the
 * array data.
 *
 *************************************************************************
 */

template<class TYPE>
OutersideData<TYPE>::OutersideData(
   const hier::Box& box,
   int depth):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{
   TBOX_ASSERT(depth > 0);

   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::Box& ghosts = getGhostBox();
      const hier::Box sidebox = SideGeometry::toSideBox(ghosts, d);
      hier::Box outersidebox = sidebox;
      outersidebox.upper(d) = sidebox.lower(d);
      d_data[d][0].initializeArray(outersidebox, depth);
      outersidebox.lower(d) = sidebox.upper(d);
      outersidebox.upper(d) = sidebox.upper(d);
      d_data[d][1].initializeArray(outersidebox, depth);
   }
}

template<class TYPE>
OutersideData<TYPE>::~OutersideData()
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
OutersideData<TYPE>::OutersideData(
   const OutersideData<TYPE>& foo):
   hier::PatchData(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<class TYPE>
void
OutersideData<TYPE>::operator = (
   const OutersideData<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
int
OutersideData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
TYPE*
OutersideData<TYPE>::getPointer(
   int side_normal,
   int side,
   int depth)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal][side].getPointer(depth);
}

template<class TYPE>
const TYPE*
OutersideData<TYPE>::getPointer(
   int side_normal,
   int side,
   int depth) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[side_normal][side].getPointer(depth);
}

template<class TYPE>
ArrayData<TYPE>&
OutersideData<TYPE>::getArrayData(
   int side_normal,
   int side)
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[side_normal][side];
}

template<class TYPE>
const ArrayData<TYPE>&
OutersideData<TYPE>::getArrayData(
   int side_normal,
   int side) const
{
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[side_normal][side];
}

template<class TYPE>
TYPE&
OutersideData<TYPE>::operator () (
   const SideIndex& i,
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
OutersideData<TYPE>::operator () (
   const SideIndex& i,
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
 * Perform a fast copy between an outerside patch data type (source) and
 * a side patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const SideData<TYPE> * const t_src =
      dynamic_cast<const SideData<TYPE> *>(&src);

   TBOX_ASSERT(t_src != NULL);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const ArrayData<TYPE>& side_array = t_src->getArrayData(axis);
      for (int loc = 0; loc < 2; loc++) {
         ArrayData<TYPE>& oside_array = d_data[axis][loc];
         oside_array.copy(side_array, oside_array.getBox());
      }
   }

}

template<class TYPE>
void
OutersideData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   SideData<TYPE>* t_dst =
      dynamic_cast<SideData<TYPE> *>(&dst);

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
OutersideData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   NULL_USE(src);
   NULL_USE(overlap);

   TBOX_ERROR("Copy with outerside as destination is not defined yet...");
}

template<class TYPE>
void
OutersideData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   SideData<TYPE>* t_dst =
      dynamic_cast<SideData<TYPE> *>(&dst);
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer(d);
      t_dst->getArrayData(d).copy(d_data[d][0], box_list, src_offset);
      t_dst->getArrayData(d).copy(d_data[d][1], box_list, src_offset);
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy from a side data object to this outerside data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::copyDepth(
   int dst_depth,
   const SideData<TYPE>& src,
   int src_depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      const ArrayData<TYPE>& src_side_array = src.getArrayData(axis);
      for (int loc = 0; loc < 2; loc++) {
         ArrayData<TYPE>& dst_oside_array = d_data[axis][loc];
         dst_oside_array.copyDepth(dst_depth,
            src_side_array,
            src_depth,
            dst_oside_array.getBox());
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy to a side data object from this outerside data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::copyDepth2(
   int dst_depth,
   SideData<TYPE>& dst,
   int src_depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      ArrayData<TYPE>& dst_side_array = dst.getArrayData(axis);
      for (int loc = 0; loc < 2; loc++) {
         const ArrayData<TYPE>& src_oside_array = d_data[axis][loc];
         dst_side_array.copyDepth(dst_depth,
            src_oside_array,
            src_depth,
            src_oside_array.getBox());
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
OutersideData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
int
OutersideData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();

   int size = 0;
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& boxlist = t_overlap->getDestinationBoxContainer(d);
      size += d_data[d][0].getDataStreamSize(boxlist, src_offset);
      size += d_data[d][1].getDataStreamSize(boxlist, src_offset);
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
OutersideData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      for (hier::BoxContainer::const_iterator b(boxes);
           b != boxes.end(); ++b) {
         const hier::Box src_box = hier::Box::shift(*b, -src_offset);
         for (int f = 0; f < 2; f++) {
            const hier::Box intersect = src_box * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].packStream(stream,
                  hier::Box::shift(intersect, src_offset),
                  src_offset);
            }
         }
      }
   }
}

template<class TYPE>
void
OutersideData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const SideOverlap* t_overlap =
      dynamic_cast<const SideOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(d);
      for (hier::BoxContainer::const_iterator b(boxes);
           b != boxes.end(); ++b) {
         for (int f = 0; f < 2; f++) {
            const hier::Box intersect = (*b) * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].unpackStream(stream, intersect, src_offset);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  outerside centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
OutersideData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth)
{
   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   for (int d = 0; d < box.getDim().getValue(); d++) {
      hier::Box lower = SideGeometry::toSideBox(box, d);
      hier::Box upper = SideGeometry::toSideBox(box, d);
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
 * Fill the outerside centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::fill(
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
OutersideData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i][0].fill(t, SideGeometry::toSideBox(box, i), d);
      d_data[i][1].fill(t, SideGeometry::toSideBox(box, i), d);
   }
}

template<class TYPE>
void
OutersideData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i][0].fillAll(t);
      d_data[i][1].fillAll(t);
   }
}

template<class TYPE>
void
OutersideData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int i = 0; i < getDim().getValue(); i++) {
      d_data[i][0].fillAll(t, SideGeometry::toSideBox(box, i));
      d_data[i][1].fillAll(t, SideGeometry::toSideBox(box, i));
   }
}

/*
 *************************************************************************
 *
 * Print routines for outerside centered arrays.
 *
 *************************************************************************
 */

template<class TYPE>
void
OutersideData<TYPE>::print(
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
OutersideData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int side_normal = 0; side_normal < getDim().getValue(); side_normal++) {
      os << "Array side normal  = " << side_normal << std::endl;
      for (int side = 0; side < 2; side++) {
         os << "side = " << ((side == 0) ? "lower" : "upper") << std::endl;
         printAxisSide(side_normal, side, box, depth, os, prec);
      }
   }
}

template<class TYPE>
void
OutersideData<TYPE>::printAxisSide(
   int side_normal,
   int side,
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      printAxisSide(side_normal, side, box, d, os, prec);
   }
}

template<class TYPE>
void
OutersideData<TYPE>::printAxisSide(
   int side_normal,
   int side,
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));
   TBOX_ASSERT((side_normal >= 0) && (side_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   const hier::Box sidebox =
      SideGeometry::toSideBox(box, side_normal);
   const hier::Box region =
      sidebox * d_data[side_normal][side].getBox();
   os.precision(prec);
   hier::Box::iterator iend(region, false);
   for (hier::Box::iterator i(region, true); i != iend; ++i) {
      os << "array" << *i << " = "
         << d_data[side_normal][side](*i, depth) << std::endl;
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
OutersideData<TYPE>::getSpecializedFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("PDAT_OUTERSIDEDATA_VERSION");
   if (ver != PDAT_OUTERSIDEDATA_VERSION) {
      TBOX_ERROR("OutersideData<DIM>::getSpecializedFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i)
         + "_1";
      array_database = database->getDatabase(array_name);
      (d_data[i][0]).getFromDatabase(array_database);

      array_name = "d_data%d_" + tbox::Utilities::intToString(i) + "_2";
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
OutersideData<TYPE>::putSpecializedToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   TBOX_ASSERT(database);

   database->putInteger("PDAT_OUTERSIDEDATA_VERSION",
      PDAT_OUTERSIDEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      std::string array_name = "d_data%d_" + tbox::Utilities::intToString(i)
         + "_1";
      array_database = database->putDatabase(array_name);
      (d_data[i][0]).putUnregisteredToDatabase(array_database);

      array_name = "d_data%d_" + tbox::Utilities::intToString(i) + "_2";
      array_database = database->putDatabase(array_name);
      (d_data[i][1]).putUnregisteredToDatabase(array_database);
   }
}

}
}

#endif
