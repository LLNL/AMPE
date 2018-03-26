/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated node centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_NodeData_C
#define included_pdat_NodeData_C

#include "SAMRAI/pdat/NodeData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int NodeData<TYPE>::PDAT_NODEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for node data objects.  The constructor
 * simply initializes data variables and sets up the array data.
 *
 *************************************************************************
 */

template<class TYPE>
NodeData<TYPE>::NodeData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_depth(depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
   TBOX_ASSERT(depth > 0);
   TBOX_ASSERT(ghosts.min() >= 0);

   const hier::Box node = NodeGeometry::toNodeBox(getGhostBox());
   d_data.initializeArray(node, depth);

}

template<class TYPE>
NodeData<TYPE>::~NodeData()
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
NodeData<TYPE>::NodeData(
   const NodeData<TYPE>& foo):
   hier::PatchData(foo.getBox(), foo.getGhostCellWidth())
{
   NULL_USE(foo);
}

template<class TYPE>
void
NodeData<TYPE>::operator = (
   const NodeData<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
int
NodeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
ArrayData<TYPE>&
NodeData<TYPE>::getArrayData()
{
   return d_data;
}

template<class TYPE>
const ArrayData<TYPE>&
NodeData<TYPE>::getArrayData() const
{
   return d_data;
}

template<class TYPE>
TYPE*
NodeData<TYPE>::getPointer(
   int depth)
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data.getPointer(depth);
}

template<class TYPE>
const TYPE*
NodeData<TYPE>::getPointer(
   int depth) const
{
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data.getPointer(depth);
}

template<class TYPE>
TYPE&
NodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data(i, depth);
}

template<class TYPE>
const TYPE&
NodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data(i, depth);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between two node centered arrays where their
 * index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const NodeData<TYPE>* t_src =
      dynamic_cast<const NodeData<TYPE> *>(&src);
   if (t_src == NULL) {
      src.copy2(*this);
   } else {
      const hier::Box box = d_data.getBox() * t_src->d_data.getBox();
      if (!box.empty()) {
         d_data.copy(t_src->d_data, box);
      }
   }
}

template<class TYPE>
void
NodeData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   NodeData<TYPE>* t_dst =
      dynamic_cast<NodeData<TYPE> *>(&dst);

   TBOX_ASSERT(t_dst != NULL);

   const hier::Box box = d_data.getBox() * t_dst->d_data.getBox();
   if (!box.empty()) {
      t_dst->d_data.copy(d_data, box);
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
NodeData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const NodeData<TYPE>* t_src =
      dynamic_cast<const NodeData<TYPE> *>(&src);
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   if ((t_src == NULL) || (t_overlap == NULL)) {
      src.copy2(*this, overlap);
   } else {
      if (t_overlap->getTransformation().getRotation() ==
          hier::Transformation::NO_ROTATE) {
         d_data.copy(t_src->d_data,
            t_overlap->getDestinationBoxContainer(),
            t_overlap->getTransformation());
      } else {
         copyWithRotation(*t_src, *t_overlap);
      }
   }
}

template<class TYPE>
void
NodeData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   NodeData<TYPE>* t_dst =
      dynamic_cast<NodeData<TYPE> *>(&dst);
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_dst != NULL);
   TBOX_ASSERT(t_overlap != NULL);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {
      t_dst->d_data.copy(d_data,
         t_overlap->getDestinationBoxContainer(),
         t_overlap->getTransformation());
   } else {
      t_dst->copyWithRotation(*this, *t_overlap);
   }
}

template<class TYPE>
void
NodeData<TYPE>::copyOnBox(
   const NodeData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);
   const hier::Box node_box = NodeGeometry::toNodeBox(box);
   d_data.copy(src.getArrayData(), node_box);
}

template<class TYPE>
void
NodeData<TYPE>::copyWithRotation(
   const NodeData<TYPE>& src,
   const NodeOverlap& overlap)
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(src.getDim());
   const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer();
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   hier::Box rotatebox(src.getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   hier::Box node_rotatebox(NodeGeometry::toNodeBox(rotatebox));

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   node_rotatebox.getBlockId(),
                                   getBox().getBlockId());

   for (hier::BoxContainer::const_iterator bi(overlap_boxes);
        bi != overlap_boxes.end(); ++bi) {
      const hier::Box& overlap_box = *bi;

      const hier::Box copybox(node_rotatebox * overlap_box);

      if (!copybox.empty()) {
         const int depth = ((getDepth() < src.getDepth()) ?
                            getDepth() : src.getDepth());

         hier::Box::iterator ciend(copybox, false);
         for (hier::Box::iterator ci(copybox, true); ci != ciend; ++ci) {

            NodeIndex dst_index(*ci, hier::IntVector::getZero(dim));
            NodeIndex src_index(dst_index);
            NodeGeometry::transform(src_index, back_trans);

            for (int d = 0; d < depth; d++) {
               d_data(dst_index, d) = src.d_data(src_index, d);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy from a node data object to this node data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::copyDepth(
   int dst_depth,
   const NodeData<TYPE>& src,
   int src_depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const hier::Box box = d_data.getBox() * src.d_data.getBox();
   if (!box.empty()) {
      d_data.copyDepth(dst_depth, src.d_data, src_depth, box);
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
NodeData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
int
NodeData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   return d_data.getDataStreamSize(t_overlap->getDestinationBoxContainer(),
      t_overlap->getSourceOffset());
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
NodeData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   if (t_overlap->getTransformation().getRotation() ==
       hier::Transformation::NO_ROTATE) {
      d_data.packStream(stream,
         t_overlap->getDestinationBoxContainer(),
         t_overlap->getTransformation());
   } else {
      packWithRotation(stream, *t_overlap);
   }
}

template<class TYPE>
void
NodeData<TYPE>::packWithRotation(
   tbox::MessageStream& stream,
   const NodeOverlap& overlap) const
{
   TBOX_ASSERT(overlap.getTransformation().getRotation() !=
      hier::Transformation::NO_ROTATE);

   const tbox::Dimension& dim(getDim());
   const hier::BoxContainer& overlap_boxes = overlap.getDestinationBoxContainer();
   const hier::Transformation::RotationIdentifier rotate =
      overlap.getTransformation().getRotation();
   const hier::IntVector& shift = overlap.getSourceOffset();

   hier::Box rotatebox(getGhostBox());
   overlap.getTransformation().transform(rotatebox);

   hier::Box node_rotatebox(NodeGeometry::toNodeBox(rotatebox));

   const hier::Transformation::RotationIdentifier back_rotate =
      hier::Transformation::getReverseRotationIdentifier(
         rotate, dim);

   hier::IntVector back_shift(dim);

   hier::Transformation::calculateReverseShift(
      back_shift, shift, rotate);

   hier::Transformation back_trans(back_rotate, back_shift,
                                   rotatebox.getBlockId(),
                                   getBox().getBlockId());

   const int depth = getDepth();

   const int size = depth * overlap_boxes.getTotalSizeOfBoxes();
   tbox::Array<TYPE> buffer(size);

   int i = 0;
   for (hier::BoxContainer::const_iterator bi(overlap_boxes);
        bi != overlap_boxes.end(); ++bi) {
      const hier::Box& overlap_box = *bi;

      const hier::Box copybox(node_rotatebox * overlap_box);

      if (!copybox.empty()) {

         for (int d = 0; d < depth; d++) {
            hier::Box::iterator ciend(copybox, false);
            for (hier::Box::iterator ci(copybox, true); ci != ciend; ++ci) {

               NodeIndex src_index(*ci, hier::IntVector::getZero(dim));
               NodeGeometry::transform(src_index, back_trans);

               buffer[i] = d_data(src_index, d);
               i++;
            }
         }
      }
   }

   stream.pack(buffer.getPointer(), size);
}

template<class TYPE>
void
NodeData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   d_data.unpackStream(stream,
      t_overlap->getDestinationBoxContainer(),
      t_overlap->getSourceOffset());
}

template<class TYPE>
void
NodeData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data.fill(t, d);
}

template<class TYPE>
void
NodeData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   d_data.fill(t, NodeGeometry::toNodeBox(box), d);
}

template<class TYPE>
void
NodeData<TYPE>::fillAll(
   const TYPE& t)
{
   d_data.fillAll(t);
}

template<class TYPE>
void
NodeData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   d_data.fillAll(t, NodeGeometry::toNodeBox(box));
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a node centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
NodeData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth,
   const hier::IntVector& ghosts)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
   TBOX_ASSERT(depth > 0);

   const hier::Box ghost_box = hier::Box::grow(box, ghosts);
   const hier::Box node_box = NodeGeometry::toNodeBox(ghost_box);
   return ArrayData<TYPE>::getSizeOfData(node_box, depth);
}

/*
 *************************************************************************
 *
 * Print node-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::print(
   const hier::Box& box,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int d = 0; d < d_depth; d++) {
      os << "Array depth = " << d << std::endl;
      print(box, d, os, prec);
   }
}

template<class TYPE>
void
NodeData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   os.precision(prec);
   NodeIterator iend(box, false);
   for (NodeIterator i(box, true); i != iend; ++i) {
      os << "array" << *i << " = "
         << d_data(*i, depth) << std::endl << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks to make sure that the class version and restart file
 * version are equal.  If so, reads in d_depth and has d_data
 * retrieve its own data from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::getSpecializedFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("PDAT_NODEDATA_VERSION");
   if (ver != PDAT_NODEDATA_VERSION) {
      TBOX_ERROR("NodeData<DIM>::getSpecializedFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   d_data.getFromDatabase(database->getDatabase("d_data"));
}

/*
 *************************************************************************
 *
 * Writes out the class version number and d_depth, Then has d_data
 * write its own data to the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
NodeData<TYPE>::putSpecializedToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{

   TBOX_ASSERT(database);

   database->putInteger("PDAT_NODEDATA_VERSION", PDAT_NODEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   d_data.putUnregisteredToDatabase(database->putDatabase("d_data"));
}

}
}

#endif
