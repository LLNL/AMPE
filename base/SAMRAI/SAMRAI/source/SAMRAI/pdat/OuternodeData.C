/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated outernode centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeData_C
#define included_pdat_OuternodeData_C

#include "SAMRAI/pdat/OuternodeData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int OuternodeData<TYPE>::PDAT_OUTERNODEDATA_VERSION = 1;

/*
 *************************************************************************
 *
 * Constructor and destructor for outernode data objects.  The
 * constructor simply initializes data variables and sets up the
 * array data.
 *
 *************************************************************************
 */

template<class TYPE>
OuternodeData<TYPE>::OuternodeData(
   const hier::Box& box,
   int depth):
   hier::PatchData(box, hier::IntVector::getZero(box.getDim())),
   d_depth(depth)
{
   TBOX_ASSERT(depth > 0);

   const tbox::Dimension& dim(box.getDim());

   for (int d = 0; d < dim.getValue(); d++) {

      hier::Box nodebox = NodeGeometry::toNodeBox(box);

      for (int dh = d + 1; dh < dim.getValue(); dh++) {

         /*
          * For dimensions higher than d, narrow the box down to avoid
          * representing edge and corner nodes multiple times.
          *
          *  i.e.    Y--Y--Y  outernodeX0 defined on nodes (0,1)
          *          |  |  |  outernodeX1 defined on nodes (2,1)
          *          X--o--X  outernodeY0 defined on node  (0,0)-(2,0)
          *          |  |  |  outernodeY1 defined on node  (0,2)-(2,2)
          *          Y--Y--Y
          *         node box
          */
         nodebox.lower(dh) += 1;
         nodebox.upper(dh) -= 1;
      }

      hier::Box outernodebox = nodebox;
      outernodebox.upper(d) = nodebox.lower(d);
      outernodebox.lower(d) = nodebox.lower(d);
      if (outernodebox.size() > 0) {
         d_data[d][0].initializeArray(outernodebox, depth);
      } else {
         d_data[d][0].invalidateArray(dim);
      }

      outernodebox = nodebox;
      outernodebox.lower(d) = nodebox.upper(d);
      outernodebox.upper(d) = nodebox.upper(d);
      d_data[d][1].initializeArray(outernodebox, depth);

   }
}

template<class TYPE>
OuternodeData<TYPE>::~OuternodeData()
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
OuternodeData<TYPE>::OuternodeData(
   const OuternodeData<TYPE>& foo):
   hier::PatchData(foo.getBox(), foo.getGhostCellWidth())

{
   NULL_USE(foo);
}

template<class TYPE>
void
OuternodeData<TYPE>::operator = (
   const OuternodeData<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
int
OuternodeData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
bool
OuternodeData<TYPE>::dataExists(
   int face_normal) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));

   return d_data[face_normal][0].isInitialized();
}

template<class TYPE>
TYPE*
OuternodeData<TYPE>::getPointer(
   int face_normal,
   int side,
   int depth)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal][side].getPointer(depth);
}

template<class TYPE>
const TYPE*
OuternodeData<TYPE>::getPointer(
   int face_normal,
   int side,
   int depth) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   return d_data[face_normal][side].getPointer(depth);
}

template<class TYPE>
ArrayData<TYPE>&
OuternodeData<TYPE>::getArrayData(
   int face_normal,
   int side)
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[face_normal][side];
}

template<class TYPE>
const ArrayData<TYPE>&
OuternodeData<TYPE>::getArrayData(
   int face_normal,
   int side) const
{
   TBOX_ASSERT((face_normal >= 0) && (face_normal < getDim().getValue()));
   TBOX_ASSERT((side == 0) || (side == 1));

   return d_data[face_normal][side];
}

template<class TYPE>
TYPE&
OuternodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth)
{
   for (int d = getDim().getValue() - 1; d >= 0; d--) {
      if (i[d] == d_data[d][0].getBox().lower()[d]) {
         return d_data[d][0](i, depth);
      }
      if (i[d] == d_data[d][1].getBox().upper()[d]) {
         return d_data[d][1](i, depth);
      }
   }

   /*
    * The following lines should only be executed if there's a bug
    * in the Outernode datatype.
    */
   TBOX_ERROR("Bad index used to access outernode data\n"
      << "Given index is not an outernode of this instance.\n");
   return d_data[0][0](i, depth);
}

template<class TYPE>
const TYPE&
OuternodeData<TYPE>::operator () (
   const NodeIndex& i,
   int depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   for (int d = getDim() - 1; d >= 0; d--) {
      if (i[d] == d_data[d][0].getBox().lower()[d]) {
         return d_data[d][0](i, depth);
      }
      if (i[d] == d_data[d][1].getBox().upper()[d]) {
         return d_data[d][1](i, depth);
      }
   }
   /*
    * The following lines should only be executed if there's a bug
    * in the Outernode datatype.
    */
   TBOX_ERROR("Bad index used to access outernode data\n"
      << "Given index is not an outernode of this instance.\n");
   return d_data[0][0](i, depth);
}

/*
 *************************************************************************
 *
 * Perform a fast copy between an outernode patch data type (source) and
 * a node patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const NodeData<TYPE> * const t_node_src =
      dynamic_cast<const NodeData<TYPE> *>(&src);
   const OuternodeData<TYPE> * const t_onode_src =
      dynamic_cast<const OuternodeData<TYPE> *>(&src);

   if (t_node_src != NULL) {
      copyFromNode(*t_node_src);
   } else if (t_onode_src != NULL) {
      copyFromOuternode(*t_onode_src);
   } else {
      TBOX_ERROR("OuternodeData<dim>::copy error!\n"
         << "Can copy only from NodeData<TYPE> or "
         << "OuternodeData<TYPE> of the same "
         << "dim and TYPE.");
   }

}

template<class TYPE>
void
OuternodeData<TYPE>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   NodeData<TYPE>* t_node_dst =
      dynamic_cast<NodeData<TYPE> *>(&dst);
   OuternodeData<TYPE>* t_onode_dst =
      dynamic_cast<OuternodeData<TYPE> *>(&dst);

   if (t_node_dst != NULL) {
      copyToNode(*t_node_dst);
   } else if (t_onode_dst != NULL) {
      copyToOuternode(*t_onode_dst);
   } else {
      TBOX_ERROR("OuternodeData<dim>::copy2 error!\n"
         << "Can copy only to NodeData<TYPE> or "
         << "OuternodeData<TYPE> of the same "
         << "dim and TYPE.");
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
OuternodeData<TYPE>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const NodeData<TYPE>* t_node_src =
      dynamic_cast<const NodeData<TYPE> *>(&src);
   const OuternodeData<TYPE>* t_onode_src =
      dynamic_cast<const OuternodeData<TYPE> *>(&src);

   if (t_node_src != NULL) {
      copyFromNode(*t_node_src, *t_overlap);
   } else if (t_onode_src != NULL) {
      copyFromOuternode(*t_onode_src, *t_overlap);
   } else {
      TBOX_ERROR("OuternodeData<dim>::copy error!\n"
         << "Can copy only from NodeData<TYPE> or "
         << "OuternodeData<TYPE> of the same "
         << "dim and TYPE.");
   }

}

template<class TYPE>
void
OuternodeData<TYPE>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   NodeData<TYPE>* t_node_dst =
      dynamic_cast<NodeData<TYPE> *>(&dst);
   OuternodeData<TYPE>* t_onode_dst =
      dynamic_cast<OuternodeData<TYPE> *>(&dst);

   if (t_node_dst != NULL) {
      copyToNode(*t_node_dst, *t_overlap);
   } else if (t_onode_dst != NULL) {
      copyToOuternode(*t_onode_dst, *t_overlap);
   } else {
      TBOX_ERROR("OuternodeData<dim>::copy2 error!\n"
         << "Can copy only to NodeData<TYPE> or "
         << "OuternodeData<TYPE> of the same "
         << "dim and TYPE.");
   }

}

/*
 *************************************************************************
 *
 * Perform a fast copy from a node data object to this outernode data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::copyDepth(
   int dst_depth,
   const NodeData<TYPE>& src,
   int src_depth)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const ArrayData<TYPE>& node_array = src.getArrayData();
   for (int d = 0; d < getDim().getValue(); d++) {
      for (int loc = 0; loc < 2; loc++) {
         ArrayData<TYPE>& onode_array = d_data[d][loc];
         onode_array.copyDepth(dst_depth,
            node_array,
            src_depth,
            onode_array.getBox());
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy to a node data object from this outernode data
 * object at the specified depths, where their index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::copyDepth2(
   int dst_depth,
   NodeData<TYPE>& dst,
   int src_depth) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   ArrayData<TYPE>& node_array = dst.getArrayData();
   for (int d = 0; d < getDim().getValue(); d++) {
      for (int loc = 0; loc < 2; loc++) {
         const ArrayData<TYPE>& onode_array = d_data[d][loc];
         node_array.copyDepth(dst_depth,
            onode_array,
            src_depth,
            onode_array.getBox());
      }
   }
}

/*
 *************************************************************************
 *
 * Add source data to the destination according to overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::sum(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const OuternodeData<TYPE>* t_onode_src =
      dynamic_cast<const OuternodeData<TYPE> *>(&src);

   // NOTE:  We assume this operation is only needed to
   //        copy and add data to another outernode data
   //        object.  If we ever need to provide this for node
   //        data or other flavors of the copy operation, we
   //        should refactor the routine similar to the way
   //        the regular copy operations are implemented.
   if (t_onode_src == NULL) {
      TBOX_ERROR("OuternodeData<dim>::sum error!\n"
         << "Can copy and add only from OuternodeData<TYPE> "
         << "of the same dim and TYPE.");
   } else {

      const hier::IntVector& src_offset = t_overlap->getSourceOffset();

      for (int src_d = 0; src_d < getDim().getValue(); src_d++) {
         for (int src_p = 0; src_p < 2; src_p++) {

            const ArrayData<TYPE>& src_array =
               t_onode_src->d_data[src_d][src_p];
            const hier::BoxContainer& box_list =
               t_overlap->getDestinationBoxContainer();

            for (int dst_d = 0; dst_d < getDim().getValue(); dst_d++) {
               for (int dst_p = 0; dst_p < 2; dst_p++) {
                  if (d_data[dst_d][dst_p].isInitialized()) {
                     d_data[dst_d][dst_p].sum(
                        src_array, box_list, src_offset);
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
 * Calculate the buffer space needed to pack/unpack messages on the box
 * region using the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
bool
OuternodeData<TYPE>::canEstimateStreamSizeFromBox() const
{
   return ArrayData<TYPE>::canEstimateStreamSizeFromBox();
}

template<class TYPE>
int
OuternodeData<TYPE>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   int size = 0;
   const hier::BoxContainer& boxlist = t_overlap->getDestinationBoxContainer();
   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
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
OuternodeData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxContainer& dst_boxes = t_overlap->getDestinationBoxContainer();
   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (hier::BoxContainer::const_iterator dst_box(dst_boxes);
        dst_box != dst_boxes.end(); ++dst_box) {
      const hier::Box src_box =
         hier::Box::shift(*dst_box, -src_offset);
      for (int d = 0; d < getDim().getValue(); d++) {
         for (int loc = 0; loc < 2; loc++) {
            const hier::Box intersect = src_box * d_data[d][loc].getBox();
            if (!intersect.empty()) {
               const hier::Box pack_box =
                  hier::Box::shift(intersect, src_offset);
               d_data[d][loc].packStream(stream, pack_box, src_offset);
            }
         }
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxContainer& dst_boxes = t_overlap->getDestinationBoxContainer();
   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (hier::BoxContainer::const_iterator dst_box(dst_boxes);
        dst_box != dst_boxes.end(); ++dst_box) {
      for (int d = 0; d < getDim().getValue(); d++) {
         for (int f = 0; f < 2; f++) {
            const hier::Box intersect =
               (*dst_box) * d_data[d][f].getBox();
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
 * Unpack data from the message stream and add to this outernode data
 * object using the index space in the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::unpackStreamAndSum(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const NodeOverlap* t_overlap =
      dynamic_cast<const NodeOverlap *>(&overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxContainer& dst_boxes = t_overlap->getDestinationBoxContainer();
   const hier::IntVector& src_offset = t_overlap->getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      for (hier::BoxContainer::const_iterator dst_box(dst_boxes);
           dst_box != dst_boxes.end(); ++dst_box) {
         for (int f = 0; f < 2; f++) {
            const hier::Box intersect =
               (*dst_box) * d_data[d][f].getBox();
            if (!intersect.empty()) {
               d_data[d][f].unpackStreamAndSum(stream, intersect, src_offset);
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Calculate the amount of memory space needed to represent the data
 * for a  outernode centered grid.
 *
 *************************************************************************
 */

template<class TYPE>
size_t
OuternodeData<TYPE>::getSizeOfData(
   const hier::Box& box,
   int depth)
{
   TBOX_ASSERT(depth > 0);

   size_t size = 0;
   for (int d = 0; d < box.getDim().getValue(); d++) {
      hier::Box loc0 = NodeGeometry::toNodeBox(box);
      hier::Box loc1 = NodeGeometry::toNodeBox(box);
      loc0.upper(d) = box.lower(d);
      loc0.lower(d) = box.lower(d);
      loc1.lower(d) = box.upper(d);
      loc1.upper(d) = box.upper(d);

      for (int dh = d + 1; dh < box.getDim().getValue(); dh++) {

         /*
          * For dimensions higher than d, narrow the box down to avoid
          * representing edge and corner nodes multiple times.
          */
         loc0.lower(dh) += 1;
         loc0.upper(dh) -= 1;
         loc1.lower(dh) += 1;
         loc1.upper(dh) -= 1;
      }
      size += ArrayData<TYPE>::getSizeOfData(loc0, depth)
         + ArrayData<TYPE>::getSizeOfData(loc1, depth);
   }
   return size;
}

/*
 *************************************************************************
 *
 * Compute the box of valid node indices given values of
 * dimension and side designating the set of data indices.
 *
 *************************************************************************
 */

template<class TYPE>
hier::Box
OuternodeData<TYPE>::getDataBox(
   int face_normal,
   int side)
{
   if (face_normal < 0 || face_normal >= getDim().getValue() || side < 0 || side > 1) {
      TBOX_ERROR("Bad values for face_normal and/or side in\n"
         "OuternodeData<dim>::getDataBox().\n");
   }

   /*
    * We start with the full box and chop it down to the databox
    * corresponding to the given face_normal and side.
    */
   hier::Box databox = NodeGeometry::toNodeBox(getBox());
   const hier::IntVector& ghosts = getGhostCellWidth();

   for (int dh = face_normal + 1; dh < getDim().getValue(); dh++) {

      /*
       * For dimensions higher than d, narrow the box down to avoid
       * representing edge and corner nodes multiple times.
       */
      databox.lower(dh) += 1;
      databox.upper(dh) -= 1;
   }

   if (side == 0) {
      databox.upper(face_normal) = databox.lower(face_normal);
      databox.lower(face_normal) = databox.lower(face_normal)
         - ghosts(face_normal);
   } else { // side == 1
      databox.lower(face_normal) = databox.upper(face_normal);
      databox.upper(face_normal) = databox.upper(face_normal)
         + ghosts(face_normal);
   }
   return databox;
}

/*
 *************************************************************************
 *
 * Fill the outernode centered box with the given value.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::fill(
   const TYPE& t,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fill(t, d);
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fill(t, d);
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   int d)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fill(t, NodeGeometry::toNodeBox(box), d);
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fill(t, NodeGeometry::toNodeBox(box), d);
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::fillAll(
   const TYPE& t)
{
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fillAll(t);
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fillAll(t);
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_data[i][0].isInitialized()) {
         d_data[i][0].fillAll(t, NodeGeometry::toNodeBox(box));
      }
      if (d_data[i][1].isInitialized()) {
         d_data[i][1].fillAll(t, NodeGeometry::toNodeBox(box));
      }
   }
}

/*
 *************************************************************************
 *
 * Perform a fast copy between an outernode patch data type (source) and
 * a node patch data type (destination) where the index spaces overlap.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::copyFromNode(
   const NodeData<TYPE>& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const ArrayData<TYPE>& node_array = src.getArrayData();
   for (int d = 0; d < getDim().getValue(); d++) {
      for (int loc = 0; loc < 2; loc++) {
         ArrayData<TYPE>& onode_array = d_data[d][loc];
         if (onode_array.isInitialized()) {
            onode_array.copy(node_array, onode_array.getBox());
         }
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::copyToNode(
   NodeData<TYPE>& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   ArrayData<TYPE>& node_array = dst.getArrayData();
   for (int d = 0; d < getDim().getValue(); d++) {
      for (int loc = 0; loc < 2; loc++) {
         if (d_data[d][loc].isInitialized()) {
            node_array.copy(d_data[d][loc], d_data[d][loc].getBox());
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
OuternodeData<TYPE>::copyFromNode(
   const NodeData<TYPE>& src,
   const NodeOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const hier::IntVector& src_offset = overlap.getSourceOffset();
   for (int d = 0; d < getDim().getValue(); d++) {
      const hier::BoxContainer& box_list = overlap.getDestinationBoxContainer();
      if (d_data[d][0].isInitialized()) {
         d_data[d][0].copy(src.getArrayData(), box_list, src_offset);
      }
      if (d_data[d][1].isInitialized()) {
         d_data[d][1].copy(src.getArrayData(), box_list, src_offset);
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::copyToNode(
   NodeData<TYPE>& dst,
   const NodeOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   const hier::IntVector& src_offset = overlap.getSourceOffset();
   const hier::BoxContainer& box_list = overlap.getDestinationBoxContainer();
   for (int d = 0; d < getDim().getValue(); d++) {
      if (d_data[d][0].isInitialized()) {
         dst.getArrayData().copy(d_data[d][0], box_list, src_offset);
      }
      if (d_data[d][1].isInitialized()) {
         dst.getArrayData().copy(d_data[d][1], box_list, src_offset);
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::copyFromOuternode(
   const OuternodeData<TYPE>& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   for (int src_d = 0; src_d < getDim().getValue(); src_d++) {
      for (int src_p = 0; src_p < 2; src_p++) {

         const ArrayData<TYPE>& src_array = src.d_data[src_d][src_p];

         for (int dst_d = 0; dst_d < getDim().getValue(); dst_d++) {
            for (int dst_p = 0; dst_p < 2; dst_p++) {
               ArrayData<TYPE>& onode_array = d_data[dst_d][dst_p];
               if (onode_array.isInitialized()) {
                  onode_array.copy(src_array, onode_array.getBox());
               }
            }
         }
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::copyFromOuternode(
   const OuternodeData<TYPE>& src,
   const NodeOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const hier::IntVector& src_offset = overlap.getSourceOffset();
   for (int src_d = 0; src_d < getDim().getValue(); src_d++) {
      for (int src_p = 0; src_p < 2; src_p++) {

         const ArrayData<TYPE>& src_array = src.d_data[src_d][src_p];
         const hier::BoxContainer& box_list = overlap.getDestinationBoxContainer();

         for (int dst_d = 0; dst_d < getDim().getValue(); dst_d++) {
            for (int dst_p = 0; dst_p < 2; dst_p++) {
               if (d_data[dst_d][dst_p].isInitialized()) {
                  d_data[dst_d][dst_p].copy(src_array, box_list, src_offset);
               }
            }
         }
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::copyToOuternode(
   OuternodeData<TYPE>& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   for (int dst_d = 0; dst_d < getDim().getValue(); dst_d++) {
      for (int dst_p = 0; dst_p < 2; dst_p++) {

         ArrayData<TYPE>& dst_array = dst.d_data[dst_d][dst_p];

         for (int src_d = 0; src_d < getDim().getValue(); src_d++) {
            for (int src_p = 0; src_p < 2; src_p++) {
               if (d_data[src_d][src_p].isInitialized()) {
                  dst_array.copy(d_data[src_d][src_p],
                     d_data[src_d][src_p].getBox());
               }
            }
         }
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::copyToOuternode(
   OuternodeData<TYPE>& dst,
   const NodeOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   const hier::IntVector& src_offset = overlap.getSourceOffset();
   const hier::BoxContainer& box_list = overlap.getDestinationBoxContainer();

   for (int dst_d = 0; dst_d < getDim().getValue(); dst_d++) {
      for (int dst_p = 0; dst_p < 2; dst_p++) {

         ArrayData<TYPE>& dst_array = dst.d_data[dst_d][dst_p];
         for (int src_d = 0; src_d < getDim().getValue(); src_d++) {
            for (int src_p = 0; src_p < 2; src_p++) {
               if (d_data[src_d][src_p].isInitialized()) {
                  dst_array.copy(d_data[src_d][src_p], box_list, src_offset);
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Print routines for outernode centered arrays.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::print(
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
OuternodeData<TYPE>::print(
   const hier::Box& box,
   int depth,
   std::ostream& os,
   int prec) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
   TBOX_ASSERT((depth >= 0) && (depth < d_depth));

   for (int axis = 0; axis < getDim().getValue(); axis++) {
      os << "Array axis = " << axis << std::endl;
      for (int side = 0; side < 2; side++) {
         os << "Side = " << ((side == 0) ? "lower" : "upper") << std::endl;
         printAxisSide(axis, side, box, depth, os, prec);
      }
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::printAxisSide(
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
      printAxisSide(face_normal, side, box, d, os, prec);
   }
}

template<class TYPE>
void
OuternodeData<TYPE>::printAxisSide(
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

   const hier::Box nodebox = NodeGeometry::toNodeBox(box);
   const hier::Box region = nodebox * d_data[face_normal][side].getBox();
   os.precision(prec);
   hier::Box::iterator iend(region, false);
   for (hier::Box::iterator i(region, true); i != iend; ++i) {
      os << "array" << *i << " = "
         << d_data[face_normal][side](*i, depth) << std::endl << std::flush;
   }
}

/*
 *************************************************************************
 *
 * Checks that class version and restart file version are equal.
 * If so, reads in d_depth from the database.
 * Then has each item in d_data read in its data from the database.
 *
 *************************************************************************
 */

template<class TYPE>
void
OuternodeData<TYPE>::getSpecializedFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("PDAT_OUTERNODEDATA_VERSION");
   if (ver != PDAT_OUTERNODEDATA_VERSION) {
      TBOX_ERROR("OuternodeData<dim>::getSpecializedFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");

   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      std::string array_name = "d_data" + tbox::Utilities::intToString(i)
         + "_1";
      if (database->keyExists(array_name)) {
         array_database = database->getDatabase(array_name);
         (d_data[i][0]).getFromDatabase(array_database);
      }

      array_name = "d_data" + tbox::Utilities::intToString(i) + "_2";
      if (database->keyExists(array_name)) {
         array_database = database->getDatabase(array_name);
         (d_data[i][1]).getFromDatabase(array_database);
      }
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
OuternodeData<TYPE>::putSpecializedToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   TBOX_ASSERT(database);

   database->putInteger("PDAT_OUTERNODEDATA_VERSION",
      PDAT_OUTERNODEDATA_VERSION);

   database->putInteger("d_depth", d_depth);

   std::string array_name;
   boost::shared_ptr<tbox::Database> array_database;
   for (int i = 0; i < getDim().getValue(); i++) {
      if (d_data[i][0].isInitialized()) {
         array_name = "d_data" + tbox::Utilities::intToString(i) + "_1";
         array_database = database->putDatabase(array_name);
         (d_data[i][0]).putUnregisteredToDatabase(array_database);
      }
      if (d_data[i][1].isInitialized()) {
         array_name = "d_data" + tbox::Utilities::intToString(i) + "_2";
         array_database = database->putDatabase(array_name);
         (d_data[i][1]).putUnregisteredToDatabase(array_database);
      }
   }
}

}
}

#endif
