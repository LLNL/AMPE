/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_IndexData_C
#define included_pdat_IndexData_C

#include "SAMRAI/pdat/IndexData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

template<class TYPE, class BOX_GEOMETRY>
const int IndexData<TYPE, BOX_GEOMETRY>::PDAT_INDEXDATA_VERSION = 1;

template<class TYPE, class BOX_GEOMETRY>
IndexDataNode<TYPE, BOX_GEOMETRY>::IndexDataNode():
   d_index(tbox::Dimension::getInvalidDimension())
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexDataNode<TYPE, BOX_GEOMETRY>::IndexDataNode(
   const hier::Index& index,
   const int offset,
   TYPE& t,
   IndexDataNode<TYPE, BOX_GEOMETRY>* n,
   IndexDataNode<TYPE, BOX_GEOMETRY>* p):
   d_index(index),
   d_offset(offset),
   d_item(&t),
   d_next(n),
   d_prev(p)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexDataNode<TYPE, BOX_GEOMETRY>::~IndexDataNode()
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexDataNode<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::getNode()
{
   return *d_node;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   const IndexData<TYPE, BOX_GEOMETRY>& index_data,
   bool begin) :
   d_index_data(const_cast<IndexData<TYPE, BOX_GEOMETRY>*>(&index_data)),
   d_node(begin ? d_index_data->d_list_head : 0)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   IndexData<TYPE, BOX_GEOMETRY>* index_data,
   IndexDataNode<TYPE, BOX_GEOMETRY>* node) :
   d_index_data(index_data),
   d_node(node)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::IndexIterator(
   const IndexIterator<TYPE, BOX_GEOMETRY>& iter) :
   d_index_data(iter.d_index_data),
   d_node(iter.d_node)
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::~IndexIterator()
{
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::operator = (
   const IndexIterator<TYPE, BOX_GEOMETRY>& iter)
{
   d_index_data = iter.d_index_data;
   d_node = iter.d_node;
   return *this;
}

template<class TYPE, class BOX_GEOMETRY>
TYPE&
IndexIterator<TYPE, BOX_GEOMETRY>::operator * ()
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
const TYPE&
IndexIterator<TYPE, BOX_GEOMETRY>::operator * () const
{
   return *d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
const hier::Index&
IndexIterator<TYPE, BOX_GEOMETRY>::getIndex() const
{
   return d_node->d_index;
}

template<class TYPE, class BOX_GEOMETRY>
TYPE*
IndexIterator<TYPE, BOX_GEOMETRY>::operator -> ()
{
   return d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
const TYPE*
IndexIterator<TYPE, BOX_GEOMETRY>::operator -> () const
{
   return d_node->d_item;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::operator ++ ()
{
   if (d_node) {
      d_node = d_node->d_next;
   }
   return *this;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::operator ++ (
   int)
{
   IndexIterator<TYPE, BOX_GEOMETRY> tmp = *this;
   if (d_node) {
      d_node = d_node->d_next;
   }
   return tmp;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>&
IndexIterator<TYPE, BOX_GEOMETRY>::operator -- ()
{
   if (d_node) {
      d_node = d_node->d_prev;
   }
   return *this;
}

template<class TYPE, class BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>
IndexIterator<TYPE, BOX_GEOMETRY>::operator -- (
   int)
{
   IndexIterator<TYPE, BOX_GEOMETRY> tmp = *this;
   if (d_node) {
      d_node = d_node->d_prev;
   }
   return tmp;
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexIterator<TYPE, BOX_GEOMETRY>::operator == (
   const IndexIterator<TYPE, BOX_GEOMETRY>& i) const
{
   return d_node == i.d_node;
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexIterator<TYPE, BOX_GEOMETRY>::operator != (
   const IndexIterator<TYPE, BOX_GEOMETRY>& i) const
{
   return d_node != i.d_node;
}

/*
 *************************************************************************
 *
 * The constructor for the irregular grid object simply initializes the
 * irregular data list to be null (this is done implicitly).
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
IndexData<TYPE, BOX_GEOMETRY>::IndexData(
   const hier::Box& box,
   const hier::IntVector& ghosts):
   hier::PatchData(box, ghosts),
   d_dim(box.getDim()),
   d_data(hier::PatchData::getGhostBox().size()),
   d_list_head(NULL),
   d_list_tail(NULL),
   d_number_items(0)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
}

template<class TYPE, class BOX_GEOMETRY>
IndexData<TYPE, BOX_GEOMETRY>::~IndexData()
{
   removeAllItems();
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

template<class TYPE, class BOX_GEOMETRY>
IndexData<TYPE, BOX_GEOMETRY>::IndexData(
   const IndexData<TYPE, BOX_GEOMETRY>& foo):
   hier::PatchData(foo.getBox(), foo.getGhostCellWidth()),
   d_dim(foo.getDim())
{

   // private and not used (but included for some compilers)
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::operator = (
   const IndexData<TYPE, BOX_GEOMETRY>& foo)
{
   // private and not used (but included for some compilers)
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * Copy into dst where src overlaps on interiors.
 *
 *************************************************************************
 */
template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy(
   const hier::PatchData& src)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const IndexData<TYPE, BOX_GEOMETRY>* t_src =
      dynamic_cast<const IndexData<TYPE, BOX_GEOMETRY> *>(&src);

   TBOX_ASSERT(t_src != NULL);

   const hier::Box& src_ghost_box = t_src->getGhostBox();
   removeInsideBox(src_ghost_box);

   typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*t_src, false);
   for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*t_src, true);
        s != send;
        ++s) {
      if (getGhostBox().contains(s.getNode().d_index)) {
         appendItem(s.getNode().d_index, *(s.getNode().d_item));
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy2(
   hier::PatchData& dst) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   dst.copy(*this);
}

/*
 *************************************************************************
 *
 * Copy data from the source into the destination according to the
 * overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy(
   const hier::PatchData& src,
   const hier::BoxOverlap& overlap)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, src);

   const IndexData<TYPE, BOX_GEOMETRY>* t_src =
      dynamic_cast<const IndexData<TYPE, BOX_GEOMETRY> *>(&src);
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);

   TBOX_ASSERT(t_src != NULL);
   TBOX_ASSERT(t_overlap != NULL);

   const hier::IntVector& src_offset(t_overlap->getSourceOffset());
   const hier::BoxContainer& box_list = t_overlap->getDestinationBoxContainer();
   const hier::Box& src_ghost_box = t_src->getGhostBox();

   for (hier::BoxContainer::const_iterator b(box_list);
        b != box_list.end(); ++b) {
      const hier::Box& dst_box = *b;
      const hier::Box src_box(hier::Box::shift(*b, -src_offset));
      removeInsideBox(dst_box);
      typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*t_src, false);
      for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*t_src, true);
           s != send;
           ++s) {
         if (src_box.contains(s.getNode().d_index)) {
            TYPE new_item;
            new_item.copySourceItem(
               s.getNode().d_index,
               src_offset,
               *(t_src->d_data[src_ghost_box.offset(s.getNode().d_index)]->
                 d_item));
            appendItem(s.getNode().d_index + src_offset, new_item);
         }
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::copy2(
   hier::PatchData& dst,
   const hier::BoxOverlap& overlap) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, dst);

   dst.copy(*this, overlap);
}

/*
 *************************************************************************
 *
 * Calculate the buffer space needed to pack/unpack messages on the box
 * region using the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
bool
IndexData<TYPE, BOX_GEOMETRY>::canEstimateStreamSizeFromBox() const
{
   return false;
}

template<class TYPE, class BOX_GEOMETRY>
int
IndexData<TYPE, BOX_GEOMETRY>::getDataStreamSize(
   const hier::BoxOverlap& overlap) const
{
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(t_overlap != NULL);

   size_t bytes = 0;
   int num_items = 0;
   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b(boxes); b != boxes.end(); ++b) {
      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*b, -(t_overlap->getSourceOffset()));
      hier::Box::iterator indexend(box, false);
      for (hier::Box::iterator index(box, true); index != indexend; ++index) {
         TYPE* item = getItem(*index);
         if (item) {
            num_items++;
            bytes += item->getDataStreamSize();
         }
      }
   }
   const int index_size = d_dim.getValue() * tbox::MessageStream::getSizeof<int>();
   bytes += (num_items * index_size + tbox::MessageStream::getSizeof<int>());
   return static_cast<int>(bytes);
}

/*
 *************************************************************************
 *
 * Pack/unpack data into/out of the message streams using the index
 * space in the overlap descriptor.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap) const
{
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   int num_items = 0;
   for (hier::BoxContainer::const_iterator b(boxes); b != boxes.end(); ++b) {
      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*b, -(t_overlap->getSourceOffset()));
      typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*this, false);
      for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*this, true);
           s != send; ++s) {
         if (box.contains(s.getNode().d_index)) {
            num_items++;
         }
      }
   }

   stream << num_items;

   for (hier::BoxContainer::const_iterator c(boxes); c != boxes.end(); ++c) {
      hier::Box box = hier::PatchData::getBox()
         * hier::Box::shift(*c, -(t_overlap->getSourceOffset()));
      typename IndexData<TYPE, BOX_GEOMETRY>::iterator tend(*this, false);
      for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator t(*this, true);
           t != tend; ++t) {
         if (box.contains(t.getNode().d_index)) {
            TYPE* item = &(*t);
            TBOX_ASSERT(item != NULL);

            int index_buf[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
            for (int i = 0; i < d_dim.getValue(); i++) {
               index_buf[i] = t.getNode().d_index(i);
            }
            stream.pack(index_buf, d_dim.getValue());
            item->packStream(stream);
         }
      }
   }

}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxOverlap& overlap)
{
   const typename BOX_GEOMETRY::Overlap * t_overlap =
      dynamic_cast<const typename BOX_GEOMETRY::Overlap *>(&overlap);
   TBOX_ASSERT(t_overlap != NULL);

   int num_items;
   stream >> num_items;

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b(boxes); b != boxes.end(); ++b) {
      removeInsideBox(*b);
   }

   int i;
   TYPE* items = new TYPE[num_items];
   for (i = 0; i < num_items; i++) {
      int index_buf[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      stream.unpack(index_buf, d_dim.getValue());
      hier::Index index(d_dim);
      for (int j = 0; j < d_dim.getValue(); j++) {
         index(j) = index_buf[j];
      }
      (items + i)->unpackStream(stream, t_overlap->getSourceOffset());
      addItem(index + (t_overlap->getSourceOffset()), items[i]);
   }
   delete[] items;
}

/*
 *************************************************************************
 *
 * List manipulation stuff.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::appendItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }

   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != NULL);

   *new_item = item;
   addItemToList(index, offset, *new_item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::appendItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   appendItemToList(index, offset, *item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::addItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != NULL);

   *new_item = item;
   addItemToList(index, offset, *new_item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::addItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   if (isElement(offset)) {
      removeItem(offset);
   }
   addItemToList(index, offset, *item);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAddItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != NULL);

   *new_item = item;

   if (node == NULL) {

      addItemToList(index, offset, *new_item);

   } else {
      delete node->d_item;

      node->d_item = new_item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAddItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   if (node == NULL) {

      addItemToList(index, offset, *item);

   } else {

      delete node->d_item;

      node->d_item = item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAppendItem(
   const hier::Index& index,
   const TYPE& item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   TYPE* new_item = new TYPE();
   TBOX_ASSERT(new_item != NULL);

   *new_item = item;

   if (node == NULL) {

      appendItemToList(index, offset, *new_item);

   } else {
      delete node->d_item;

      node->d_item = new_item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::replaceAppendItemPointer(
   const hier::Index& index,
   TYPE* item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   if (node == NULL) {

      appendItemToList(index, offset, *item);

   } else {

      delete node->d_item;

      node->d_item = item;
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeItem(
   const hier::Index& index)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   int offset = hier::PatchData::getGhostBox().offset(index);
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   removeItem(offset);
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeItem(
   const int offset)
{
   TBOX_ASSERT(offset >= 0 && offset <= hier::PatchData::getGhostBox().size());

   IndexDataNode<TYPE, BOX_GEOMETRY>* node = d_data[offset];

   TBOX_ASSERT(node);

   removeNodeFromList(node);

   delete node->d_item;
   delete node;

   d_data[offset] = NULL;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::addItemToList(
   const hier::Index& index,
   const int offset,
   TYPE& item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   IndexDataNode<TYPE, BOX_GEOMETRY>* new_node =
      new IndexDataNode<TYPE, BOX_GEOMETRY>(index,
                                            offset,
                                            item,
                                            d_list_head,
                                            NULL);

   if (d_list_head) {
      d_list_head->d_prev = new_node;
   }

   d_list_head = new_node;

   if (!d_list_tail) {
      d_list_tail = new_node;
   }

   d_data[offset] = new_node;

   d_number_items++;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::appendItemToList(
   const hier::Index& index,
   const int offset,
   TYPE& item)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   IndexDataNode<TYPE, BOX_GEOMETRY>* new_node =
      new IndexDataNode<TYPE, BOX_GEOMETRY>(index,
                                            offset,
                                            item,
                                            NULL,
                                            d_list_tail);

   if (d_list_tail) {
      d_list_tail->d_next = new_node;
   }

   d_list_tail = new_node;

   if (!d_list_head) {
      d_list_head = new_node;
   }

   d_data[offset] = new_node;

   d_number_items++;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeNodeFromList(
   IndexDataNode<TYPE, BOX_GEOMETRY>* node)
{
   if ((d_list_head == node) && (d_list_tail == node)) {
      d_list_head = d_list_tail = NULL;

   } else if (d_list_head == node) {
      d_list_head = node->d_next;
      node->d_next->d_prev = NULL;

   } else if (d_list_tail == node) {
      d_list_tail = node->d_prev;
      node->d_prev->d_next = NULL;

   } else {
      node->d_next->d_prev = node->d_prev;
      node->d_prev->d_next = node->d_next;
   }

   d_data[node->d_offset] = NULL;

   d_number_items--;
}

template<class TYPE, class BOX_GEOMETRY>
int
IndexData<TYPE, BOX_GEOMETRY>::getNumberOfItems() const
{
   return d_number_items;
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeInsideBox(
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   typename IndexData<TYPE, BOX_GEOMETRY>::iterator l(*this, true);
   typename IndexData<TYPE, BOX_GEOMETRY>::iterator lend(*this, false);

   while (l != lend) {
      if (box.contains(l.getNode().d_index)) {
         hier::Index index(l.getNode().d_index);
         ++l;
         removeItem(index);
      } else {
         ++l;
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeOutsideBox(
   const hier::Box& box)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);

   typename IndexData<TYPE, BOX_GEOMETRY>::iterator l(*this, true);
   typename IndexData<TYPE, BOX_GEOMETRY>::iterator lend(*this, false);

   while (l != lend) {
      if (!box.contains(l.getNode().d_index)) {
         hier::Index index(l.getNode().d_index);
         ++l;
         removeItem(index);
      } else {
         ++l;
      }
   }
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeGhostItems()
{
   removeOutsideBox(hier::PatchData::getBox());
}

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::removeAllItems()
{
   removeInsideBox(hier::PatchData::getGhostBox());
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexData<TYPE, BOX_GEOMETRY>::isElement(
   const hier::Index& index) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);
   TBOX_ASSERT(hier::PatchData::getGhostBox().contains(index));

   return d_data[hier::PatchData::getGhostBox().offset(index)] != NULL;
}

template<class TYPE, class BOX_GEOMETRY>
bool
IndexData<TYPE, BOX_GEOMETRY>::isElement(
   int offset) const
{
   return d_data[offset] != NULL;
}

/*
 *************************************************************************
 *
 * Just checks to make sure that the class version is the same
 * as the restart file version number.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::getSpecializedFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("PDAT_INDEXDATA_VERSION");
   if (ver != PDAT_INDEXDATA_VERSION) {
      TBOX_ERROR("IndexData::getSpecializedFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   int item_count = 0;
   bool item_found = true;

   do {
      std::string index_keyword = "index_data_" + tbox::Utilities::intToString(
            item_count,
            6);

      if (database->isDatabase(index_keyword)) {

         boost::shared_ptr<tbox::Database> item_db(
            database->getDatabase(index_keyword));

         tbox::Array<int> index_array =
            item_db->getIntegerArray(index_keyword);
         hier::Index index(d_dim);
         for (int j = 0; j < d_dim.getValue(); j++) {
            index(j) = index_array[j];
         }

         TYPE item;
         item.getFromDatabase(item_db);

         appendItem(index, item);

      } else {
         item_found = false;
      }

      item_count++;

   } while (item_found);

}

/*
 *************************************************************************
 *
 * Just writes out the class version number to the database.
 *
 *************************************************************************
 */

template<class TYPE, class BOX_GEOMETRY>
void
IndexData<TYPE, BOX_GEOMETRY>::putSpecializedToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   TBOX_ASSERT(database);

   database->putInteger("PDAT_INDEXDATA_VERSION", PDAT_INDEXDATA_VERSION);

   int item_count = 0;
   typename IndexData<TYPE, BOX_GEOMETRY>::iterator send(*this, false);
   for (typename IndexData<TYPE, BOX_GEOMETRY>::iterator s(*this, true);
        s != send; ++s) {

      std::string index_keyword = "index_data_" + tbox::Utilities::intToString(
            item_count,
            6);
      hier::Index index = s.getNode().d_index;
      tbox::Array<int> index_array(d_dim.getValue());
      for (int i = 0; i < d_dim.getValue(); i++) {
         index_array[i] = index(i);
      }

      boost::shared_ptr<tbox::Database> item_db(
         database->putDatabase(index_keyword));

      item_db->putIntegerArray(index_keyword, index_array);

      TYPE* item = getItem(index);

      item->putUnregisteredToDatabase(item_db);

      item_count++;
   }
}

template<class TYPE, class BOX_GEOMETRY>
TYPE*
IndexData<TYPE, BOX_GEOMETRY>::getItem(
   const hier::Index& index) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, index);

   TYPE* item;
   if (!isElement(index)) {
      item = NULL;
   } else {
      item = d_data[hier::PatchData::getGhostBox().offset(index)]->d_item;
   }

   return item;
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
