/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Templated array data structure supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayData_C
#define included_pdat_ArrayData_C

#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/ArrayDataOperationUtilities.h"
#include "SAMRAI/pdat/CopyOperation.h"
#include "SAMRAI/pdat/SumOperation.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace pdat {

template<class TYPE>
const int ArrayData<TYPE>::PDAT_ARRAYDATA_VERSION = 1;

template<class TYPE>
bool
ArrayData<TYPE>::canEstimateStreamSizeFromBox()
{
   if ((typeid(TYPE) == typeid(bool))
       | (typeid(TYPE) == typeid(char))
       | (typeid(TYPE) == typeid(double))
       | (typeid(TYPE) == typeid(float))
       | (typeid(TYPE) == typeid(int))
       | (typeid(TYPE) == typeid(dcomplex))) {
      return true;
   } else {
      return false;
   }
}

template<class TYPE>
size_t
ArrayData<TYPE>::getSizeOfData(
   const hier::Box& box,
   unsigned int depth)
{
   return tbox::MemoryUtilities::align(box.size() * depth * sizeof(TYPE));
}

/*
 *************************************************************************
 *
 * The main constructor allocates data for the given box and depth.  It
 * does not initialize the memory.  The destructor automatically
 * deallocates memory via the array destructor.
 *
 *************************************************************************
 */

template<class TYPE>
ArrayData<TYPE>::ArrayData(
   const hier::Box& box,
   unsigned int depth):
   d_depth(depth),
   d_offset(box.size()),
   d_box(box),
   d_array(d_depth * d_offset)
{
   TBOX_ASSERT(depth > 0);

#ifdef DEBUG_INITIALIZE_UNDEFINED
   undefineData();
#endif
}

template<class TYPE>
ArrayData<TYPE>::~ArrayData()
{
}

template<class TYPE>
bool
ArrayData<TYPE>::isInitialized() const
{
   return d_depth * d_offset > 0;
}

template<class TYPE>
const hier::Box&
ArrayData<TYPE>::getBox() const
{
   return d_box;
}

template<class TYPE>
unsigned int
ArrayData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
size_t
ArrayData<TYPE>::getOffset() const
{
   return d_offset;
}

template<class TYPE>
size_t
ArrayData<TYPE>::getIndex(
   const hier::Index& i,
   unsigned int d) const
{
   TBOX_ASSERT((d < d_depth));

   size_t index = d_box.offset(i) + d * d_offset;

   TBOX_ASSERT((index < d_depth * d_offset));

   return index;
}

template<class TYPE>
TYPE *
ArrayData<TYPE>::getPointer(
   unsigned int d)
{
   TBOX_ASSERT((d < d_depth));

   return &d_array[d * d_offset];
}

template<class TYPE>
const TYPE *
ArrayData<TYPE>::getPointer(
   unsigned int d) const
{
   TBOX_ASSERT((d < d_depth));

   return &d_array[d * d_offset];
}

template<class TYPE>
TYPE&
ArrayData<TYPE>::operator () (
   const hier::Index& i,
   unsigned int d)
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((d < d_depth));

   return d_array[getIndex(i, d)];
}

template<class TYPE>
const TYPE&
ArrayData<TYPE>::operator () (
   const hier::Index& i,
   unsigned int d) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, i);
   TBOX_ASSERT((d < d_depth));

   return d_array[getIndex(i, d)];
}

/*
 *************************************************************************
 *
 * Copy data between two array data objects on a specified box domain.
 * Don't use C++ indexing member functions, since compilers are probably
 * too stupid to do strength reduction on the loops to get performance.
 *
 * If the source box, destination box, and copy box are the same and the
 * source and destination have the same depth, then perform a fast copy
 * of all data.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::copy(
   const ArrayData<TYPE>& src,
   const hier::Box& box)
{
   TBOX_ASSERT_OBJDIM_EQUALITY3(*this, src, box);

   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ((d_depth == src.d_depth) &&
       (d_box.isSpatiallyEqual(src.d_box)) &&
       (box.isSpatiallyEqual(d_box))) {

      TYPE * const dst_ptr = &d_array[0];
      const TYPE * const src_ptr = &src.d_array[0];
      const size_t n = d_offset * d_depth;
      for (size_t i = 0; i < n; ++i) {
         copyop(dst_ptr[i], src_ptr[i]);
      }

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector src_shift(box.getDim(), 0);

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
         doArrayDataOperationOnBox(*this,
            src,
            copybox,
            src_shift,
            dst_start_depth,
            src_start_depth,
            num_depth,
            copyop);

      }

   }

}

/*
 *************************************************************************
 *
 * Copy data from source ArrayData object to this (destination)
 * ArrayData object on given box domain.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::copy(
   const ArrayData<TYPE>& src,
   const hier::Box& box,
   const hier::IntVector& src_shift)
{

   if (src_shift == hier::IntVector::getZero(box.getDim())) {

      copy(src, box);

   } else {

      const hier::Box copybox =
         box * d_box * hier::Box::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         CopyOperation<TYPE> copyop;

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
         doArrayDataOperationOnBox(*this,
            src,
            copybox,
            src_shift,
            dst_start_depth,
            src_start_depth,
            num_depth,
            copyop);

      }

   }

}

template<class TYPE>
void
ArrayData<TYPE>::copy(
   const ArrayData<TYPE>& src,
   const hier::Box& box,
   const hier::Transformation& transformation)
{
   if (transformation.getRotation() == hier::Transformation::NO_ROTATE
       && transformation.getOffset() == hier::IntVector::getZero(box.getDim())
       && transformation.getBeginBlock() == transformation.getEndBlock()) {

      copy(src, box);

   } else {

      hier::Box transformed_src(src.d_box);
      transformation.transform(transformed_src);
      const hier::Box copybox(
         box * d_box * transformed_src);

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         CopyOperation<TYPE> copyop;

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
         doArrayDataOperationOnBox(*this,
            src,
            copybox,
            transformation.getOffset(),
            dst_start_depth,
            src_start_depth,
            num_depth,
            copyop);

      }

   }

}

/*
 *************************************************************************
 *
 * Copy over the boxlist by calling the single-box copy for each box in
 * the boxlist.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::copy(
   const ArrayData<TYPE>& src,
   const hier::BoxContainer& boxes,
   const hier::IntVector& src_shift)
{
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      copy(src, *b, src_shift);
   }
}

template<class TYPE>
void
ArrayData<TYPE>::copy(
   const ArrayData<TYPE>& src,
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation)
{
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      copy(src, *b, transformation);
   }
}

/*
 *************************************************************************
 *
 * Copy data between two array data objects on a specified box domain.
 *
 * If the source box, destination box, and copy box are the same and the
 * source and destination have the same depth, then perform a fast copy
 * of all data.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::copyDepth(
   unsigned int dst_depth,
   const ArrayData<TYPE>& src,
   unsigned int src_depth,
   const hier::Box& box)
{
   TBOX_ASSERT((dst_depth <= d_depth));
   TBOX_ASSERT((src_depth <= src.d_depth));

   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ((d_box.isSpatiallyEqual(src.d_box)) && (box.isSpatiallyEqual(d_box))) {

      TYPE * const dst_ptr = &d_array[0];
      const TYPE * const src_ptr = &src.d_array[0];

      TYPE * const dst_ptr_d = dst_ptr + dst_depth * d_offset;
      const TYPE * const src_ptr_d = src_ptr + src_depth * d_offset;
      for (size_t i = 0; i < d_offset; ++i) {
         copyop(dst_ptr_d[i], src_ptr_d[i]);
      }

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = dst_depth;
         const unsigned int src_start_depth = src_depth;
         const unsigned int num_depth = 1;
         const hier::IntVector src_shift(box.getDim(), 0);

         ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
         doArrayDataOperationOnBox(*this,
            src,
            copybox,
            src_shift,
            dst_start_depth,
            src_start_depth,
            num_depth,
            copyop);

      }

   }

}

/*
 *************************************************************************
 *
 * Add data from source ArrayData object to this (destination)
 * ArrayData object on given box region.
 *
 * If the source box, destination box, and copy box are the same and the
 * source and destination have the same depth, then perform a fast sum
 * on all data rather than performing explicit looping operations.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::sum(
   const ArrayData<TYPE>& src,
   const hier::Box& box)
{

   SumOperation<TYPE> sumop;

   /*
    * Do a fast copy and add if all data aligns with copy region
    */

   if ((d_depth == src.d_depth) &&
       (d_box.isSpatiallyEqual(src.d_box)) &&
       (box.isSpatiallyEqual(d_box))) {

      TYPE * const dst_ptr = &d_array[0];
      const TYPE * const src_ptr = &src.d_array[0];
      const size_t n = d_offset * d_depth;
      for (size_t i = 0; i < n; ++i) {
         sumop(dst_ptr[i], src_ptr[i]);
      }

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector src_shift(box.getDim(), 0);

         ArrayDataOperationUtilities<TYPE, SumOperation<TYPE> >::
         doArrayDataOperationOnBox(*this,
            src,
            copybox,
            src_shift,
            dst_start_depth,
            src_start_depth,
            num_depth,
            sumop);

      }

   }

}

/*
 *************************************************************************
 *
 * Add data from source ArrayData object to this (destination)
 * ArrayData object on region described by given box and offset.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::sum(
   const ArrayData<TYPE>& src,
   const hier::Box& box,
   const hier::IntVector& src_shift)
{

   if (src_shift == hier::IntVector::getZero(box.getDim())) {

      sum(src, box);

   } else {

      const hier::Box copybox =
         box * d_box * hier::Box::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const unsigned int dst_start_depth = 0;
         const unsigned int src_start_depth = 0;
         const unsigned int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         SumOperation<TYPE> sumop;

         ArrayDataOperationUtilities<TYPE, SumOperation<TYPE> >::
         doArrayDataOperationOnBox(*this,
            src,
            copybox,
            src_shift,
            dst_start_depth,
            src_start_depth,
            num_depth,
            sumop);

      }

   }

}

/*
 *************************************************************************
 *
 * Add data from source ArrayData object to this (destination)
 * ArrayData object on regions described by given boxes and offset.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::sum(
   const ArrayData<TYPE>& src,
   const hier::BoxContainer& boxes,
   const hier::IntVector& src_shift)
{
   for (hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {
      sum(src, *b, src_shift);
   }
}

template<class TYPE>
size_t
ArrayData<TYPE>::getDataStreamSize(
   const hier::BoxContainer& boxes,
   const hier::IntVector& source_shift) const
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(source_shift);
#endif

   TBOX_ASSERT_OBJDIM_EQUALITY2(*this, source_shift);

   const size_t nelements = boxes.getTotalSizeOfBoxes();

   if (typeid(TYPE) == typeid(bool)) {
      return tbox::MessageStream::getSizeof<bool>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(char)) {
      return tbox::MessageStream::getSizeof<char>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(dcomplex)) {
      return tbox::MessageStream::getSizeof<dcomplex>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(double)) {
      return tbox::MessageStream::getSizeof<double>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(float)) {
      return tbox::MessageStream::getSizeof<float>(d_depth * nelements);
   } else if (typeid(TYPE) == typeid(int)) {
      return tbox::MessageStream::getSizeof<int>(d_depth * nelements);
   }

   TBOX_ERROR("ArrayData::getDataStreamSize() -- Invalid type" << std::endl);
   return 0;
}

/*
 *************************************************************************
 *
 * Pack data into the message stream.  Both packing routines add one
 * level of copy into a temporary buffer to reduce the number of calls
 * to the abstract stream packing routines.  These definitions will only
 * work for the standard built-in types of bool, char, double, float,
 * and int.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::Box& dest_box,
   const hier::IntVector& src_shift) const
{

   const size_t size = d_depth * dest_box.size();
   std::vector<TYPE> buffer(size);

   packBuffer(&buffer[0], hier::Box::shift(dest_box, -src_shift));

   stream.pack(&buffer[0], size);

}

template<class TYPE>
void
ArrayData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::IntVector& src_shift) const
{

   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   std::vector<TYPE> buffer(size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      packBuffer(&buffer[ptr], hier::Box::shift(*b, -src_shift));
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);

   stream.pack(&buffer[0], size);

}

template<class TYPE>
void
ArrayData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::Box& dest_box,
   const hier::Transformation& transformation) const
{

   const size_t size = d_depth * dest_box.size();
   std::vector<TYPE> buffer(size);

   hier::Box pack_box(dest_box);
   transformation.inverseTransform(pack_box);
   packBuffer(&buffer[0], pack_box);
//      hier::Box::shift(dest_box, -src_shift));

   stream.pack(&buffer[0], size);

}

template<class TYPE>
void
ArrayData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::Transformation& transformation) const
{

   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   std::vector<TYPE> buffer(size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      hier::Box pack_box(*b);
      transformation.inverseTransform(pack_box);
      packBuffer(&buffer[ptr], pack_box);
//         hier::Box::shift(*b, -src_shift));
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);

   stream.pack(&buffer[0], size);

}

/*
 *************************************************************************
 *
 * Unpack data from the message stream.  Both unpacking routines add one
 * level of copy into a temporary buffer to reduce the number of calls
 * to the abstract stream packing routines.  These definitions will only
 * work for the standard built-in types of bool, char, double, float,
 * and int.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::Box& dest_box,
   const hier::IntVector& src_shift)
{

   NULL_USE(src_shift);

   const size_t size = d_depth * dest_box.size();
   std::vector<TYPE> buffer(size);

   stream.unpack(&buffer[0], size);
   unpackBuffer(&buffer[0], dest_box);

}

template<class TYPE>
void
ArrayData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::IntVector& src_shift)
{

   NULL_USE(src_shift);

   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   std::vector<TYPE> buffer(size);

   stream.unpack(&buffer[0], size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      unpackBuffer(&buffer[ptr], *b);
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);
}

/*
 *************************************************************************
 *
 * Unpack data from the message stream and add to this array data object.
 * Both unpacking routines add one level of copy into a temporary buffer
 * to reduce the number of calls to the abstract stream packing routines.
 * These definitions will only work for the standard built-in types of
 * bool, char, double, float, and int.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::unpackStreamAndSum(
   tbox::MessageStream& stream,
   const hier::Box& dest_box,
   const hier::IntVector& src_shift)
{

   NULL_USE(src_shift);

   const size_t size = d_depth * dest_box.size();
   std::vector<TYPE> buffer(size);

   stream.unpack(&buffer[0], size);
   unpackBufferAndSum(&buffer[0], dest_box);

}

template<class TYPE>
void
ArrayData<TYPE>::unpackStreamAndSum(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::IntVector& src_shift)
{

   NULL_USE(src_shift);

   const size_t size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   std::vector<TYPE> buffer(size);

   stream.unpack(&buffer[0], size);

   size_t ptr = 0;
   for (hier::BoxContainer::const_iterator b = dest_boxes.begin();
        b != dest_boxes.end(); ++b) {
      unpackBufferAndSum(&buffer[ptr], *b);
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);
}

/*
 *************************************************************************
 *
 * Fill all or portions of the array with the specified data value.
 * The templated TYPE must define the assignment operator.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::fillAll(
   const TYPE& t)
{
   if (!d_box.empty()) {
      TYPE* ptr = &d_array[0];
      const size_t n = d_depth * d_offset;
      for (size_t i = 0; i < n; ++i) {
         ptr[i] = t;
      }
   }
}

template<class TYPE>
void
ArrayData<TYPE>::fillAll(
   const TYPE& t,
   const hier::Box& box)
{
   for (tbox::Dimension::dir_t d = 0; d < d_depth; ++d) {
      fill(t, box, d);
   }
}

template<class TYPE>
void
ArrayData<TYPE>::fill(
   const TYPE& t,
   const unsigned int d)
{
   TBOX_ASSERT((d < d_depth));

   if (!d_box.empty()) {
      TYPE* ptr = &d_array[d * d_offset];
      const size_t n = d_offset;
      for (size_t i = 0; i < n; ++i) {
         ptr[i] = t;
      }
   }
}

template<class TYPE>
void
ArrayData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   const unsigned int d)
{
   TBOX_ASSERT((d < d_depth));

   const hier::Box ispace = d_box * box;

   if (!ispace.empty()) {

      const tbox::Dimension& dim = box.getDim();

      int box_w[SAMRAI::MAX_DIM_VAL];
      int dst_w[SAMRAI::MAX_DIM_VAL];
      int dim_counter[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t i = 0; i < dim.getValue(); ++i) {
         box_w[i] = ispace.numberCells(i);
         dst_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int num_d0_blocks = static_cast<int>(ispace.size() / box_w[0]);

      size_t dst_counter = d_box.offset(ispace.lower()) + d * d_offset;

      size_t dst_b[SAMRAI::MAX_DIM_VAL];
      for (tbox::Dimension::dir_t nd = 0; nd < dim.getValue(); ++nd) {
         dst_b[nd] = dst_counter;
      }

      TYPE * const dst_ptr = &d_array[0];

      for (int nb = 0; nb < num_d0_blocks; ++nb) {

         for (int i0 = 0; i0 < box_w[0]; ++i0) {
            dst_ptr[dst_counter + i0] = t;
         }
         int dim_jump = 0;

         for (int j = 1; j < dim.getValue(); ++j) {
            if (dim_counter[j] < box_w[j] - 1) {
               ++dim_counter[j];
               dim_jump = j;
               break;
            } else {
               dim_counter[j] = 0;
            }
         }

         if (dim_jump > 0) {
            int dst_step = 1;
            for (int k = 0; k < dim_jump; ++k) {
               dst_step *= dst_w[k];
            }
            dst_counter = dst_b[dim_jump - 1] + dst_step;

            for (int m = 0; m < dim_jump; ++m) {
               dst_b[m] = dst_counter;
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Checks to make sure that class and restart file version numbers are
 * equal.  If so, reads in d_depth, d_offset, and d_box from the
 * database.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::getFromRestart(
   const boost::shared_ptr<tbox::Database>& restart_db)
{
   TBOX_ASSERT(restart_db);

   int ver = restart_db->getInteger("PDAT_ARRAYDATA_VERSION");
   if (ver != PDAT_ARRAYDATA_VERSION) {
      TBOX_ERROR("ArrayData::getFromRestart error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = restart_db->getInteger("d_depth");
   d_offset = restart_db->getInteger("d_offset");
   d_box = restart_db->getDatabaseBox("d_box");

   restart_db->getVector("d_array", d_array);
}

/*
 *************************************************************************
 *
 * Write out the class version number, d_depth, d_offset, d_box, and
 * d_array to the restart database.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
   TBOX_ASSERT(restart_db);

   restart_db->putInteger("PDAT_ARRAYDATA_VERSION", PDAT_ARRAYDATA_VERSION);

   restart_db->putInteger("d_depth", d_depth);
   restart_db->putInteger("d_offset", static_cast<int>(d_offset));
   restart_db->putDatabaseBox("d_box", d_box);

   restart_db->putVector("d_array", d_array);
}

template<class TYPE>
const tbox::Dimension&
ArrayData<TYPE>::getDim() const
{
   return d_box.getDim();
}

template<class TYPE>
bool
ArrayData<TYPE>::isValid()
{
   return !d_box.empty();
}

/*
 *************************************************************************
 *
 * Set all array data to undefined values.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::undefineData()
{
   fillAll(tbox::MathUtilities<TYPE>::getSignalingNaN());
}

/*
 *************************************************************************
 *
 * Private member functions to pack and unpack data on the specified box
 * (for all components) into/from the buffer.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::packBuffer(
   TYPE* buffer,
   const hier::Box& box) const
{
   TBOX_ASSERT((box * d_box).isSpatiallyEqual(box));

   bool src_is_buffer = false;

   CopyOperation<TYPE> copyop;

   ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
   doArrayDataBufferOperationOnBox(*this,
      buffer,
      box,
      src_is_buffer,
      copyop);

}

template<class TYPE>
void
ArrayData<TYPE>::unpackBuffer(
   const TYPE* buffer,
   const hier::Box& box)
{
   TBOX_ASSERT((box * d_box).isSpatiallyEqual(box));

   bool src_is_buffer = true;

   CopyOperation<TYPE> copyop;

   ArrayDataOperationUtilities<TYPE, CopyOperation<TYPE> >::
   doArrayDataBufferOperationOnBox(*this,
      buffer,
      box,
      src_is_buffer,
      copyop);
}

/*
 *************************************************************************
 *
 * Private member function to unpack data on the specified box
 * (all components) from the buffer and add to this array data object.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::unpackBufferAndSum(
   const TYPE* buffer,
   const hier::Box& box)
{
   TBOX_ASSERT((box * d_box).isSpatiallyEqual(box));

   bool src_is_buffer = true;

   SumOperation<TYPE> sumop;

   ArrayDataOperationUtilities<TYPE, SumOperation<TYPE> >::
   doArrayDataBufferOperationOnBox(*this,
      buffer,
      box,
      src_is_buffer,
      sumop);
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
