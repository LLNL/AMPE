/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
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
   int depth)
{
   return tbox::MemoryUtilities::align(box.size() * depth * sizeof(TYPE));
}

/*
 *************************************************************************
 *
 * Default constructor creates an object that absolutely should
 * not be used until it is initialized using initializeArray().
 *
 *************************************************************************
 */
template<class TYPE>
ArrayData<TYPE>::ArrayData():
   d_dim(tbox::Dimension::getInvalidDimension()),
   d_depth(0),
   d_offset(0),
   d_array(0)
{
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
   int depth):
   d_dim(box.getDim()),
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

/*
 *************************************************************************
 *
 * The const constructor and assignment operator are not actually used
 * but are defined here for compilers that require an implementation for
 * every declaration.
 *
 *************************************************************************
 */

template<class TYPE>
ArrayData<TYPE>::ArrayData(
   const ArrayData<TYPE>& foo):
   d_dim(foo.d_dim),
   d_depth(0),
   d_offset(0),
   d_box(d_dim),
   d_array(0)
{
   /*
    * Throw an assertion if this method is used.
    */
   TBOX_ASSERT(true);
   NULL_USE(foo);
}

template<class TYPE>
void
ArrayData<TYPE>::operator = (
   const ArrayData<TYPE>& foo)
{
   /*
    * Throw an assertion if this method is used.
    */
   TBOX_ASSERT(true);
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * Initialize the array using the specified box, depth.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::initializeArray(
   const hier::Box& box,
   int depth)
{
   TBOX_ASSERT(depth > 0);

   d_dim = box.getDim();
   d_depth = depth;
   d_offset = box.size();
   d_box = box;

   /*
    * If array is a different size then create a new one, otherwise
    * just use existing array; there is no reason to
    * deallocate/allocate if existing is the correct size.
    *
    * Also note that the array is not initialized to anything.  This is
    * to avoid the performance hit of initializing but is potentially
    * bad.  In debug mode will explicitly initialize the array.
    *
    * Performance timing of some sample applications (LinAdv) showed
    * the impact of initialization was significant.  This is probably
    * not the safest programming practice.
    */
   if (d_array.size() != depth * d_offset) {
      d_array = tbox::Array<TYPE>(depth * d_offset,
                                  tbox::Array<TYPE>::UNINITIALIZED);
   }

#ifdef DEBUG_INITIALIZE_UNDEFINED
   undefineData();
#endif
}

template<class TYPE>
bool
ArrayData<TYPE>::isInitialized() const
{
   return d_depth > 0;
}

template<class TYPE>
const hier::Box&
ArrayData<TYPE>::getBox() const
{
   return d_box;
}

template<class TYPE>
int
ArrayData<TYPE>::getDepth() const
{
   return d_depth;
}

template<class TYPE>
int
ArrayData<TYPE>::getOffset() const
{
   return d_offset;
}

template<class TYPE>
TYPE*
ArrayData<TYPE>::getPointer(
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_array.getPointer(d * d_offset);
}

template<class TYPE>
const TYPE*
ArrayData<TYPE>::getPointer(
   int d) const
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   return d_array.getPointer(d * d_offset);
}

template<class TYPE>
TYPE&
ArrayData<TYPE>::operator () (
   const hier::Index& i,
   int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   const int index = d_box.offset(i) + d * d_offset;

   TBOX_ASSERT((index >= 0) && (index < d_depth * d_offset));

   return d_array[index];
}

template<class TYPE>
const TYPE&
ArrayData<TYPE>::operator () (
   const hier::Index& i,
   int d) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, i);

   TBOX_ASSERT((d >= 0) && (d < d_depth));

   const int index = d_box.offset(i) + d * d_offset;

   TBOX_ASSERT((index >= 0) && (index < d_depth * d_offset));

   return d_array[index];
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
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, src, box);

   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ((d_depth == src.d_depth) &&
       (d_box.isSpatiallyEqual(src.d_box)) &&
       (box.isSpatiallyEqual(d_box))) {

      TYPE * const dst_ptr = d_array.getPointer();
      const TYPE * const src_ptr = src.d_array.getPointer();
      const int n = d_offset * d_depth;
      for (int i = 0; i < n; i++) {
         copyop(dst_ptr[i], src_ptr[i]);
      }

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector src_shift(d_dim, 0);

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

   if (src_shift == hier::IntVector::getZero(d_dim)) {

      copy(src, box);

   } else {

      const hier::Box copybox =
         box * d_box * hier::Box::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

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
       && transformation.getOffset() == hier::IntVector::getZero(d_dim)) {

      copy(src, box);

   } else {

      hier::Box transformed_src(src.d_box);
      transformation.transform(transformed_src);
      const hier::Box copybox(
         box * d_box * transformed_src);

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

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
   for (hier::BoxContainer::const_iterator b(boxes); b != boxes.end(); ++b) {
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
   for (hier::BoxContainer::const_iterator b(boxes); b != boxes.end(); ++b) {
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
   int dst_depth,
   const ArrayData<TYPE>& src,
   int src_depth,
   const hier::Box& box)
{
   TBOX_ASSERT((0 <= dst_depth) && (dst_depth <= d_depth));
   TBOX_ASSERT((0 <= src_depth) && (src_depth <= src.d_depth));

   CopyOperation<TYPE> copyop;

   /*
    * Do a fast copy of data if all data aligns with copy region
    */

   if ((d_box.isSpatiallyEqual(src.d_box)) && (box.isSpatiallyEqual(d_box))) {

      TYPE * const dst_ptr = d_array.getPointer();
      const TYPE * const src_ptr = src.d_array.getPointer();

      TYPE * const dst_ptr_d = dst_ptr + dst_depth * d_offset;
      const TYPE * const src_ptr_d = src_ptr + src_depth * d_offset;
      for (int i = 0; i < d_offset; i++) {
         copyop(dst_ptr_d[i], src_ptr_d[i]);
      }

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const int dst_start_depth = dst_depth;
         const int src_start_depth = src_depth;
         const int num_depth = 1;
         const hier::IntVector src_shift(d_dim, 0);

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

      TYPE * const dst_ptr = d_array.getPointer();
      const TYPE * const src_ptr = src.d_array.getPointer();
      const int n = d_offset * d_depth;
      for (int i = 0; i < n; i++) {
         sumop(dst_ptr[i], src_ptr[i]);
      }

   } else {

      const hier::Box copybox = box * d_box * src.d_box;

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
         const hier::IntVector src_shift(d_dim, 0);

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

   if (src_shift == hier::IntVector::getZero(d_dim)) {

      sum(src, box);

   } else {

      const hier::Box copybox =
         box * d_box * hier::Box::shift(src.d_box, src_shift);

      if (!copybox.empty()) {

         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

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
   for (hier::BoxContainer::const_iterator b(boxes); b != boxes.end(); ++b) {
      sum(src, *b, src_shift);
   }
}

template<class TYPE>
int
ArrayData<TYPE>::getDataStreamSize(
   const hier::BoxContainer& boxes,
   const hier::IntVector& source_shift) const
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(source_shift);
#endif

   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, source_shift);

   const int nelements = boxes.getTotalSizeOfBoxes();

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

   const int size = d_depth * dest_box.size();
   tbox::Array<TYPE> buffer(size);

   packBuffer(buffer.getPointer(),
      hier::Box::shift(dest_box, -src_shift));

   stream.pack(buffer.getPointer(), size);

}

template<class TYPE>
void
ArrayData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::IntVector& src_shift) const
{

   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Array<TYPE> buffer(size);

   int ptr = 0;
   for (hier::BoxContainer::const_iterator b(dest_boxes);
        b != dest_boxes.end(); ++b) {
      packBuffer(buffer.getPointer(ptr),
         hier::Box::shift(*b, -src_shift));
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);

   stream.pack(buffer.getPointer(), size);

}

template<class TYPE>
void
ArrayData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::Box& dest_box,
   const hier::Transformation& transformation) const
{

   const int size = d_depth * dest_box.size();
   tbox::Array<TYPE> buffer(size);

   hier::Box pack_box(dest_box);
   transformation.inverseTransform(pack_box);
   packBuffer(buffer.getPointer(), pack_box);
//      hier::Box::shift(dest_box, -src_shift));

   stream.pack(buffer.getPointer(), size);

}


template<class TYPE>
void
ArrayData<TYPE>::packStream(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::Transformation& transformation) const
{

   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Array<TYPE> buffer(size);

   int ptr = 0;
   for (hier::BoxContainer::const_iterator b(dest_boxes);
        b != dest_boxes.end(); ++b) {
      hier::Box pack_box(*b);
      transformation.inverseTransform(pack_box);
      packBuffer(buffer.getPointer(ptr), pack_box);
//         hier::Box::shift(*b, -src_shift));
      ptr += d_depth * b->size();
   }

   TBOX_ASSERT(ptr == size);

   stream.pack(buffer.getPointer(), size);

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

   const int size = d_depth * dest_box.size();
   tbox::Array<TYPE> buffer(size);

   stream.unpack(buffer.getPointer(), size);
   unpackBuffer(buffer.getPointer(), dest_box);

}

template<class TYPE>
void
ArrayData<TYPE>::unpackStream(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::IntVector& src_shift)
{

   NULL_USE(src_shift);

   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Array<TYPE> buffer(size);

   stream.unpack(buffer.getPointer(), size);

   int ptr = 0;
   for (hier::BoxContainer::const_iterator b(dest_boxes);
        b != dest_boxes.end(); ++b) {
      unpackBuffer(buffer.getPointer(ptr), *b);
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

   const int size = d_depth * dest_box.size();
   tbox::Array<TYPE> buffer(size);

   stream.unpack(buffer.getPointer(), size);
   unpackBufferAndSum(buffer.getPointer(), dest_box);

}

template<class TYPE>
void
ArrayData<TYPE>::unpackStreamAndSum(
   tbox::MessageStream& stream,
   const hier::BoxContainer& dest_boxes,
   const hier::IntVector& src_shift)
{

   NULL_USE(src_shift);

   const int size = d_depth * dest_boxes.getTotalSizeOfBoxes();
   tbox::Array<TYPE> buffer(size);

   stream.unpack(buffer.getPointer(), size);

   int ptr = 0;
   for (hier::BoxContainer::const_iterator b(dest_boxes);
        b != dest_boxes.end(); ++b) {
      unpackBufferAndSum(buffer.getPointer(ptr), *b);
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
      TYPE* ptr = d_array.getPointer();
      const int n = d_depth * d_offset;
      for (int i = 0; i < n; i++) {
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
   for (int d = 0; d < d_depth; d++) {
      fill(t, box, d);
   }
}

template<class TYPE>
void
ArrayData<TYPE>::fill(
   const TYPE& t,
   const int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   if (!d_box.empty()) {
      TYPE* ptr = d_array.getPointer(d * d_offset);
      const int n = d_offset;
      for (int i = 0; i < n; i++) {
         ptr[i] = t;
      }
   }
}

template<class TYPE>
void
ArrayData<TYPE>::fill(
   const TYPE& t,
   const hier::Box& box,
   const int d)
{
   TBOX_ASSERT((d >= 0) && (d < d_depth));

   const hier::Box ispace = d_box * box;

   if (!ispace.empty()) {

      int box_w[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      int dst_w[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      int dim_counter[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      for (int i = 0; i < d_dim.getValue(); i++) {
         box_w[i] = ispace.numberCells(i);
         dst_w[i] = d_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int num_d0_blocks = ispace.size() / box_w[0];

      int dst_counter = d_box.offset(ispace.lower()) + d * d_offset;

      int dst_b[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      for (int nd = 0; nd < d_dim.getValue(); nd++) {
         dst_b[nd] = dst_counter;
      }

      TYPE * const dst_ptr = d_array.getPointer();

      for (int nb = 0; nb < num_d0_blocks; nb++) {

         for (int i0 = 0; i0 < box_w[0]; i0++) {
            dst_ptr[dst_counter + i0] = t;
         }
         int dim_jump = 0;

         for (int j = 1; j < d_dim.getValue(); j++) {
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
            for (int k = 0; k < dim_jump; k++) {
               dst_step *= dst_w[k];
            }
            dst_counter = dst_b[dim_jump - 1] + dst_step;

            for (int m = 0; m < dim_jump; m++) {
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
 * database.  Then calls getSpecializedFromDatabase() to read in the
 * actual data.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::getFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   TBOX_ASSERT(database);

   int ver = database->getInteger("PDAT_ARRAYDATA_VERSION");
   if (ver != PDAT_ARRAYDATA_VERSION) {
      TBOX_ERROR("ArrayData::getFromDatabase error...\n"
         << " : Restart file version different than class version" << std::endl);
   }

   d_depth = database->getInteger("d_depth");
   d_offset = database->getInteger("d_offset");
   d_box = database->getDatabaseBox("d_box");

   getSpecializedFromDatabase(database);
}

/*
 *************************************************************************
 *
 * Write out the class version number, d_depth, d_offset, and d_box
 * to the database.  Then calls putSpecializedToDatabase() to write
 * in the actual data.
 *
 *************************************************************************
 */

template<class TYPE>
void
ArrayData<TYPE>::putUnregisteredToDatabase(
   const boost::shared_ptr<tbox::Database>& database,
   bool data_only) const
{
   TBOX_ASSERT(database);

   if (!data_only) {
      database->putInteger("PDAT_ARRAYDATA_VERSION", PDAT_ARRAYDATA_VERSION);

      database->putInteger("d_depth", d_depth);
      database->putInteger("d_offset", d_offset);
      database->putDatabaseBox("d_box", d_box);
   }

   putSpecializedToDatabase(database);
}

template<class TYPE>
void
ArrayData<TYPE>::putSpecializedToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   database->putArray("d_array", d_array);
}

template<class TYPE>
void
ArrayData<TYPE>::getSpecializedFromDatabase(
   const boost::shared_ptr<tbox::Database>& database)
{
   database->getArray("d_array", d_array);
}

template<class TYPE>
const tbox::Dimension&
ArrayData<TYPE>::getDim() const
{
   return d_dim;
}

template<class TYPE>
void
ArrayData<TYPE>::invalidateArray(
   const tbox::Dimension& dim) {
   d_dim = dim;
   d_depth = 0;
   d_offset = 0;
   d_box = hier::Box::getEmptyBox(dim);
}

template<class TYPE>
bool
ArrayData<TYPE>::isValid()
{
   return !d_box.isEmpty();
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
