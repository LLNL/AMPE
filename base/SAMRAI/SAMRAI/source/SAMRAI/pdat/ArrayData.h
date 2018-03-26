/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated array data structure supporting patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayData
#define included_pdat_ArrayData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayDataIterator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/tbox/MessageStream.h"

#include <typeinfo>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class ArrayData<TYPE> is a basic templated array structure defined
 * over the index space of a box (with a specified depth) that provides
 * the support for the various standard array-based patch data subclasses.
 *
 * The data storage is in (i,...,k,d) order, where i,...,k indicates
 * spatial indices and the d indicates the component at that location.
 * Memory allocation is in column-major ordering (e.g., Fortran style)
 * so that the leftmost index runs fastest in memory.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.  Note that a number of
 * functions only work for standard built-in types (bool, char, double,
 * float, and int).  To use this class with other user-defined types,
 * many of these functions will need to be specialized, especially those
 * that deal with message packing and unpacking.
 */

template<class TYPE>
class ArrayData
{
public:
   /*!
    * Static member function that returns tru when the amount of buffer space in a
    * message stream can be estimated from box only.  For built-in types (bool, char,
    * double, float, int, and dcomplex), this routine returns true.  For other
    * data types (template paramters) that may require special handling,
    * a different implementation must be provided.
    */
   static bool
   canEstimateStreamSizeFromBox();

   /*!
    * Static member function that returns the amount of memory space needed to
    * store data of given depth on a box.
    *
    * Note that this function is only defined for the standard data types:
    * bool, char, double, float, int, and dcomplex.  It must be provided for other
    * template parameter types.
    *
    * @return size_t value indicating the amount of memory space needed for the data.
    *
    * @param box   Const reference to box object describing the spatial extents
    *              of the array data index region of interest.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    */
   static size_t
   getSizeOfData(
      const hier::Box& box,
      int depth);

   /*!
    * The default constructor creates an empty array data object.
    * The initializeArray() member function must be called before the
    * array can be used.
    */
   ArrayData();

   /*!
    * Construct an array data object.
    *
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space associated with the array data object.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    */
   ArrayData(
      const hier::Box& box,
      int depth);

   /*!
    * The destructor for an array data object releases all memory allocated
    * for the array elements.
    */
   ~ArrayData();

   /*!
    * Initialize the size of array data and depth.  This method is
    * somewhat poorly named as the data is NOT initialized to anything
    * to avoid the performance cost of the data initialization.
    *
    * Use undefineData to initialize the data.
    *
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space associated with the array data object.
    * @param depth Integer number of data values at each spatial location in
    *              the array.
    */
   void
   initializeArray(
      const hier::Box& box,
      int depth);

   /*!
    * @brief Returns true when the array has been properly initialized
    * and storage has been allocated; otherwise, return false.
    *
    * Note: Only arrays that have been initialized can do anything useful.
    * Initialize an uninitialized array by calling the initializeArray() method.
    */
   bool
   isInitialized() const;

   /*!
    * Set the array data to an ``undefined'' state appropriate for the data type.
    * For example, for float and double, this means setting data to signaling NaNs
    * that cause a floating point assertion when used in a numerical expression
    * without being set to valid values.
    */
   void
   undefineData();

   /*!
    * Return the box over which the array is defined.
    */
   const hier::Box&
   getBox() const;

   /*!
    * Return the depth (e.g., the number of data values at each spatial
    * location) of this array.
    */
   int
   getDepth() const;

   /*!
    * Return the offset (e.g., the number of data values for each
    * depth component) of this array.
    */
   int
   getOffset() const;

   /*!
    * Get a non-const pointer to the beginning of the given depth
    * component of this data array.
    */
   TYPE *
   getPointer(
      const int d = 0);

   /*!
    * Get a const pointer to the beginning of the given depth
    * component of this data array.
    */
   const TYPE *
   getPointer(
      const int d = 0) const;

   /*!
    * Return reference to value in this array associated with the given
    * box index and depth component.
    */
   TYPE&
   operator () (
      const hier::Index& i,
      const int d);

   /*!
    * Return const reference to value in this array associated with the given
    * box index and depth component.
    */
   const TYPE&
   operator () (
      const hier::Index& i,
      const int d) const;

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the copy operation.
    *              Note: the box is in either the source or destination index space
    *                    (which are assumed to be the same).
    */
   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::Box& box);

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified box with the destination box and
    * shifted source box to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the copy operation.
    *              Note: the box is in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::Box& box,
      const hier::IntVector& src_shift);

   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::Box& box,
      const hier::Transformation& transformation);

   /*!
    * Copy data from the source array data object to this array data object
    * on the specified index space regions.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified boxes with the destination box and
    * shifted source box to find the regions of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions over which to perform the copy operation.
    *              Note: the boxes are in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::BoxContainer& boxes,
      const hier::IntVector& src_shift);

   void
   copy(
      const ArrayData<TYPE>& src,
      const hier::BoxContainer& boxes,
      const hier::Transformation& transformation);

   /*!
    * Copy given source depth of source array data object to given destination
    * depth of this array data object on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param dst_depth Integer depth of destination array.
    * @param src       Const reference to source array data object.
    * @param src_depth Integer depth of source array.
    * @param box       Const reference to box object describing the spatial extents
    *                  of the index space region over which to perform the copy operation.
    *                  Note: the box is in either the source or destination index space
    *                        (which are assumed to be the same).
    */
   void
   copyDepth(
      int dst_depth,
      const ArrayData<TYPE>& src,
      int src_depth,
      const hier::Box& box);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source and destination
    * box regions require no shifting to make them consistent.  This routine
    * will intersect the specified box with the source and destination boxes
    * to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the sum operation.
    *              Note: the box is in either the source or destination index space
    *                    (which are assumed to be the same).
    */
   void
   sum(
      const ArrayData<TYPE>& src,
      const hier::Box& box);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space region.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified box with the destination box and
    * shifted source box to find the region of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param box   Const reference to box object describing the spatial extents
    *              of the index space region over which to perform the sum operation.
    *              Note: the box is in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void
   sum(
      const ArrayData<TYPE>& src,
      const hier::Box& box,
      const hier::IntVector& src_shift);

   /*!
    * Add data from the source array data object to this array data object
    * on the specified index space regions.
    *
    * Note that this routine assumes that the source array box region must
    * be shifted to be consistent with the destination (this) array box region.
    * This routine will intersect the specified boxes with the destination box and
    * shifted source box to find the regions of intersection.
    *
    * @param src   Const reference to source array data object.
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions over which to perform the sum operation.
    *              Note: the boxes are in the destination index space.
    * @param src_shift Const reference to shift vector used to put the source
    *              array data box into the index space region of this array data object.
    */
   void
   sum(
      const ArrayData<TYPE>& src,
      const hier::BoxContainer& boxes,
      const hier::IntVector& src_shift);

   /*!
    * Calculate the number of bytes needed to stream the data living
    * in the specified box domains.  This routine is only defined for
    * the built-in types of bool, char, double, float, int, and dcomplex.  For
    * all other types, a specialized implementation must be provided.
    *
    * @param boxes Const reference to box list describing the spatial extents
    *              of the index space regions of interest.
    *              Note: the boxes are assumed to be in the index space of this
    *              array data object.
    * @param src_shift Const reference to vector used to shift the given
    *              boxes into the index space region of this array data object.
    *              Note: this argument is currently ignored.
    */
   int
   getDataStreamSize(
      const hier::BoxContainer& boxes,
      const hier::IntVector& src_shift) const;

   /*!
    * Pack data living on the specified index region into the stream.
    *
    * Note that this routine assumes that the given box region must
    * be shifted to be consistent with the source (this) array box region.
    *
    * @param stream Reference to stream into which to pack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_shift Const reference to vector used to shift the given
    *              box into the index space region of this (source) array data
    *              object.
    *
    * Note: The shifted box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void
   packStream(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::IntVector& src_shift) const;

   void
   packStream(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::Transformation& src_shift) const;

   /*!
    * Pack data living on the specified index regions into the stream.
    *
    * Note that this routine assumes that the given box regions must
    * be shifted to be consistent with the source (this) array box region.
    *
    * @param stream Reference to stream into which to pack data.
    * @param dest_boxes Const reference to boxes describing the spatial extents
    *              of the destination index space regions of interest.
    * @param src_shift Const reference to vector used to shift the given
    *              boxes into the index space region of this (source) array data
    *              object.
    *
    * Note: The shifted boxes must lie completely within the index space of this
    * array.  If compiled with assertions enabled, the routine will abort if
    * the shifted boxes are not contained in the index space of this array.
    */
   void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::IntVector& src_shift) const;

   void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::Transformation& transformation) const;

   /*!
    * Unpack data from the stream into the index region specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_offset Const reference to vector used to offset
    *              box into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::IntVector& src_offset);

   /*!
    * Unpack data from the stream into the index regions specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_boxes Const reference to box list describing the spatial extents
    *              of the destination index space regions of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              boxes into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given boxes must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if some box is not contained in the index space of this array.
    */
   void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::IntVector& src_offset);

   /*!
    * Unpack data from the stream and add to the array in the index region specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_box Const reference to box describing the spatial extent
    *              of the destination index space region of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              box into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given box must lie completely within the index space of this
    * array data object.  When assertion checking is active, the routine will
    * abort if the box is not contained in the index space of this array.
    */
   void
   unpackStreamAndSum(
      tbox::MessageStream& stream,
      const hier::Box& dest_box,
      const hier::IntVector& src_offset);

   /*!
    * Unpack data from the stream and ad to the array in the index region specified.
    *
    * @param stream Reference to stream from which to unpack data.
    * @param dest_boxes Const reference to box list describing the spatial extents
    *              of the destination index space regions of interest.
    * @param src_offset Const reference to vector used to offset the given
    *              boxes into the index space region of some (source) array data
    *              object. Currently, this argument is ignored.
    *
    * Note: The given boxes must lie completely within the index space of this
    * array.  If compiled with assertions enabled, the routine will abort if
    * some box is not contained in the index space of this array.
    */
   void
   unpackStreamAndSum(
      tbox::MessageStream& stream,
      const hier::BoxContainer& dest_boxes,
      const hier::IntVector& src_offset);

   /*!
    * Fill all array values with value t.
    */
   void
   fillAll(
      const TYPE& t);

   /*!
    * Fill all array values within the box with value t.
    */
   void
   fillAll(
      const TYPE& t,
      const hier::Box& box);

   /*!
    * Fill all array values associated with depth component d with the value t.
    */
   void
   fill(
      const TYPE& t,
      const int d = 0);

   /*!
    * Fill all array values associated with depth component d
    * within the box with the value t.
    */
   void
   fill(
      const TYPE& t,
      const hier::Box& box,
      const int d = 0);

   /*!
    * Check to make sure that the class version and restart file
    * version are equal.  If so, read in data from database.  This
    * routine calls getSpecializedFromDatabase() to read in the
    * proper data type.
    *
    * Assertions:  database must be a non-null pointer.
    */
   void
   getFromDatabase(
      const boost::shared_ptr<tbox::Database>& database);

   /*!
    * Write out array data object data to database.  This
    * routine calls putSpecializedToDatabase() to read in the
    * proper data type.  The default behavior (boolean argument is
    * false) is to put all data members in database.  Otherwise, only
    * the array contents are written out.
    *
    * Assertions:  database must be a non-null pointer.
    */
   void
   putUnregisteredToDatabase(
      const boost::shared_ptr<tbox::Database>& database,
      bool data_only = false) const;

   /*!
    * Use specialized template method to get the correct behavior
    * when reading in the array of data.
    */
   void
   getSpecializedFromDatabase(
      const boost::shared_ptr<tbox::Database>& database);

   /*!
    * Use specialized template method to get the correct behavior
    * when writing out the array of data.
    */
   void
   putSpecializedToDatabase(
      const boost::shared_ptr<tbox::Database>& database) const;

   /**
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const;

   /**
    * @brief Invalidate an array as opposed to initializing it.
    *
    * The box associated with the array will be set to empty so
    * intersections will be the empty set.
    *
    */
   void
   invalidateArray(
      const tbox::Dimension& dim);

   /**
    * @brief Returns true if the array is valid.
    *
    * Invalid state can be set using the invalidateArray method.
    */
   bool
   isValid();

   /*!
    * The array data iterator iterates over the elements of a box
    * associated with an ArrayData object.  This typedef is
    * convenient link to the ArrayDataIterator class.
    */
   typedef ArrayDataIterator iterator;

private:
   ArrayData(
      const ArrayData<TYPE>&);          // not implemented
   void
   operator = (
      const ArrayData<TYPE>&);                  // not implemented

   /*
    * Static integer constant describing this class's version number.
    */
   static const int PDAT_ARRAYDATA_VERSION;

   /*
    * Private member functions to pack/unpack data to/from buffer.
    *
    * Note: box of this array data object must completely contain given box.
    */
   void
   packBuffer(
      TYPE* buffer,
      const hier::Box& box) const;
   void
   unpackBuffer(
      const TYPE* buffer,
      const hier::Box& box);

   /*
    * Private member functions to unpack data from buffer and add to
    * this array data object.
    *
    * Note: box of this array data object must completely contain given box.
    */
   void
   unpackBufferAndSum(
      const TYPE* buffer,
      const hier::Box& box);

   tbox::Dimension d_dim;

   int d_depth;
   int d_offset;
   hier::Box d_box;
   tbox::Array<TYPE> d_array;
};

}
}

#include "SAMRAI/pdat/ArrayData.C"

#endif
