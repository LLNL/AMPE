/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated face centered patch data type
 *
 ************************************************************************/

#ifndef included_pdat_FaceData
#define included_pdat_FaceData

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/PIO.h"

#include <boost/shared_ptr.hpp>
#include <iostream>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief Class FaceData<DIM> provides an implementation for data defined
 * at cell faces on AMR patches.  It is derived from the hier::PatchData
 * interface common to all SAMRAI patch data types.  Given a CELL-centered
 * AMR index space box, a face data object represents data of some template
 * TYPE and depth on the faces of the cells in the box.  Here, depth
 * indicates the number of data values at each face index location.  The
 * FaceGeometry class provides the translation between the standard SAMRAI
 * cell-centered AMR index space and face-centered data.
 *
 * Face data is stored in DIM arrays, each of which contains the
 * data for the faces normal to a corresponding coordinate direction.
 * In addition, the array indices are permuted so that the fastest array
 * dimension is the same as the normal direction of the face.
 * Memory allocation is in column-major ordering (e.g., Fortran
 * style) so that the leftmost index runs fastest in memory.
 * For example, a three-dimensional face data object created over a
 * CELL-centered AMR index space [l0:u0,l1:u1,l2:u2] allocates three data
 * arrays dimensioned as follows:
 * \verbatim
 *
 * face normal 0
 *   [ l0 : u0+1 ,
 *     l1 : u1 ,
 *     l2 : u2 , d ]   ,
 *
 * face normal 1
 *   [ l1 : u1+1 ,
 *     l2 : u2 ,
 *     l0 : u0 , d ]   ,
 *
 * face normal 2
 *   [ l2 : u2+1 ,
 *     l0 : u0 ,
 *     l1 : u1 , d ]   ,
 *
 * \endverbatim
 * Here the face normal directions 0, 1, 2 can be thought of as the x, y, and z
 * face normal directions, respectively, and d is the depth index (i.e., number
 * of values at each face index location).  Other spatial dimensions are
 * represented similarly.
 *
 * The data type TYPE must define a default constructor (that takes no
 * arguments) and also the assignment operator.
 *
 * IMPORTANT: The SideData<DIM> class provides the same storage
 * as this face data class, except that the coordinate directions of the
 * individual arrays are not permuted.
 *
 * @see pdat::ArrayData
 * @see hier::PatchData
 * @see pdat::FaceDataFactory
 * @see pdat::FaceIndex
 * @see pdat::FaceIterator
 * @see pdat::FaceGeometry
 */

template<class TYPE>
class FaceData:public hier::PatchData
{
public:
   /*!
    * @brief Calculate the amount of memory needed to represent face-
    * centered data over a CELL-centered AMR index space box.
    *
    * This function assumes that the amount of memory
    * needed for TYPE is sizeof(TYPE).  If this is not the case, then a
    * specialized function must be defined.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            face data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *               of the ghost cell region around the box over which
    *               the face data will be allocated.
    */
   static size_t
   getSizeOfData(
      const hier::Box& box,
      int depth,
      const hier::IntVector& ghosts);

   /*!
    * @brief The constructor for a face data object.
    *
    * @param box const Box reference describing the interior of the
    *            standard CELL-centered index box over which the
    *            face data object will be created.
    * @param depth gives the number of components for each
    *              spatial location in the array.
    * @param ghosts const IntVector reference indicating the width
    *               of the ghost cell region around the box over which
    *               the face data will be allocated.
    */
   FaceData(
      const hier::Box& box,
      int depth,
      const hier::IntVector& ghosts);

   /*!
    * @brief The virtual destructor for a face data object.
    */
   virtual ~FaceData<TYPE>();

   /*!
    * @brief Return the depth (e.g., the number of components in each spatial
    * location) of the array.
    */
   int
   getDepth() const;

   /*!
    * @brief Get a pointer to the beginning of a particular face normal and
    * depth component of the face centered array.
    */
   TYPE *
   getPointer(
      int face_normal,
      int depth = 0);

   /*!
    * @brief Get a const pointer to the beginning of a particular face normal
    * and depth component of the face centered array.
    */
   const TYPE *
   getPointer(
      int face_normal,
      int depth = 0) const;

   /*!
    * @brief Return a reference to the data entry corresponding
    * to a given face index and depth.
    */
   TYPE&
   operator () (
      const FaceIndex& i,
      int depth = 0);

   /*!
    * @brief Return a const reference to the data entry corresponding
    * to a given face index and depth.
    */
   const TYPE&
   operator () (
      const FaceIndex& i,
      int depth = 0) const;

   /*!
    * @brief Return a reference to the array data object for the
    * given face normal of the face centered data object.
    */
   ArrayData<TYPE>&
   getArrayData(
      int face_normal);

   /*!
    * @brief Return a const reference to the array data object for the
    * given face normal of the face centered data object.
    */
   const ArrayData<TYPE>&
   getArrayData(
      int face_normal) const;

   /*!
    * @brief A fast copy from source to destination (i.e., this)
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, source data must be
    * an FaceData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void
   copy(
      const hier::PatchData& src);

   /*!
    * @brief A fast copy from source (i.e., this) to destination
    * patch data object.
    *
    * Data is copied where there is overlap in the underlying index space.
    * The copy is performed on the interior plus the ghost cell width (for
    * both the source and destination).  Currently, destination data must be
    * an FaceData of the same DIM and TYPE.  If not, then an unrecoverable
    * error results.
    */
   virtual void
   copy2(
      hier::PatchData& dst) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given overlap.
    *
    * Currently, source data must be FaceData of the same DIM and TYPE
    * and the overlap must be a FaceOverlap of the same DIM. If not,
    * then an unrecoverable error results.
    */
   virtual void
   copy(
      const hier::PatchData& src,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Copy data from source (i.e., this) to destination
    * patch data object on the given overlap.
    *
    * Currently, destination data must be FaceData of the same DIM and TYPE
    * and the overlap must be a FaceOverlap of the same DIM.  If not,
    * then an unrecoverable error results.
    */
   virtual void
   copy2(
      hier::PatchData& dst,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Copy data from source to destination (i.e., this)
    * patch data object on the given CELL-centered AMR index box.
    */
   void
   copyOnBox(
      const FaceData<TYPE>& src,
      const hier::Box& box);

   /*!
    * @brief Fast copy (i.e., source and this face data objects are
    * defined over the same box) to this destination face data object
    * from the given source face data object at the specified depths.
    */
   void
   copyDepth(
      int dst_depth,
      const FaceData<TYPE>& src,
      int src_depth);

   /*!
    * @brief Return true if the patch data object can estimate the
    * stream size required to fit its data using only index
    * space information (i.e., a box).
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual bool
   canEstimateStreamSizeFromBox() const;

   /*!
    * @brief Return the number of bytes needed to stream the data
    * in this patch data object lying in the specified box overlap
    * region.
    *
    * This routine is defined for the standard types (bool, char,
    * double, float, int, and dcomplex).
    */
   virtual int
   getDataStreamSize(
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Pack data in this patch data object lying in the specified
    * box overlap region into the stream.  The overlap must be an
    * FaceOverlap of the same DIM.
    */
   virtual void
   packStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap) const;

   /*!
    * @brief Unpack data from stream into this patch data object over
    * the specified box overlap region. The overlap must be an
    * FaceOverlap of the same DIM.
    */
   virtual void
   unpackStream(
      tbox::MessageStream& stream,
      const hier::BoxOverlap& overlap);

   /*!
    * @brief Fill all values at depth d with the value t.
    */
   void
   fill(
      const TYPE& t,
      int d = 0);

   /*!
    * @brief Fill all values at depth d within the box with the value t.
    */
   void
   fill(
      const TYPE& t,
      const hier::Box& box,
      int d = 0);

   /*!
    * @brief Fill all depth components with value t.
    */
   void
   fillAll(
      const TYPE& t);

   /*!
    * @brief Fill all depth components within the box with value t.
    */
   void
   fillAll(
      const TYPE& t,
      const hier::Box& box);

   /*!
    * @brief Print all face data values residing in the specified box.
    * If the depth of the array is greater than one, all depths are printed.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void
   print(
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all face data values at the given array depth in
    * the specified box.
    *
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os   reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void
   print(
      const hier::Box& box,
      int depth,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all face centered data values for specified face normal
    * direction residing in the specified box.  If the depth of the data is
    * greater than one, all depths are printed.
    *
    * @param face_normal  integer face normal coordinate direction,
    *              must satisfy 0 <= face_normal < DIM
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void
   printAxis(
      int face_normal,
      const hier::Box& box,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * @brief Print all face centered data values for specified face normal
    * direction residing in the specified box.  If the depth of the data is
    * greater than one, all depths are printed.
    *
    * @param face_normal  integer face normal coordinate direction,
    *              must satisfy 0 <= face_normal < DIM
    * @param box  const reference to box over whioch to print data. Note box
    *        is assumed to reside in standard cell-centered index space
    *        and will be converted to face index space.
    * @param depth integer depth component, must satisfy
    *              0 <= depth < actual depth of data array
    * @param os    reference to output stream.
    * @param prec integer precision for printing floating point numbers
    *        (i.e., TYPE = float, double, or dcomplex). The default
    *        is 12 decimal places for double and complex floating point numbers,
    *        and the default is 6 decimal places floats.  For other types, this
    *        value is ignored.
    */
   void
   printAxis(
      int face_normal,
      const hier::Box& box,
      int depth,
      std::ostream& os = tbox::plog,
      int prec = 12) const;

   /*!
    * Check that class version and restart file version are equal.  If so,
    * read data members from the database.
    *
    * Assertions: database must be non-null pointer.
    */
   virtual void
   getSpecializedFromDatabase(
      const boost::shared_ptr<tbox::Database>& database);

   /*!
    * Write out the class version number and other data members to
    * the database.
    *
    * Assertions: database must be non-null pointer.
    */
   virtual void
   putSpecializedToDatabase(
      const boost::shared_ptr<tbox::Database>& database) const;

   /*!
    * The face iterator iterates over the elements on one face normal of a face
    * centered box geometry.  This typedef is a convenience for using the
    * FaceIterator class.
    */
   typedef FaceIterator iterator;

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int PDAT_FACEDATA_VERSION;

   FaceData(
      const FaceData<TYPE>&);           // not implemented
   void
   operator = (
      const FaceData<TYPE>&);                           // not implemented

   void
   copyWithRotation(
      const FaceData<TYPE>& src,
      const FaceOverlap& overlap);

   void
   packWithRotation(
      tbox::MessageStream& stream,
      const FaceOverlap& overlap) const;

   int d_depth;
   ArrayData<TYPE> d_data[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

};

}
}

#include "SAMRAI/pdat/FaceData.C"

#endif
