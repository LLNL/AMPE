/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Iterator for edge centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_EdgeIterator
#define included_pdat_EdgeIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/EdgeIndex.h"
#include "SAMRAI/hier/Box.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class EdgeIterator is an iterator that provides methods for
 * stepping through the index space associated with a edge centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
 * \verbatim
 * hier::Box box;
 * ...
 * EdgeIterator cend(box, axis, false);
 * for (EdgeIterator c(box, axis, true); c != cend; ++c) {
 *    // use index c of the box
 * }
 * \endverbatim
 * Note that the edge iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::EdgeData
 * @see pdat::EdgeGeometry
 * @see pdat::EdgeIndex
 */

class EdgeIterator
{
public:
   /**
    * Constructor for the edge iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   EdgeIterator(
      const hier::Box& box,
      const int axis,
      bool begin);

   /**
    * Copy constructor for the edge iterator
    */
   EdgeIterator(
      const EdgeIterator& iterator);

   /**
    * Assignment operator for the edge iterator.
    */
   EdgeIterator&
   operator = (
      const EdgeIterator& iterator)
   {
      d_index = iterator.d_index;
      d_box = iterator.d_box;
      return *this;
   }

   /**
    * Destructor for the edge iterator.
    */
   ~EdgeIterator();

   /**
    * Extract the edge index corresponding to the iterator position in the box.
    */
   const EdgeIndex&
   operator * () const
   {
      return d_index;
   }

   /**
    * Extract a pointer to the edge index corresponding to the iterator
    * position in the box.
    */
   const EdgeIndex*
   operator -> () const
   {
      return &d_index;
   }

   /**
    * Pre-increment the iterator to point to the next index in the box.
    */
   EdgeIterator&
   operator ++ ();

   /**
    * Post-increment the iterator to point to the next index in the box.
    */
   EdgeIterator
   operator ++ (
      int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const EdgeIterator& iterator) const
   {
      return d_index == iterator.d_index;
   }

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const EdgeIterator& iterator) const
   {
      return d_index != iterator.d_index;
   }

private:
   EdgeIndex d_index;
   hier::Box d_box;
};

}
}

#endif
