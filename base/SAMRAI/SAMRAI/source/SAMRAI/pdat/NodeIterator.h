/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Iterator for node centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_NodeIterator
#define included_pdat_NodeIterator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/hier/Box.h"

namespace SAMRAI {
namespace pdat {

/**
 * Class NodeIterator is an iterator that provides methods for
 * stepping through the index space associated with a node centered box.
 * The indices are enumerated in column-major (e.g., Fortran) order.
 * The iterator should be used as follows:
 * \verbatim
 * hier::Box box;
 * ...
 * NodeIterator cend(box, false);
 * for (NodeIterator c(box, true); c != cend; ++c) {
 *    // use index c of the box
 * }
 * \endverbatim
 * Note that the node iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see pdat::NodeData
 * @see pdat::NodeGeometry
 * @see pdat::NodeIndex
 */

class NodeIterator
{
public:
   /**
    * Constructor for the node iterator.  The iterator will enumerate
    * the indices in the argument box.
    */
   explicit NodeIterator(
      const hier::Box& box,
      bool begin);

   /**
    * Copy constructor for the node iterator
    */
   NodeIterator(
      const NodeIterator& iterator);

   /**
    * Assignment operator for the node iterator.
    */
   NodeIterator&
   operator = (
      const NodeIterator& iterator)
   {
      d_index = iterator.d_index;
      d_box = iterator.d_box;
      return *this;
   }

   /**
    * Destructor for the node iterator.
    */
   ~NodeIterator();

   /**
    * Extract the node index corresponding to the iterator position in the box.
    */
   const NodeIndex&
   operator * () const
   {
      return d_index;
   }

   /**
    * Extract a pointer to the node index corresponding to the iterator
    * position in the box.
    */
   const NodeIndex*
   operator -> () const
   {
      return &d_index;
   }

   /**
    * Pre-increment the iterator to point to the next index in the box.
    */
   NodeIterator&
   operator ++ ();

   /**
    * Post-increment the iterator to point to the next index in the box.
    */
   NodeIterator
   operator ++ (
      int);

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const NodeIterator& iterator) const
   {
      return d_index == iterator.d_index;
   }

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const NodeIterator& iterator) const
   {
      return d_index != iterator.d_index;
   }

private:
   NodeIndex d_index;
   hier::Box d_box;
};

}
}

#endif
