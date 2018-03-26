/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Box representing a portion of the AMR index space
 *
 ************************************************************************/

#ifndef included_hier_Box
#define included_hier_Box

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/BoxId.h"
#include "SAMRAI/hier/Transformation.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iostream>

/*
 * Forward declaration, which is questionable with respect to SAMRAI 
 * package ordering.   This is needed since pdat::ArrayData class
 * needs to access private Box default constructor.
 *
 * It would be good to come up with an alternative to this.
 */
namespace SAMRAI {

namespace pdat {
template<class TYPE>
class ArrayData;
}

}

namespace SAMRAI {
namespace hier {

class BoxIterator;

/**
 * Class Box represents a n-dimensional box in the AMR index
 * space.  It is defined by lower and upper bounds given by index objects.
 * The box semantics assumes that the box is cell-centered.  A cell-centered
 * convention implies that the index set covered by the box includes both
 * the lower and upper bounds.
 *
 * The Box contains identifying information in its state with a BlockId and
 * a BoxId.  If the Box is part of a single-block mesh, the BlockId should be
 * zero.  If the mesh is multi-block, the BlockId will have a value identifying
 * on which block the Box exists.  The BoxId contains information about the
 * MPI rank associated with the Box as well as a PeriodicId to handle periodic
 * shifts.  Since periodic conditions and multiblock meshes cannot be mixed in
 * SAMRAI, the BlockId and the PeriodicId associated with a Box cannot both be
 * nonzero.
 *
 * @see hier::BoxIterator
 * @see hier::Index
 * @see hier::BlockId
 * @see hier::BoxId
 * @see hier::PeriodicId
 */

class Box
{
public:
   /*!
    * @brief Creates an empty box with invalid BlockId and BoxId values.
    *
    * @param[in]  dim
    */
   explicit Box(
      const tbox::Dimension& dim);

   /*!
    * Create a box describing the index space between lower and upper.  The
    * box is assumed to be cell centered and include all elements between lower
    * and upper, including the end points.
    *
    * @param[in] lower     Lower extent
    * @param[in] upper     Upper extent
    * @param[in] block_id  Block where the Box exists
    */
   Box(
      const Index& lower,
      const Index& upper,
      const BlockId& block_id);

   /*!
    * @brief Copy constructor 
    */
   Box(
      const Box& box);

   /*!
    * Construct a Box from a DatabaseBox.
    */
   explicit Box(
      const tbox::DatabaseBox& box);

   /*!
    * @brief Initializing constructor.
    *
    * @param[in] box
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] periodic_id Describes the periodic shift.  If
    * periodic_id is non-zero, specify the Box in the position shifted
    * according to the @c periodic_id.  The default argument for @c
    * periodic_id corresponds to the zero-shift.
    */
   Box(
      const Box& box,
      const LocalId& local_id,
      const int owner_rank,
      const PeriodicId& periodic_id = PeriodicId::zero());

   /*!
    * @brief Constructor with undefined box.
    *
    * The box can be initialized using any of the initialize()
    * methods or by assignment.
    *
    * @param[in] dim
    *
    * @param[in] id
    *
    * @param[in] periodic_id
    */
   /*
    * TODO: Constructors initializing boxes are only used to construct
    * temporary objects for finding other Boxes in a
    * stl::set<Box>.  We need another way to do it and get rid
    * of these constructors.
    */
   Box(
      const tbox::Dimension& dim,
      const GlobalId& id,
      const PeriodicId& periodic_id = PeriodicId::zero());

   /*!
    * @brief Constructor with undefined box and a BoxId.
    *
    * The box can be initialized using any of the initialize()
    * methods or by assignment.
    *
    * @param[in] dim
    *
    * @param[in] mapped_box_id
    */
   /*
    * TODO: Constructors initializing boxes are only used to construct
    * temporary objects for finding other Boxes in a
    * stl::set<Box>.  We need another way to do it and get rid
    * of these constructors.
    */
   Box(
      const tbox::Dimension& dim,
      const BoxId& mapped_box_id);

   /*!
    * @brief "Copy" constructor allowing change in PeriodicId.
    *
    * @param[in] other Make a copy (but not an exact copy) of this
    * Box.
    *
    * @param[in] periodic_id Periodic shift number to use instead of
    * the shift in @c other.  The box will be set to the real box
    * shifted to the position specified by this value.
    *
    * @param[in] refinement_ratio The index space where the Box
    * lives.
    *
    * @see initialize( const Box&, const int, const IntVector&);
    */
   Box(
      const Box& other,
      const PeriodicId& periodic_id,
      const IntVector& refinement_ratio);

   /*!
    * @brief The destructor for Box.
    */
   ~Box();

   /*!
    * @brief Set all the attributes of the Box.
    *
    * @param[in] box
    *
    * @param[in] local_id
    *
    * @param[in] owner_rank
    *
    * @param[in] periodic_id The periodic shift number.  If
    * this is not zero, specify @c box in the shifted position.  The
    * default argument for @c periodic_id corresponds to the
    * zero-shift.
    */
   void
   initialize(
      const Box& box,
      const LocalId& local_id,
      const int owner_rank,
      const PeriodicId& periodic_id = PeriodicId::zero())
   {
      d_lo = box.d_lo;
      d_hi = box.d_hi;
      d_block_id = box.d_block_id;
      if (!d_id_locked) {
         d_id.initialize(local_id, owner_rank, periodic_id);
      } else {
         TBOX_ERROR("Attempted to change BoxId that is locked in an ordered BoxContainer.");
      }
   }

   /*!
    * @brief Set all the attributes identical to that of a reference
    * Box, but with a different PeriodicId.
    *
    * @param[in] other Initialize to this Box, but with the shift
    * given by @c periodic_id.
    *
    * @param[in] periodic_id PeriodicId number to use instead of the
    * shift in @c other.  The box will be set to the real box shifted
    * to the position specified by this value.
    *
    * @param[in] refinement_ratio The index space where the Box
    * lives.
    */
   void
   initialize(
      const Box& other,
      const PeriodicId& periodic_id,
      const IntVector& refinement_ratio);

   //! @brief Set the BoxId.
   void
   setId(
      const BoxId& box_id)
   {
      d_id = box_id;
   }

   //! @brief Get the BoxId.
   const BoxId&
   getId() const
   {
      return d_id;
   }

   //! @brief Set the BlockId.
   void
   setBlockId(
      const BlockId& block_id)
   {
      d_block_id = block_id;
   }

   //! @brief Get the BlockId.
   const BlockId&
   getBlockId() const
   {
      return d_block_id;
   }

   //! @brief Get the LocalId.
   const LocalId&
   getLocalId() const
   {
      return d_id.getLocalId();
   }

   //! @brief Get the GlobalId.
   const GlobalId&
   getGlobalId() const
   {
      return d_id.getGlobalId();
   }

   //! @brief Get the owner rank.
   int
   getOwnerRank() const
   {
      return d_id.getOwnerRank();
   }

   /*!
    * @brief Get the periodic shift number.
    *
    * @see PeriodicShiftCatalog.
    */
   const PeriodicId&
   getPeriodicId() const
   {
      return d_id.getPeriodicId();
   }

   //! @brief Whether the Box is a periodic image.
   bool
   isPeriodicImage() const
   {
      return d_id.isPeriodicImage();
   }

   bool
   isIdEqual(
      const Box& other) const
   {
      return d_id == other.d_id;
   }

   struct id_equal {
      bool
      operator () (const Box& b1, const Box& b2) const
      {
         return b1.isIdEqual(b2);
      }

      bool
      operator () (const Box* b1, const Box* b2) const
      {
         return b1->isIdEqual(*b2);
      }
   };

   struct id_less {
      bool
      operator () (const Box& b1, const Box& b2) const
      {
         return b1.getId() < b2.getId();
      }

      bool
      operator () (const Box* b1, const Box* b2) const
      {
         return b1->getId() < b2->getId();
      }
   };

   //@{
   //! @name Packing and unpacking Boxes

   /*!
    * @brief Give number of ints required for putting a Box in
    * message passing buffer.
    *
    * This number is independent of instance but dependent on
    * dimension.
    *
    * @see putToIntBuffer(), getFromIntBuffer().
    */
   static int
   commBufferSize(
      const tbox::Dimension& dim)
   {
      return BoxId::commBufferSize() + (2 * dim.getValue()) + 1;
   }

   /*!
    * @brief Put self into a int buffer.
    *
    * This is the opposite of getFromIntBuffer().  Number of ints
    * written is given by commBufferSize().
    */
   void
   putToIntBuffer(
      int * buffer) const;

   /*!
    * @brief Set attributes according to data in int buffer.
    *
    * This is the opposite of putToIntBuffer().  Number of ints read
    * is given by commBufferSize().
    */
   void
   getFromIntBuffer(
      const int * buffer);

   /*!
    * @brief Put self into a MessageStream.
    *
    * This is the opposite of getFromMessageStream().
    *
    * @param msg [in] Stream.  Must be in write mode.
    */
   void
   putToMessageStream(
      tbox::MessageStream &msg) const;

   /*!
    * @brief Set attributes according to data in MessageStream.
    *
    * This is the opposite of putToMessageStream().
    *
    * @param msg [in] Stream.  Must be in read mode.
    */
   void
   getFromMessageStream(
      tbox::MessageStream &msg);

   //@}

   /*!
    * @brief assignment operator
    *
    * An assignment to an uninitialized box is allowed but assigning
    * from an uninitialized box will result in an assertion failure
    */
   Box&
   operator = (
      const Box& box);

   /*!
    * @brief Return a non-const lower index of the box.
    */
   Index&
   lower()
   {
      return d_lo;
   }

   /*!
    * @brief Return a non-const upper index of the box.
    */
   Index&
   upper()
   {
      return d_hi;
   }

   /*!
    * @brief Return a const lower index of the box.
    */
   const Index&
   lower() const
   {
      return d_lo;
   }

   /*!
    * @brief Return a const upper index of the box.
    */
   const Index&
   upper() const
   {
      return d_hi;
   }

   /*!
    * @brief Return the i'th component (non-const) of the lower index.
    */
   int&
   lower(
      const int i)
   {
      return d_lo(i);
   }

   /*!
    * @brief Return the i'th component (non-const) of the upper index.
    */
   int&
   upper(
      const int i)
   {
      return d_hi(i);
   }

   /*!
    * @brief Return the i'th component (const) of the lower index.
    */
   const int&
   lower(
      const int i) const
   {
      return d_lo(i);
   }

   /*!
    * @brief Return the i'th component (const) of the upper index.
    */
   const int&
   upper(
      const int i) const
   {
      return d_hi(i);
   }

   /*!
    * @brief Set the state of the box to empty.
    */
   void
   setEmpty()
   {
      const tbox::Dimension& dim(getDim());
      d_lo = Index(dim, tbox::MathUtilities<int>::getMax());
      d_hi = Index(dim, tbox::MathUtilities<int>::getMin());
   }

   /*!
    * @brief Return whether the box is ``empty''.
    *
    * isEmpty() is preferred to match "is" standard syntax for
    * boolean methods.
    *
    * @see isEmpty()
    */
   bool
   empty() const
   {
      for (int i = 0; i < getDim().getValue(); i++) {
         if (d_hi(i) < d_lo(i)) {
            return true;
         }
      }
      return false;
   }

   /*!
    * @brief Return whether the box is ``empty''.
    *
    * A box is empty if any of the lower bounds is greater than the
    * corresponding upper bound.  An empty box has a size of zero.
    */
   bool
   isEmpty() const
   {
      for (int i = 0; i < getDim().getValue(); i++) {
         if (d_hi(i) < d_lo(i)) {
            return true;
         }
      }
      return false;
   }

   /*!
    * @brief Return the number of cells (an integer) represented by the box in
    * the given coordinate direction.
    */
   int
   numberCells(
      const int i) const
   {
      if (empty()) {
         return 0;
      } else {
         return d_hi(i) - d_lo(i) + 1;
      }
   }

   /*!
    * @brief Return the number of cells (a vector of integers) represented by
    * the box in every coordinate direction.
    */
   IntVector
   numberCells() const
   {
      if (empty()) {
         return IntVector::getZero(getDim());
      } else {
         return d_hi - d_lo + 1;
      }
   }

   /*!
    * @brief Calculate the number of indices represented by the box.
    * 
    * If the box is empty, then the number of index points within the box is
    * zero.
    */
   int
   size() const
   {
      int mysize = 0;
      if (!empty()) {
         mysize = 1;
         for (int i = 0; i < getDim().getValue(); i++) {
            mysize *= (d_hi(i) - d_lo(i) + 1);
         }
      }
      return mysize;
   }

   /*!
    * @brief Return the dimension of the box that is longest.
    */
   int
   longestDimension() const;

   /*!
    * @brief Given an index, calculate the offset of the index into the box.
    *
    * This function assumes column-major (e.g., Fortran) ordering of
    * the indices within the box.  This operation is a convenience
    * function for array indexing operations.
    */
   int
   offset(
      const Index& p) const
   {
      int myoffset = 0;
      for (int i = getDim().getValue() - 1; i > 0; i--) {
         myoffset = (d_hi(i - 1) - d_lo(i - 1) + 1) * (p(i) - d_lo(i) + myoffset);
      }
      myoffset += p(0) - d_lo(0);
      return myoffset;
   }

   /*!
    * @brief Given an offset, calculate the index of the offset into the box.
    *
    * This function assumes column-major (e.g., Fortran) ordering of
    * the indices within the box.  This operation is a convenience
    * function for array indexing operations.
    */
   Index
   index(
      const int offset) const;

   /*!
    * @brief Check whether an index lies within the bounds of the box.
    */
   bool
   contains(
      const Index& p) const
   {
      for (int i = 0; i < getDim().getValue(); i++) {
         if ((p(i) < d_lo(i)) || (p(i) > d_hi(i))) {
            return false;
         }
      }
      return true;
   }

   /*!
    * @brief Check whether a given box lies within the bounds of the box.
    *
    * If @c b is empty, always return true.
    */
   bool
   contains(
      const Box& b) const;

   /*!
    * @brief Check whether two boxes represent the same portion of index space
    * on the same block.
    */
   bool
   isSpatiallyEqual(
      const Box& box) const
   {
      return ((d_lo == box.d_lo) && (d_hi == box.d_hi) &&
              (d_block_id == box.d_block_id)) || (empty() && box.empty());
   }
   struct box_equality {
      bool
      operator () (const Box& b1, const Box& b2) const
      {
         return b1.isSpatiallyEqual(b2);
      }
   };

   /*!
    * @brief Calculate the intersection of the index spaces of two boxes.
    *
    * The intersection with an empty box always yields an empty box.  If the
    * two boxes have different BlockIds and are both non-empty, an error will
    * occur.
    */
   Box
   operator * (
      const Box& box) const;

   /*!
    * @brief Calculate the intersection of the index spaces of this and the
    * given box.
    *
    * The intersection with an empty box always yields an empty box.  If the
    * two boxes have different BlockIds and are both non-empty, an error will
    * occur.
    */
   Box
   &operator *= (
      const Box& box);

   /*!
    * @brief Box intersection.
    * 
    * Calculate the intersection of the index spaces of two boxes.  The
    * intersection with an empty box always yields an empty box.  If the
    * two boxes have different BlockIds and are both non-empty, an error will
    * occur.
    */
   void
   intersect(
      const Box& other,
      Box& result) const;

   /*!
    * @brief check if boxes intersect.
    *
    * Returns true if two boxes have a non-empty intersection.
    */
   bool
   intersects(
      const Box& box) const;

   /*!
    * @brief Calculate the bounding box of two boxes.
    *
    * Note that this is not the union of the two boxes (since union
    * is not closed over boxes), but rather the smallest box that
    * contains both boxes.
    *
    * If one box is empty and the other is non-empty, the non-empty box will
    * be returned.  If both boxes are non-empty and have different BlockIds,
    * an assertion failure will occur.
    */
   Box
   operator + (
      const Box& box) const;

   /*!
    * @brief Increase this box to become the bounding box of this box and the
    * argument box.
    *
    * If this box is empty, this box will become equal to the argument box.
    *
    * If the argument box is empty, this box will be unchanged.
    *
    * If the two boxes are both non-empty and have different BlockIds, and
    * assertion failure will occur.
    */
   Box&
   operator += (
      const Box& box);

   /*!
    * @brief Return true if this box can be coalesced with the argument box,
    * and set this box to the union of the boxes.
    * 
    * Otherwise, return false and leave boxes as is.  Two boxes may be
    * coalesced if their union is a box (recall that index set union is not
    * closed over boxes).  If one box is empty and the other is non-empty, then
    * the return value is true and this box becomes equal to the non-empty box.
    * If both boxes are empty, the return value is true and this box remains
    * empty.
    */
   bool
   coalesceWith(
      const Box& box);

   /*!
    * @brief Grow this box by the specified ghost cell width.
    *
    * The lower bound decremented by the width, and the upper bound is
    * incremented by the width.  All dimensions are grown by the corresponding
    * component in the IntVector; ghost cell widths may be different in each
    * dimension.  Negative ghost cell widths will shrink the box.
    */
   void
   grow(
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);
      if (!empty()) {
         d_lo -= ghosts;
         d_hi += ghosts;
      }
   }

   /*!
    * @brief Grow this box by the specified ghost cell width in the given
    * coordinate direction.
    *
    * The lower bound is decremented by the width, and the upper bound is
    * incremented by the width.  Note that negative ghost cell widths
    * will shrink the box.
    */
   void
   grow(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));
      if (!empty()) {
         d_lo(direction) -= ghosts;
         d_hi(direction) += ghosts;
      }
   }

   /*!
    * @brief Similar to grow() functions. However, box is only grown in lower
    * directions (i.e., only lower index is changed).
    */
   void
   growLower(
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);
      if (!empty()) {
         d_lo -= ghosts;
      }
   }

   /*!
    * @brief Similar to grow() functions. However, box is only grown in lower
    * bound of given direction in index space.
    */
   void
   growLower(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));
      if (!empty()) {
         d_lo(direction) -= ghosts;
      }
   }

   /*!
    * @brief Similar to grow() function. However, box is only grown in upper
    * directions (i.e., only upper index is changed).
    */
   void
   growUpper(
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ghosts);
      if (!empty()) {
         d_hi += ghosts;
      }
   }

   /*!
    * @brief Similar to grow() functions. However, box is only grown in upper
    * bound of given direction in index space.
    */
   void
   growUpper(
      const int direction,
      const int ghosts)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));
      if (!empty()) {
         d_hi(direction) += ghosts;
      }
   }

   /*!
    * @brief Similar to growUpper() and growLower() functions. However, box is
    * lengthened (never shortened).
    *
    * The sign of @c ghosts refer to whether the box is lengthened in
    * the upper or lower side.
    */
   void
   lengthen(
      const int direction,
      const int ghosts);

   /*!
    * @brief Similar to growUpper() and growLower() functions. However, box is
    * shortened (never lengthened).
    *
    * The sign of @c ghosts refer to whether the box is shortened in
    * the upper or lower side.
    */
   void
   shorten(
      const int direction,
      const int ghosts);

   /*!
    * @brief Shift this box by the specified offset.
    *
    * The box will be located at (lower+offset, upper+offset).
    */
   void
   shift(
      const IntVector& offset)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, offset);
      d_lo += offset;
      d_hi += offset;
   }

   /*!
    * @brief Similar to shift() function above, but shift occurs only in
    * specified direction in index space.
    *
    * The box will located at (lower+offset, upper+offset) in that direction.
    */
   void
   shift(
      const int direction,
      const int offset)
   {
      TBOX_ASSERT((direction >= 0) && (direction < getDim().getValue()));
      d_lo(direction) += offset;
      d_hi(direction) += offset;
   }

   /*!
    * Rotate 90 degrees around origin.
    */
   void
   rotate(
      const Transformation::RotationIdentifier rotation_ident);

   /*!
    * @brief Refine the index space of a box by specified vector ratio.
    *
    * Each component of the box is multiplied by the refinement ratio,
    * then @c (ratio-1) is added to the upper corner.
    */
   void
   refine(
      const IntVector& ratio);

   /*!
    * @brief Coarsen the index space of a box by specified vector ratio.
    *
    * Each component is divided by the specified coarsening ratio and rounded
    * (if necessary) such that the coarsened box contains the cells that
    * are the parents of the refined box.  In other words, refining a
    * coarsened box will always yield a box that is equal to or larger
    * than the original box.
    */
   void
   coarsen(
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);
      for (int i = 0; i < getDim().getValue(); i++) {
         d_lo(i) = coarsen(d_lo(i), ratio(i));
         d_hi(i) = coarsen(d_hi(i), ratio(i));
      }
   }

   /*!
    * @brief This assignment operator constructs a Box given a DatabaseBox.
    */
   Box&
   operator = (
      const tbox::DatabaseBox& box)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, box);
#ifdef BOX_TELEMETRY
      // Increment the cumulative assigned count only.
      ++s_cumulative_assigned_ct;
#endif
      return Box_from_DatabaseBox(box);
   }

   /*!
    * @brief Sets a Box from a tbox::DatabaseBox and returns a reference to
    * the Box.
    */
   Box&
   Box_from_DatabaseBox(
      const tbox::DatabaseBox& box)
   {
      set_Box_from_DatabaseBox(box);
      return *this;
   }

   /*!
    * @brief Sets a Box from a DatabaseBox.
    */
   void
   set_Box_from_DatabaseBox(
      const tbox::DatabaseBox& box);

   /*!
    * @brief Returns a tbox::DatabaseBox generated from a Box.
    */
   tbox::DatabaseBox
   DatabaseBox_from_Box() const;

   /*!
    * Type conversion from Box to DatabaseBox
    */
   operator tbox::DatabaseBox()
   {
      return DatabaseBox_from_Box();
   }

   /*!
    * Type conversion from Box to DatabaseBox
    */
   operator tbox::DatabaseBox() const
   {
      return DatabaseBox_from_Box();
   }

   /*!
    * @brief Static function to grow a box by the specified vector ghost cell
    * width.
    */
   static Box
   grow(
      const Box& box,
      const IntVector& ghosts)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, ghosts);
      Box tmp = box;
      tmp.grow(ghosts);
      return tmp;
   }

   /*!
    * @brief Static function to shift a box by the specified offset.
    */
   static Box
   shift(
      const Box& box,
      const IntVector& offset)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, offset);
      return Box(box.lower() + offset, box.upper() + offset, box.d_block_id);
   }

   /*!
    * @brief Static function to refine the index space of a box by the
    * specified refinement ratio.
    */
   static Box
   refine(
      const Box& box,
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);
      Box tmp = box;
      tmp.refine(ratio);
      return tmp;
   }

   /*!
    * @brief Static function to coarsen the index space of a box by the
    * specified coarsening ratio.
    */
   static Box
   coarsen(
      const Box& box,
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(box, ratio);
      Box tmp = box;
      tmp.coarsen(ratio);
      return tmp;
   }

   /*!
    * @brief Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_lo.getDim();
   }

   /*!
    * @brief Lock the BoxId of this Box so that it may not be changed in any
    * way.
    */
   void
   lockId()
   {
      d_id_locked = true;
   }

   /*!
    * @brief Read the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      Box& box);

   /*!
    * @brief Output the box description in the form [L,U], where L and U are the
    * lower and upper bounds of the box.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Box& box);

   /*!
    * @brief Return an empty Box of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Box&
   getEmptyBox(
      const tbox::Dimension& dim)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(dim);
      return *(s_emptys[dim.getValue() - 1]);
   }

   /*!
    * @brief Returns a Box that represents the maximum allowed index extents
    * for a given dimension.   The "universe" that can be represented.
    */
   static const Box&
   getUniverse(
      const tbox::Dimension& dim)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(dim);
      return *(s_universes[dim.getValue() - 1]);
   }

   /**
    * A box iterator iterates over the elements of a box.  This class is
    * defined elsewhere, and the typedef is used to point to that class.
    */
   typedef BoxIterator iterator;

   template<class>
   friend class pdat::ArrayData;
   friend class BoxIterator;
   friend class std::vector<Box>;

#ifdef BOX_TELEMETRY
   // These are to optionally track the cumulative number of Boxes constructed,
   // the cumulative number of Box assignments, and the high water mark of
   // Boxes in existance at any given time.
   static int s_cumulative_constructed_ct;

   static int s_cumulative_assigned_ct;

   static int s_active_ct;

   static int s_high_water;
#endif

private:
   /**
    * The default constructor creates an uninitialized box.
    *
    * This should never be invoked, it will cause assertions
    */
   Box();

   static int
   coarsen(
      const int index,
      const int ratio)
   {
      return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
   }

   static bool
   coalesceIntervals(
      const int* lo1,
      const int* hi1,
      const int* lo2,
      const int* hi2,
      const int dim);

   void
   rotateAboutAxis(
      const int axis,
      const int num_rotations);


   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   Index d_lo;
   Index d_hi;
   BlockId d_block_id;
   BoxId d_id;
   bool d_id_locked;

   /*
    * Array of empty boxes for each dimension.  Preallocated
    * as a performance enhancement to avoid constructing
    * them in multiple places.
    */

   static Box* s_emptys[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*
    * Box that represents the maximum allowed index extents,
    * the "universe" that can be represented.
    */
   static Box* s_universes[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;
};

/**
 * Class BoxIterator is an iterator that provides methods for
 * stepping through the index space associated with a box.  The indices
 * are enumerated in column-major (e.g., Fortran) order.  The iterator
 * should be used as follows:
 * \verbatim
 * Box box;
 * ...
 * Box::iterator bend(box, false);
 * for (Box::iterator b(box, true); b != bend; ++b) {
 *    // use index b of the box
 * }
 * \endverbatim
 * Note that the box iterator may not compile to efficient code, depending
 * on your compiler.  Many compilers are not smart enough to optimize the
 * looping constructs and indexing operations.
 *
 * @see hier::Index
 * @see hier::Box
 */

class BoxIterator
{
public:
   /**
    * Constructor for the box iterator.  The iterator will enumerate the
    * indices in the argument box.
    */
   explicit BoxIterator(
      const Box& box,
      bool begin);

   /**
    * Copy constructor for the box iterator.
    */
   BoxIterator(
      const BoxIterator& iterator);

   /**
    * Assignment operator for the box iterator.
    */
   BoxIterator&
   operator = (
      const BoxIterator& iterator)
   {
      d_index = iterator.d_index;
      d_box = iterator.d_box;
      return *this;
   }

   /**
    * Destructor for the box iterator.
    */
   ~BoxIterator();

   /**
    * Return the current index in the box.  This operation is undefined
    * if the iterator is past the last Index in the box.
    */
   const Index&
   operator * () const
   {
      return d_index;
   }

   /**
    * Return a pointer to the current index in the box.  This operation is
    * undefined if the iterator is past the last Index in the box.
    */
   const Index*
   operator -> () const
   {
      return &d_index;
   }

   /**
    * Post-increment the iterator to point to the next index in the box.
    */
   BoxIterator
   operator ++ (
      int)
   {
      BoxIterator tmp = *this;
      d_index(0)++;
      for (int i = 0; i < (d_index.getDim().getValue() - 1); i++) {
         if (d_index(i) > d_box.upper(i)) {
            d_index(i) = d_box.lower(i);
            d_index(i + 1)++;
         } else
            break;
      }
      return tmp;
   }

   /**
    * Pre-increment the iterator to point to the next index in the box.
    */
   BoxIterator&
   operator ++ ()
   {
      d_index(0)++;
      for (int i = 0; i < (d_index.getDim().getValue() - 1); i++) {
         if (d_index(i) > d_box.upper(i)) {
            d_index(i) = d_box.lower(i);
            d_index(i + 1)++;
         } else
            break;
      }
      return *this;
   }

   /**
    * Test two iterators for equality (same index value).
    */
   bool
   operator == (
      const BoxIterator& iterator) const
   {
      return d_index == iterator.d_index;
   }

   /**
    * Test two iterators for inequality (different index values).
    */
   bool
   operator != (
      const BoxIterator& iterator) const
   {
      return d_index != iterator.d_index;
   }

private:
   /*
    * Unimplemented default constructor.
    */
   BoxIterator();

   Index d_index;

   Box d_box;

};

}
}


#endif
