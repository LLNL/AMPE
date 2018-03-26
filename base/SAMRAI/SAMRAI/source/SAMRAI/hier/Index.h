/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface for the AMR Index object
 *
 ************************************************************************/

#ifndef included_hier_Index
#define included_hier_Index

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

/**
 * Class Index implements a simple n-dimensional integer vector in the
 * AMR index space.  Index is used as lower and upper bounds when
 * creating a box and also when iterating over the cells in a box.  An
 * index is essentially an integer vector but it carries along the
 * notion of indexing into AMR's abstract index space.
 *
 * @see hier::Box
 * @see hier::BoxIterator
 * @see hier::IntVector
 */

class Index:public IntVector
{
public:
   /**
    * Creates an uninitialized vector.
    */
   explicit Index(
      const tbox::Dimension& dim);

   /**
    * Construct an index with all components equal to the argument.
    */
   Index(
      const tbox::Dimension& dim,
      const int i);

   /**
    * Construct a two-dimensional index with the value (i,j).
    */
   Index(
      const int i,
      const int j);

   /**
    * Construct a three-dimensional index with the value (i,j,k).
    */
   Index(
      const int i,
      const int j,
      const int k);

   /**
    * Construct an n-dimensional index with the values copied
    * from the integer tbox::Array i of size n.
    */
   explicit Index(
      const tbox::Array<int>& i);

   /**
    * The copy constructor creates an index equal to the argument.
    */
   Index(
      const Index& rhs);

   /**
    * Construct an index equal to the argument IntVector.
    */
   explicit Index(
      const IntVector& rhs);

   /**
    * Construct an index equal to the argument array.
    */
   Index(
      const tbox::Dimension& dim,
      const int array[]);

   /**
    * The assignment operator sets the index equal to the argument.
    *
    * An assignment to an uninitialized Index is allowed but assigning
    * from an uninitialized Index will result in an assert.
    */
   Index&
   operator = (
      const Index& rhs)
   {
#ifdef DEBUG_CHECK_ASSERTIONS
      /*
       * Allow assignment of to an uninitialized
       * but do not allow assignment from an
       * uninitialized.
       */
      if (getDim().isValid()) {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      } else {
         TBOX_DIM_ASSERT_CHECK_DIM(rhs.getDim());
      }
#endif
      IntVector::operator = (rhs);
      return *this;
   }

   /**
    * The assignment operator sets the index equal to the argument IntVector.
    */
   Index&
   operator = (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector::operator = (rhs);
      return *this;
   }

   /**
    * The index destructor does nothing interesting.
    */
   virtual ~Index();

   /**
    * Plus-equals operator for an index and an integer vector.
    */
   Index&
   operator += (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for an index and an integer vector.
    */
   Index
   operator + (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      Index tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus-equals operator for an index and an integer.
    */
   Index&
   operator += (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      IntVector::operator += (rhs);
      return *this;
   }

   /**
    * Plus operator for an index and an integer.
    */
   Index
   operator + (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      Index tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for an index and an integer vector.
    */
   Index&
   operator -= (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for an index and an integer vector.
    */
   Index
   operator - (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      Index tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for an index and an integer.
    */
   Index&
   operator -= (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      IntVector::operator -= (rhs);
      return *this;
   }

   /**
    * Minus operator for an index and an integer.
    */
   Index
   operator - (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      Index tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for an index and an integer vector.
    */
   Index&
   operator *= (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for an index and an integer vector.
    */
   Index
   operator * (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      Index tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for an index and an integer.
    */
   Index&
   operator *= (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      IntVector::operator *= (rhs);
      return *this;
   }

   /**
    * Times operator for an index and an integer.
    */
   Index
   operator * (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      Index tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Assign-quotient operator for an index and an integer vector.
    */
   Index&
   operator /= (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector::operator /= (rhs);
      return *this;
   }

   /**
    * Quotient operator for an index and an integer vector.
    */
   Index
   operator / (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      Index tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /**
    * Assign-quotient operator for an index and an integer.
    */
   Index&
   operator /= (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      IntVector::operator /= (rhs);
      return *this;
   }

   /**
    * Quotient operator for an index and an integer.
    */
   Index
   operator / (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(getDim());
      Index tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /*!
    * @brief Coarsen the Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    */
   Index&
   coarsen(
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, ratio);
      for (int d = 0; d < getDim().getValue(); ++d) {
         (*this)(d) = coarsen((*this)(d), ratio(d));
      }
      return *this;
   }

   /*!
    * @brief Return an Index of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getZeroIndex(
      const tbox::Dimension& dim)
   {
      return *(s_zeros[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an Index of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getOneIndex(
      const tbox::Dimension& dim)
   {
      return *(s_ones[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an Index with minimum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMinIndex(
      const tbox::Dimension& dim)
   {
      return *(s_mins[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an Index with maximum index values for the
    * specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const Index&
   getMaxIndex(
      const tbox::Dimension& dim)
   {
      return *(s_maxs[dim.getValue() - 1]);
   }

   /*!
    * @brief Coarsen an Index by a given ratio.
    *
    * For positive indices, this is the same as dividing by the ratio.
    */
   static Index
   coarsen(
      const Index& index,
      const IntVector& ratio)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(index, ratio);
      tbox::Dimension dim(index.getDim());
      Index tmp(dim);
      for (int d = 0; d < dim.getValue(); ++d) {
         tmp(d) = coarsen(index(d), ratio(d));
      }
      return tmp;
   }

private:
   friend class std::vector<Index>;

   /**
    * The default constructor for Index creates an uninitialized index.
    */
   Index();

   static int
   coarsen(
      const int index,
      const int ratio)
   {
      return index < 0 ? (index + 1) / ratio - 1 : index / ratio;
   }

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    *
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    *
    */
   static void
   finalizeCallback();

   static Index* s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static Index* s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static Index* s_maxs[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static Index* s_mins[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
