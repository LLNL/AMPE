/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A N-dimensional integer vector
 *
 ************************************************************************/

#ifndef included_hier_IntVector
#define included_hier_IntVector

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"

#include <vector>
#include <iostream>

namespace SAMRAI {

namespace hier {

/**
 * Class IntVector implements a simple N-dimensional integer
 * vector.  This class is the base class for most of the simple indexing
 * classes.
 *
 */

class IntVector
{
public:
   /**
    * Creates an uninitialized vector.
    */
   explicit IntVector(
      const tbox::Dimension& dim);

   /**
    * Construct an integer vector with all components equal to the argument.
    */
   IntVector(
      const tbox::Dimension& dim,
      const int i);

   /**
    * Construct a n-dimensional integer vector with the value with
    * values provided by the array.
    *
    * Dimension inferred from array size.
    */
   explicit IntVector(
      const tbox::Array<int>& a);

   /**
    * Construct a n-dimensional integer vector with the value with
    * values provided by the array.
    *
    */
   IntVector(
      const tbox::Dimension& dim,
      const int array[]);

   /**
    * Construct an integer vector equal to the argument.
    */
   IntVector(
      const IntVector& rhs);

   /**
    * The assignment operator sets the integer vector equal to the argument.
    *
    * An assignment to an uninitialized Index is allowed but assigning
    * from an uninitialized Index will result in an assert.
    */
   IntVector&
   operator = (
      const IntVector& rhs)
   {
      /*
       * Allow assignment of to an uninitialized but do not allow
       * assignment from an uninitialized.
       */
      if (d_dim.isValid()) {
         TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      } else {
         TBOX_DIM_ASSERT_CHECK_DIM(rhs.getDim());
         d_dim = rhs.getDim();
      }
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] = rhs.d_vector[i];
      }
#ifdef DEBUG_INITIALIZE_UNDEFINED
      for (int i = d_dim.getValue(); i < tbox::Dimension::MAXIMUM_DIMENSION_VALUE; i++) {
         d_vector[i] = tbox::MathUtilities<int>::getMin();
      }
#endif
      return *this;
   }

   /**
    * The integer vector destructor does nothing interesting.
    */
   virtual ~IntVector();

   /**
    * Return the specified component of the vector.  No bounds checking.
    */
   int&
   operator [] (
      const int i)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      return d_vector[i];
   }

   /**
    * Return the specified component of the vector as a const integer.
    * No bounds checking.
    */
   const int&
   operator [] (
      const int i) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      return d_vector[i];
   }


   /**
    * Return the specified component of the vector.  No bounds checking.
    */
   int&
   operator () (
      const int i)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      return d_vector[i];
   }

   /**
    * Return the specified component of the vector as a const integer.
    * No bounds checking.
    */
   const int&
   operator () (
      const int i) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      TBOX_ASSERT(i >= 0 && i < d_dim.getValue());
      return d_vector[i];
   }

   /**
    * Plus-equals operator for two integer vectors.
    */
   IntVector&
   operator += (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] += rhs.d_vector[i];
      }
      return *this;
   }

   /**
    * Plus operator for two integer vectors.
    */
   IntVector
   operator + (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector tmp(*this);
      tmp += rhs;
      return tmp;
   }

   /**
    * Plus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator += (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] += rhs;
      }
      return *this;
   }

   /**
    * Plus operator for an integer vector and an integer.
    */
   IntVector
   operator + (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      IntVector tmp = *this;
      tmp += rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for two integer vectors.
    */
   IntVector&
   operator -= (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] -= rhs.d_vector[i];
      }
      return *this;
   }

   /**
    * Minus operator for two integer vectors.
    */
   IntVector
   operator - (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Minus-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator -= (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] -= rhs;
      }
      return *this;
   }

   /**
    * Minus operator for an integer vector and an integer.
    */
   IntVector
   operator - (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      IntVector tmp = *this;
      tmp -= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for two integer vectors.
    */
   IntVector&
   operator *= (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] *= rhs.d_vector[i];
      }
      return *this;
   }

   /**
    * Times operator for two integer vectors.
    */
   IntVector
   operator * (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Times-equals operator for an integer vector and an integer.
    */
   IntVector&
   operator *= (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] *= rhs;
      }
      return *this;
   }

   /**
    * Times operator for an integer vector and an integer.
    */
   IntVector
   operator * (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      IntVector tmp = *this;
      tmp *= rhs;
      return tmp;
   }

   /**
    * Assign-quotient operator for two integer vectors.
    */
   IntVector&
   operator /= (
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] /= rhs.d_vector[i];
      }
      return *this;
   }

   /**
    * Component-wise ceilingDivide quotient (integer divide with rounding up).
    */
   void
   ceilingDivide(
      const IntVector& denominator)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, denominator);
      for (int i = 0; i < getDim().getValue(); i++) {
         /*
          * This is the formula for integer divide, rounding away from
          * zero.  It is meant as an extension of the ceilingDivide quotient of
          * 2 positive integers.
          *
          * The ceilingDivide is the integer divide plus 0, -1 or 1 representing
          * the results of rounding.
          * - Add zero if there's no remainder to round.
          * - Round remainder to 1 if numerator and denominator has same sign.
          * - Round remainder to -1 if numerator and denominator has opposite sign.
          */
         d_vector[i] = (d_vector[i] / denominator[i]) +
            ((d_vector[i] % denominator[i]) ?
               ((d_vector[i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
      }
   }

   /**
    * Component-wise ceilingDivide quotient (integer divide with rounding up).
    */
   static IntVector
   ceilingDivide(
      const IntVector& numerator,
      const IntVector& denominator)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(numerator, denominator);
      IntVector rval(numerator.getDim());
      for (int i = 0; i < numerator.getDim().getValue(); i++) {
         /*
          * This is the formula for integer divide, rounding away from
          * zero.  It is meant as an extension of the ceilingDivide quotient of
          * 2 positive integers.
          *
          * The ceilingDivide is the integer divide plus 0, -1 or 1 representing
          * the results of rounding.
          * - Add zero if there's no remainder to round.
          * - Round remainder to 1 if numerator and denominator has same sign.
          * - Round remainder to -1 if numerator and denominator has opposite sign.
          */
         rval[i] = (numerator[i] / denominator[i]) +
            ((numerator[i] % denominator[i]) ?
               ((numerator[i] > 0) == (denominator[i] > 0) ? 1 : -1) : 0);
      }
      return rval;
   }

   /**
    * Quotient operator for two integer vectors.
    */
   IntVector
   operator / (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      IntVector tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /**
    * Assign-quotient operator for an integer vector and an integer.
    */
   IntVector&
   operator /= (
      const int rhs)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      for (int i = 0; i < d_dim.getValue(); i++) {
         d_vector[i] /= rhs;
      }
      return *this;
   }

   /**
    * Quotient operator for an integer vector and an integer.
    */
   IntVector
   operator / (
      const int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      IntVector tmp = *this;
      tmp /= rhs;
      return tmp;
   }

   /**
    * Unary minus to negate an integer vector.
    */
   IntVector
   operator - () const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      IntVector tmp(d_dim);
      for (int i = 0; i < d_dim.getValue(); i++) {
         tmp.d_vector[i] = -d_vector[i];
      }
      return tmp;
   }

   /**
    * Returns true if all components are equal to a given integer.
    */
   bool
   operator == (
      int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      bool result = true;
      for (int i = 0; result && (i < d_dim.getValue()); i++) {
         result = d_vector[i] == rhs;
      }
      return result;
   }

   /**
    * Returns true if some components are not equal to a given integer.
    */
   bool
   operator != (
      int rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      bool result = true;
      for (int i = 0; result && (i < d_dim.getValue()); i++) {
         result = d_vector[i] != rhs;
      }
      return result;
   }

   /**
    * Returns true if two vector objects are equal.  All components
    * must be the same for equality.
    */
   bool
   operator == (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      bool result = true;
      for (int i = 0; result && (i < d_dim.getValue()); i++) {
         result = result && (d_vector[i] == rhs.d_vector[i]);
      }
      return result;
   }

   /**
    * Returns true if two vector objects are not equal.  Any of
    * the components may be different for inequality.
    */
   bool
   operator != (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      return !(*this == rhs);
   }

   /**
    * Returns true if each integer in vector is less than
    * corresponding integer in comparison vector.
    */
   bool
   operator < (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      bool result = true;
      for (int i = 0; result && (i < d_dim.getValue()); i++) {
         result = result && (d_vector[i] < rhs.d_vector[i]);
      }
      return result;
   }

   /**
    * Returns true if each integer in vector is less or equal to
    * corresponding integer in comparison vector.
    */
   bool
   operator <= (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      bool result = true;
      for (int i = 0; result && (i < d_dim.getValue()); i++) {
         result = result && (d_vector[i] <= rhs.d_vector[i]);
      }
      return result;
   }

   /**
    * Returns true if each integer in vector is greater than
    * corresponding integer in comparison vector.
    */
   bool
   operator > (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      bool result = true;
      for (int i = 0; result && (i < d_dim.getValue()); i++) {
         result = result && (d_vector[i] > rhs.d_vector[i]);
      }
      return result;
   }

   /**
    * Returns true if each integer in vector is greater or equal to
    * corresponding integer in comparison vector.
    */
   bool
   operator >= (
      const IntVector& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      bool result = true;
      for (int i = 0; result && (i < d_dim.getValue()); i++) {
         result = result && (d_vector[i] >= rhs.d_vector[i]);
      }
      return result;
   }

   /**
    * Return the component-wise minimum of two integer vector objects.
    */
   void
   min(
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      for (int i = 0; i < d_dim.getValue(); i++) {
         if (rhs.d_vector[i] < d_vector[i]) {
            d_vector[i] = rhs.d_vector[i];
         }
      }
   }

   /**
    * Return the minimum entry in an integer vector.
    */
   int
   min() const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      int min = d_vector[0];
      for (int i = 1; i < d_dim.getValue(); i++) {
         if (d_vector[i] < min) {
            min = d_vector[i];
         }
      }
      return min;
   }

   /**
    * Return the component-wise maximum of two integer vector objects.
    */
   void
   max(
      const IntVector& rhs)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, rhs);
      for (int i = 0; i < d_dim.getValue(); i++) {
         if (rhs.d_vector[i] > d_vector[i]) {
            d_vector[i] = rhs.d_vector[i];
         }
      }
   }

   /**
    * Return the maximum entry in an integer vector.
    */
   int
   max() const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(d_dim);
      int max = d_vector[0];
      for (int i = 1; i < d_dim.getValue(); i++) {
         if (d_vector[i] > max) {
            max = d_vector[i];
         }
      }
      return max;
   }

   /**
    * Utility function to take the minimum of two integer vector objects.
    */
   static IntVector
   min(
      const IntVector& a,
      const IntVector& b)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(a, b);
      IntVector tmp = a;
      tmp.min(b);
      return tmp;
   }

   /**
    * Utility function to take the maximum of two integer vector objects.
    */
   static IntVector
   max(
      const IntVector& a,
      const IntVector& b)
   {
      TBOX_DIM_ASSERT_CHECK_ARGS2(a, b);
      IntVector tmp = a;
      tmp.max(b);
      return tmp;
   }

   /**
    * Return the product of the entries in the integer vector.
    */
   int
   getProduct() const
   {
      int prod = 1;
      for (int i = 0; i < d_dim.getValue(); i++) {
         prod *= d_vector[i];
      }
      return prod;
   }

   /**
    * Store the object state to the specified database
    * with the provided name.
    *
    */
   virtual void
   putUnregisteredToDatabase(
      tbox::Database& database,
      const std::string& name) const;

   /**
    * Restores the object state from the specified database
    * with the provided name.
    *
    */
   virtual void
   getFromDatabase(
      tbox::Database& database,
      const std::string& name);

   /**
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

   /*!
    * @brief Return an IntVector of zeros of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const IntVector&
   getZero(
      const tbox::Dimension& dim)
   {
      return *(s_zeros[dim.getValue() - 1]);
   }

   /*!
    * @brief Return an IntVector of ones of the specified dimension.
    *
    * Can be used to avoid object creation overheads.
    */
   static const IntVector&
   getOne(
      const tbox::Dimension& dim)
   {
      return *(s_ones[dim.getValue() - 1]);
   }

   /*
    * @brief Sort the given IntVector the smallest to the largest value.
    *
    * Set the ith entry of this to the position of the ith largest
    * value in the given IntVector.
    */
   void
   sortIntVector(
      const IntVector& values);

   /**
    * Read an integer vector from an input stream.  The format for
    * the input is (i0,...,in) for an n-dimensional vector.
    */
   friend std::istream&
   operator >> (
      std::istream& s,
      IntVector& rhs);

   /**
    * Write an integer vector into an output stream.  The format for
    * the output is (i0,...,in) for an n-dimensional vector.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const IntVector& rhs);

   friend class std::vector<IntVector>;

protected:
   /**
    * Default ctor for IntVector is protected to disallow normal use.
    * This is needed by the poorly designed STL container library.
    *
    *
    */
   IntVector();

private:
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
    */
   static void
   finalizeCallback();

   tbox::Dimension d_dim;

   int d_vector[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static IntVector* s_zeros[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   static IntVector* s_ones[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

};

}
}

#endif
