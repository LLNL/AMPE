/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Dimension class for abstracting dimension
 *
 ************************************************************************/

#ifndef included_tbox_Dimension
#define included_tbox_Dimension

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Utilities.h"

#include <iostream>
#include <limits>

/*
 * Forward declarations, which are questionable with respect to SAMRAI 
 * package ordering.   These are needed since pdat::ArrayData and 
 * hier::IntVector classes need to access private Dimension assignment
 * constructor.
 *
 * It would be good to come up with an alternative to this.
 */
namespace SAMRAI {

namespace hier {
class IntVector;
}

namespace pdat {
template<class TYPE>
class ArrayData;
}

}

namespace SAMRAI {
namespace tbox {

class DatabaseBox;

/**
 * Class Dimension is used to represent the dimension of a SAMRAI
 * object.
 *
 * The maximum dimension is set at compile time using a flag to the
 * configure script.  This is used to allocate arrays in some lower
 * level classes such as IntVector.  If dynamic memory allocation is
 * used the performance impact is significant; a maximum dimension
 * allows for stack based memory allocation in performance critical
 * classes at the expense of wasting storage for objects with
 * dimension less than the maximum dimension.
 *
 * A class is used rather than a simple short or integer to provide
 * enhanced type safety.
 *
 */

class Dimension
{
public:
   /**
    * Constructor for Dimension, object is built using the specified dimension
    *
    * Note that the constructor is "explicit" thus making automatic
    * type conversions from integers impossible.  This is intentionally to
    * avoid unintended conversions.
    *
    * When dimensional assertion checking is active an assert is
    * thrown when dim < 1 or dim > getMaxDimension() value specified when
    * the library is configured (defaults to 3).  dim also cannot be
    * the getInvalidDimension() (the largest unsigned short value).
    *
    */
   explicit Dimension(
      const unsigned short& dim);

   /**
    * Construct a dimension equal to the argument.
    */
   Dimension(
      const Dimension& dimension);

   /**
    * Returns true if Dimension is valid.
    *
    * A valid Dimension != 0; != getInvalidDimension(),
    * and <= getMaxDimension().
    *
    */
   bool
   isValid() const
   {
      return (d_dim != 0) && (d_dim <= Dimension::getMaxDimValue());
   }

   /**
    * Returns true if Dimension is initialized (not set
    * to getInvalidDimension()).
    *
    */
   bool
   isInitialized() const
   {
      return d_dim != Dimension::getInvalidDimValue();
   }

   /**
    * Equality operator.
    */
   bool
   operator == (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim == rhs.d_dim;
   }

   /**
    * Inequality operator.
    */
   bool
   operator != (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim != rhs.d_dim;
   }

   /**
    * Greater than operator.
    */
   bool
   operator > (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim > rhs.d_dim;
   }

   /**
    * Greater than or equal operator.
    */
   bool
   operator >= (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim >= rhs.d_dim;
   }

   /**
    * Less than operator.
    */
   bool
   operator < (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim < rhs.d_dim;
   }

   /**
    * Less than or equal operator.
    */
   bool
   operator <= (
      const Dimension& rhs) const
   {
      TBOX_DIM_ASSERT_CHECK_DIM(*this);
      TBOX_DIM_ASSERT_CHECK_DIM(rhs);
      return d_dim <= rhs.d_dim;
   }

   /**
    * Returns the dimension of the Dimension as an unsigned short.
    *
    * The method is provided to allow sizing of arrays based on the
    * dimension and for iteration.  In general this should not be
    * used for comparisons, the Dimension comparison operations are
    * better suited for that purpose.
    */
   unsigned short
   getValue() const
   {
      return d_dim;
   }

   /**
    * Returns the maximum dimension for the currently compiled library
    * as an unsigned short.
    *
    * When the SAMRAI library is compiled a maximum dimension allowed
    * is specified (the default is 3).  This method is typically used
    * to allocate arrays.
    *
    *  double array[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
    *
    * The value must be >= 1 and < numeric_limits<unsigned short>::max()
    */
   static const unsigned short MAXIMUM_DIMENSION_VALUE =
      SAMRAI_MAXIMUM_DIMENSION;
   static unsigned short
   getMaxDimValue()
   {
      return SAMRAI_MAXIMUM_DIMENSION;
   }

   /**
    * Returns the maximum dimension for the currently compiled library
    * as a Dimension object.
    *
    * When the SAMRAI library is compiled a maximum dimension allowed
    * is specified (the default is 3).  This method is typically used
    * to allocate arrays.
    *
    */
   static const Dimension&
   getMaxDimension()
   {
      static Dimension dim(SAMRAI_MAXIMUM_DIMENSION);
      return dim;
   }

   /**
    * An invalid dimension value as a Dimension object.
    */
   static const Dimension&
   getInvalidDimension()
   {
      static Dimension invalidDim(Dimension::getInvalidDimValue());
      return invalidDim;
   }

   /**
    * An invalid dimension value as an unsigned short.
    *
    * Currently this value is numeric_limits<unsigned short>::max() but
    * use this symbol as it is more readable.
    *
    */
   static unsigned short
   getInvalidDimValue()
   {
      static unsigned short invalid =
         std::numeric_limits<unsigned short>::max();
      return invalid;
   }

   /*
    * Classes that are friends of dimension in order to access th
    * private ctor which builds invalid dimensions.
    *
    * This is obviously not a very good design but so far
    * a better solution has been elusive.   Allowing
    * any code to create invalid dimensions seemed too
    * error prone.
    */
   template<class>
   friend class pdat::ArrayData;
   friend class hier::IntVector;
   friend class DatabaseBox;

   /**
    * Output operator for debugging and error messages.
    */
   friend std::ostream&
   operator << (
      std::ostream& s,
      const Dimension& rhs);

private:
   /**
    * @brief Create an invalid dimension object.
    *
    * This ctor is private to prevent a default constructor call.
    * Currently Dimension objects must always created with a dimension
    * specified for normal code.  Several special classes are allowed
    * and are declared to be friends to access this ctor.
    */
   Dimension();

   /**
    * Assignment operator is private to prevent dimensions
    * from being assigned.  This was done to improve type
    * safety.
    */
   Dimension&
   operator = (
      const Dimension& rhs)
   {
      d_dim = rhs.d_dim;
      return *this;
   }

   unsigned short d_dim;

   static Dimension s_maximum_dimension;
};

}
}

#endif
