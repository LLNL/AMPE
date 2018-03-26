/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A box structure representing a portion of the AMR index space
 *
 ************************************************************************/

#ifndef included_tbox_DatabaseBox
#define included_tbox_DatabaseBox

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Dimension.h"

#ifndef DatabaseBox_MAX_DIM
#define DatabaseBox_MAX_DIM 3
#else
Macro overloaded : DatabaseBox_MAX_DIM
#endif

namespace SAMRAI {
namespace tbox {

/*!
 * @brief POD data for class DatabaseBox
 *
 * The data in DatabaseBox need to reside in a POD class so that
 * HDF5's HOFFSET macro works.  (According to ANSI C++ standard,
 * it does not have to work with non-POD data.)
 */
struct DatabaseBox_POD {
   int d_dimension;
   int d_lo[DatabaseBox_MAX_DIM];
   int d_hi[DatabaseBox_MAX_DIM];
};

/**
 * Class DatabaseBox represents a one, two, or three dimensional box in the
 * AMR index space.  It is defined by lower and upper bounds given by integer
 * arrays.
 *
 * This box is an auxilliary data structure used by the database routines to
 * manipulate boxes.  This box type removes cyclic dependencies among the
 * database routines (which need a box) and the box (which needs the database
 * routines).  The box classes in the hierarchy package convert this box
 * structure into the standard SAMRAI box class used by the AMR algorithms.
 *
 * @internal This class should have @em NO data except for d_data.
 * See d_data for details.
 */

class DatabaseBox
{
public:
   /**
    * The default constructor creates a zero dimension empty box.
    */
   DatabaseBox();

   /**
    * Create a box of the specified dimension describing the index space
    * between lower and upper.  The dimension argument must be zero, one,
    * two, or three.
    */
   DatabaseBox(
      const Dimension& dim,
      const int * lower,
      const int * upper);

   /**
    * The copy constructor copies the index space of the argument box.
    */
   DatabaseBox(
      const DatabaseBox& box);

   /**
    * The assignment operator copies the index space of the argument box.
    */
   DatabaseBox&
   operator = (
      const DatabaseBox& box);

   /**
    * The destructor does nothing interesting.
    */
   ~DatabaseBox();

   /**
    * Return whether the box is empty.  A box is empty if it has dimension
    * zero or if any of the upper dimensions is less than its corresponding
    * lower dimension.
    */
   bool
   empty() const;

   /**
    * Return the dimension of this object.
    */
   const Dimension&
   getDim() const
   {
      return d_dim;
   }

   void
   setDim(
      const Dimension& dim)
   {
      TBOX_DIM_ASSERT_CHECK_DIM(dim);
      TBOX_ASSERT(dim.getValue() <= DatabaseBox_MAX_DIM);
      d_dim = Dimension(dim);
      d_data.d_dimension = d_dim.getValue();
   }

   /**
    * Return the specified component (non-const) of the lower index of the box.
    */
   int&
   lower(
      const int i)
   {
      TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
      return d_data.d_lo[i];
   }

   /**
    * Return the specified component (non-const) of the upper index of the box.
    */
   int&
   upper(
      const int i)
   {
      TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
      return d_data.d_hi[i];
   }

   /**
    * Return the specified component (const) of the lower index of the box.
    */
   int
   lower(
      const int i) const
   {
      TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
      return d_data.d_lo[i];
   }

   /**
    * Return the specified component (const) of the upper index of the box.
    */
   int
   upper(
      const int i) const
   {
      TBOX_ASSERT((i >= 0) && (i < d_data.d_dimension));
      return d_data.d_hi[i];
   }


   /**
    * Check whether two boxes represent the same portion of index space.
    */
   int
   operator == (
      const DatabaseBox& box) const;

   /**
    * @brief All data members in a POD type.
    *
    * Due to the need to compute offsets for data members and that
    * offsets cannot be computed for non-POD data, we place all
    * data members in a POD struct and own an object of that
    * struct.
    *
    * Data members are public so that the HDFDatabase need not
    * mirror this structure in defining a compound type for HDF.
    */
   DatabaseBox_POD d_data;

   Dimension d_dim;
};

}
}

#endif
