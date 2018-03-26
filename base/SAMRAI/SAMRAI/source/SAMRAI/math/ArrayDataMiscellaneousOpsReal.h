/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Miscellaneous templated operations for real array data
 *
 ************************************************************************/

#ifndef included_math_ArrayDataMiscellaneousOpsReal
#define included_math_ArrayDataMiscellaneousOpsReal

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/ArrayData.h"

namespace SAMRAI {
namespace math {

/**
 * Class ArrayDataMiscellaneousOpsReal provides various operations that
 * may be applied to arrays of real (double and float) numerical data
 * values maintained using pdat::ArrayData<DIM> objects.  These operations are
 * sufficiently different from basic arithmetic and norm operations that we
 * chose to implement them in a separate class.  However, as in the case of the  * more common operations, the intent of this class is to provide a single
 * implementation of the operations as they are needed by objects that
 * manipulate standard array-based patch data types (i.e., cell-centered,
 * face-centered, node-centered).  Each operation is implemented in two
 * different ways.  The choice of operation is based on whether control volume
 * information is to be used to weight the contribution of each data entry
 * to the calculation.  The use of control volumes is important when these
 * operations are used in vector kernels where the data resides over multiple
 * levels of spatial resolution in an AMR hierarchy.  The actual index
 * region on which each operation occurs is the intersection of this box
 * and the boxes of all the pdat::ArrayData<DIM> objects involved.
 *
 * Since these operations are used only by the vector kernels for the KINSOL
 * and CVODE solver packages at this time, they are intended to be instantiated
 * for the standard built-in types double and float (since those solvers only
 * treat double and float data).  To extend this class to other data types or
 * to include other operations, the member functions must be specialized or the  * new operations must be added.
 *
 * @see pdat::ArrayData
 */

template<class TYPE>
class ArrayDataMiscellaneousOpsReal
{
public:
   /**
    * Empty constructor and destructor.
    */
   ArrayDataMiscellaneousOpsReal();

   ~ArrayDataMiscellaneousOpsReal<TYPE>();

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the index region, where \f$cvol_i > 0\f$.  Otherwise return 0.
    */
   int
   computeConstrProdPosWithControlVolume(
      const pdat::ArrayData<TYPE>& data1,
      const pdat::ArrayData<TYPE>& data2,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Return 1 if \f$\|data2_i\| > 0\f$ and \f$data1_i * data2_i \leq 0\f$, for
    * any \f$i\f$ in the index region.  Otherwise return 0.
    */
   int
   computeConstrProdPos(
      const pdat::ArrayData<TYPE>& data1,
      const pdat::ArrayData<TYPE>& data2,
      const hier::Box& box) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1\f$
    * if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.
    */
   void
   compareToScalarWithControlVolume(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src,
      const TYPE& alpha,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Set \f$dst_i = 1\f$ if \f$\|src_i\| > \alpha\f$, and \f$dst_i = 0\f$ otherwise.
    */
   void
   compareToScalar(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src,
      const TYPE& alpha,
      const hier::Box& box) const;

   /**
    * Wherever \f$cvol_i > 0\f$ in the index region, set \f$dst_i = 1/src_i\f$ if
    * \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.  If \f$dst_i = 0\f$ anywhere,
    * 0 is the return value.  Otherwise 1 is returned.
    */
   int
   testReciprocalWithControlVolume(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src,
      const pdat::ArrayData<double>& cvol,
      const hier::Box& box) const;

   /**
    * Set \f$dst_i = 1/src_i\f$ if \f$src_i \neq 0\f$, and \f$dst_i = 0\f$ otherwise.
    * If \f$dst_i = 0\f$ anywhere, 0 is the return value.  Otherwise 1 is returned.
    */
   int
   testReciprocal(
      pdat::ArrayData<TYPE>& dst,
      const pdat::ArrayData<TYPE>& src,
      const hier::Box& box) const;

   /*!
    * @brief Compute max of "conditional" quotients of two arrays.
    *
    * Return the maximum of pointwise "conditional" quotients of the numerator
    * and denominator.
    *
    * The "conditional" quotient is defined as |numerator/denominator|
    * where the denominator is nonzero.  Otherwise, it is defined as
    * |numerator|.
    *
    * @b Note: This method is currently intended to support the
    * PETSc-2.1.6 vector wrapper only.  Please do not use it!
    */
   TYPE
   maxPointwiseDivide(
      const pdat::ArrayData<TYPE>& numer,
      const pdat::ArrayData<TYPE>& denom,
      const hier::Box& box) const;

   /*!
    * @brief Compute min of quotients of two arrays.
    *
    * Return the minimum of pointwise quotients of the numerator
    * and denominator.
    *
    * The quotient is defined as (numerator/denominator)
    * where the denominator is nonzero.  When the denominator is zero, the
    * entry is skipped.  If the denominator is always zero, the value of
    * tbox::IEEE::getDBL_MAX() is returned (see @ref SAMRAI::tbox::IEEE).
    *
    * @b Note: This method is currently intended to support the
    * SUNDIALS vector wrapper only.  Please do not use it!
    */
   TYPE
   minPointwiseDivide(
      const pdat::ArrayData<TYPE>& numer,
      const pdat::ArrayData<TYPE>& denom,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   ArrayDataMiscellaneousOpsReal(
      const ArrayDataMiscellaneousOpsReal<TYPE>&);
   void
   operator = (
      const ArrayDataMiscellaneousOpsReal<TYPE>&);

};

}
}

#include "SAMRAI/math/ArrayDataMiscellaneousOpsReal.C"

#endif
