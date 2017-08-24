/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Basic templated cell-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataBasicOps
#define included_math_PatchCellDataBasicOps

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/ArrayDataBasicOps.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Complex.h"

#include "boost/shared_ptr.hpp"

namespace SAMRAI {
namespace math {

/**
 * Class PatchCellDataBasicOps provides access to a collection
 * of basic numerical operations that may be applied to numerical cell-
 * centered patch data.  These operations include simple arithmetic
 * operations as well as min and max, etc.  This class provides a single
 * implementation of these operations that may be used to manipulate any
 * cell-centered patch data object.   Each member function accepts a box
 * argument indicating the
 * region of index space on which the operation should be performed.  The
 * operation will be performed on the intersection of this box and those
 * boxes corresponding to the patch data objects involved.
 *
 * These operations typically apply only to the numerical standard built-in
 * types, such as double, float, and int, and the complex type (which may or
 * may not be a built-in type depending on the C++ compiler).  Thus, this
 * templated class should only be used to instantiate objects with those
 * types as the template parameter.  None of the operations are implemented
 * for any other type.
 *
 * @see ArrayDataBasicOps
 */

template<class TYPE>
class PatchCellDataBasicOps
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchCellDataBasicOps();

   virtual ~PatchCellDataBasicOps<TYPE>();

   /**
    * Set dst = alpha * src, elementwise.
    *
    * @pre dst && src
    * @pre (dst->getDim() == src->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   scale(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const TYPE& alpha,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src,
      const hier::Box& box) const;

   /**
    * Set dst = src + alpha, elementwise.
    *
    * @pre dst && src
    * @pre (dst->getDim() == src->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   addScalar(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src,
      const TYPE& alpha,
      const hier::Box& box) const;

   /**
    * Set dst = src1 + src2, elementwise.
    *
    * @pre dst && src1 && src2
    * @pre (dst->getDim() == src1->getDim()) &&
    *      (dst->getDim() == src2->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   add(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src1,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src2,
      const hier::Box& box) const;

   /**
    * Set dst = src1 - src2, elementwise.
    *
    * @pre dst && src1 && src2
    * @pre (dst->getDim() == src1->getDim()) &&
    *      (dst->getDim() == src2->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   subtract(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src1,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src2,
      const hier::Box& box) const;

   /**
    * Set dst = src1 * src2, elementwise.
    *
    * @pre dst && src1 && src2
    * @pre (dst->getDim() == src1->getDim()) &&
    *      (dst->getDim() == src2->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   multiply(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src1,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src2,
      const hier::Box& box) const;

   /**
    * Set dst = src1 / src2, elementwise.  No check for division by zero.
    *
    * @pre dst && src1 && src2
    * @pre (dst->getDim() == src1->getDim()) &&
    *      (dst->getDim() == src2->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   divide(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src1,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src2,
      const hier::Box& box) const;

   /**
    * Set dst = 1 / src, elementwise.  No check for division by zero.
    *
    * @pre dst && src
    * @pre (dst->getDim() == src->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   reciprocal(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src,
      const hier::Box& box) const;

   /**
    * Set dst = alpha * src1 + beta * src2, elementwise.
    *
    * @pre dst && src1 && src2
    * @pre (dst->getDim() == src1->getDim()) &&
    *      (dst->getDim() == src2->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   linearSum(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const TYPE& alpha,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src1,
      const TYPE& beta,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src2,
      const hier::Box& box) const;

   /**
    * Set dst = alpha * src1 + src2, elementwise.
    *
    * @pre dst && src1 && src2
    * @pre (dst->getDim() == src1->getDim()) &&
    *      (dst->getDim() == src2->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   axpy(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const TYPE& alpha,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src1,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src2,
      const hier::Box& box) const;

   /**
    * Set dst = alpha * src1 - src2, elementwise.
    *
    * @pre dst && src1 && src2
    * @pre (dst->getDim() == src1->getDim()) &&
    *      (dst->getDim() == src2->getDim()) && (dst->getDim() == box.getDim())
    */
   void
   axmy(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const TYPE& alpha,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src1,
      const boost::shared_ptr<pdat::CellData<TYPE> >& src2,
      const hier::Box& box) const;

   /**
    * Return the minimum patch data component entry  When the data is
    * complex, the result is the data element with the smallest norm.
    *
    * @pre data
    * @pre data->getDim() == box.getDim()
    */
   TYPE
   min(
      const boost::shared_ptr<pdat::CellData<TYPE> >& data,
      const hier::Box& box) const;

   /**
    * Return the maximum patch data component entry  When the data is
    * complex, the result is the data element with the largest norm.
    *
    * @pre data
    * @pre data->getDim() == box.getDim()
    */
   TYPE
   max(
      const boost::shared_ptr<pdat::CellData<TYPE> >& data,
      const hier::Box& box) const;

   /**
    * Set patch data to random values.  See the operations in the
    * ArrayDataBasicOps class for details on the generation
    * of the random values for each data type.
    *
    * @pre dst
    * @pre dst->getDim() == box.getDim()
    */
   void
   setRandomValues(
      const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
      const TYPE& width,
      const TYPE& low,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   PatchCellDataBasicOps(
      const PatchCellDataBasicOps&);
   PatchCellDataBasicOps&
   operator = (
      const PatchCellDataBasicOps&);

   ArrayDataBasicOps<TYPE> d_array_ops;

};

}
}

#include "SAMRAI/math/PatchCellDataBasicOps.C"

#endif
