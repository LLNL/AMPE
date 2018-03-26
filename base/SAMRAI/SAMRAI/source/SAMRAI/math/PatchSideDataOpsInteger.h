/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataOpsInteger
#define included_math_PatchSideDataOpsInteger

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/math/ArrayDataNormOpsInteger.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/shared_ptr.hpp>
#include <iostream>

namespace SAMRAI {
namespace math {

/**
 * Class PatchSideDataOpsInteger provides a collection of operations
 * that may be used to manipulate integer side-centered patch data.  The
 * operations include basic arithmetic, min, max, etc.  With the assertion
 * of a few basic routines, this class inherits its interface (and
 * thus its functionality) from the base class PatchSideDataBasicOps
 * from which it is derived.
 *
 * A more extensive set of operations is implemented for real (double and
 * float) and complex patch data in the classes PatchSideDataOpsReal
 * and PatchSideDataOpsComplex, repsectively.
 *
 * @see math::PatchSideDataBasicOps
 */

class PatchSideDataOpsInteger:
   public PatchSideDataBasicOps<int>
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchSideDataOpsInteger();

   virtual ~PatchSideDataOpsInteger();

   /**
    * Return the number of data values for the side-centered data object
    * in the given box.  Note that it is assumed that the box refers to
    * the cell-centered index space corresponding to the patch hierarchy.
    */
   int
   numberOfEntries(
      const boost::shared_ptr<pdat::SideData<int> >& data,
      const hier::Box& box) const;

   /**
    * Copy dst data to src data over given box.
    */
   void
   copyData(
      const boost::shared_ptr<pdat::SideData<int> >& dst,
      const boost::shared_ptr<pdat::SideData<int> >& src,
      const hier::Box& box) const;

   /**
    * Swap pointers for patch data objects.  Objects are checked for
    * consistency of depth, box, and ghost box.
    */
   void
   swapData(
      const boost::shared_ptr<hier::Patch>& patch,
      const int data1_id,
      const int data2_id) const;

   /**
    * Print data entries over given box to given output stream.
    */
   void
   printData(
      const boost::shared_ptr<pdat::SideData<int> >& data,
      const hier::Box& box,
      std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void
   setToScalar(
      const boost::shared_ptr<pdat::SideData<int> >& dst,
      const int& alpha,
      const hier::Box& box) const
   {
      TBOX_ASSERT(dst);
      TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);
      dst->fillAll(alpha, box);
   }

   /**
    * Set destination component to absolute value of source component.
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    */
   void
   abs(
      const boost::shared_ptr<pdat::SideData<int> >& dst,
      const boost::shared_ptr<pdat::SideData<int> >& src,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   PatchSideDataOpsInteger(
      const PatchSideDataOpsInteger&);
   void
   operator = (
      const PatchSideDataOpsInteger&);

   ArrayDataNormOpsInteger d_array_ops;

};

}
}

#endif
