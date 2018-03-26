/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataOpsReal
#define included_math_PatchSideDataOpsReal

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/math/PatchSideDataMiscellaneousOpsReal.h"
#include "SAMRAI/math/PatchSideDataNormOpsReal.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/PIO.h"

#include <boost/shared_ptr.hpp>
#include <iostream>

namespace SAMRAI {
namespace math {

/**
 * Class PatchSideDataOpsReal provides a collection of operations
 * to manipulate float and double numerical side-centered patch data.  The
 * operations include basic arithmetic, norms and ordering, and assorted
 * miscellaneous operations.  With the assertion of a few basic routines,
 * this class inherits its interface (and thus its functionality) from the
 * base classes PatchSideDataBasicOps, PatchSideDataNormOpsReal,
 * and PatchSideDataMiscellaneousOpsReal from which it is derived.  The
 * name of each of these base classes is indicative of the set of
 * side-centered patch data operations that it provides.
 *
 * Note that this templated class should only be used to instantiate
 * objects with double or float as the template parameter.  A similar set of
 * operations is implemented for complex and integer patch data in the classes
 * PatchSideDataOpsComplex and PatchSideDataOpsInteger,
 * repsectively.
 *
 * @see math::PatchSideDataBasicOps
 * @see math::PatchSideDataMiscellaneousOpsReal
 * @see math::PatchSideDataNormOpsReal
 */

template<class TYPE>
class PatchSideDataOpsReal:
   public PatchSideDataBasicOps<TYPE>,
   public PatchSideDataMiscellaneousOpsReal<TYPE>,
   public PatchSideDataNormOpsReal<TYPE>
{
public:
   /**
    * Empty constructor and destructor.
    */
   PatchSideDataOpsReal();

   virtual ~PatchSideDataOpsReal<TYPE>();

   /**
    * Copy dst data to src data over given box.
    */
   void
   copyData(
      const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
      const boost::shared_ptr<pdat::SideData<TYPE> >& src,
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
      const boost::shared_ptr<pdat::SideData<TYPE> >& data,
      const hier::Box& box,
      std::ostream& s = tbox::plog) const;

   /**
    * Initialize data to given scalar over given box.
    */
   void
   setToScalar(
      const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
      const TYPE& alpha,
      const hier::Box& box) const;

private:
   // The following are not implemented:
   PatchSideDataOpsReal(
      const PatchSideDataOpsReal<TYPE>&);
   void
   operator = (
      const PatchSideDataOpsReal<TYPE>&);

};

}
}

#include "SAMRAI/math/PatchSideDataOpsReal.C"

#endif
