/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for face-centered complex data.
 *
 ************************************************************************/

#ifndef included_pdat_FaceComplexLinearTimeInterpolateOp
#define included_pdat_FaceComplexLinearTimeInterpolateOp

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"

#include <boost/shared_ptr.hpp>
#include <string>

namespace SAMRAI {
namespace pdat {

/**
 * Class FaceComplexLinearTimeInterpolateOp implements standard
 * linear time interpolation for face-centered complex patch data.
 * It is derived from the hier::TimeInterpolateOperator base class.
 * The interpolation uses FORTRAN numerical routines.
 *
 * @see hier::TimeInterpolateOperator
 */

class FaceComplexLinearTimeInterpolateOp:
   public hier::TimeInterpolateOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   FaceComplexLinearTimeInterpolateOp();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~FaceComplexLinearTimeInterpolateOp();

   /**
    * Perform linear time interpolation between two face-centered complex
    * patch data sources and place result in the destination patch data.
    * Time interpolation is performed on the intersection of the destination
    * patch data and the input box.  The time to which data is interpolated
    * is provided by the destination data.
    */
   void
   timeInterpolate(
      hier::PatchData& dst_data,
      const hier::Box& where,
      const hier::PatchData& src_data_old,
      const hier::PatchData& src_data_new) const;

private:
};

}
}
#endif
