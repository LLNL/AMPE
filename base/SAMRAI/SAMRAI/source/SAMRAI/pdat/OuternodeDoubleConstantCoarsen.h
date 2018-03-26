/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Constant averaging operator for node-centered double data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeDoubleConstantCoarsen
#define included_pdat_OuternodeDoubleConstantCoarsen

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <boost/shared_ptr.hpp>
#include <string>

namespace SAMRAI {
namespace pdat {

/*!
 * @brief
 * Class OuternodeDoubleConstantCoarsen implements constant
 * averaging (i.e., injection) for outernode-centered double patch data defined
 * over a  mesh.
 *
 * It is derived from the hier::CoarsenOperator base class.
 * The numerical operations for the averaging use FORTRAN numerical
 * routines.
 *
 * @see hier::CoarsenOperator
 */

class OuternodeDoubleConstantCoarsen:
   public hier::CoarsenOperator
{
public:
   /*!
    * @brief Uninteresting default constructor.
    */
   explicit OuternodeDoubleConstantCoarsen(
      const tbox::Dimension& dim);

   /*!
    * @brief Uninteresting virtual destructor.
    */
   virtual ~OuternodeDoubleConstantCoarsen();

   /*!
    * @brief Give the operator priority.
    *
    * The priority of outernode-centered constant averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int
   getOperatorPriority() const;

   /*!
    * @brief Give the operator stencil width.
    *
    * The stencil width of the constant averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector
   getStencilWidth() const;

   /*!
    * @brief Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the outernode-centered double constant
    * averaging operator.
    *
    * Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    */
   void
   coarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const int dst_component,
      const int src_component,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;
};

}
}

#endif
