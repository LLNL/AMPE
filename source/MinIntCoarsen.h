// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
//
#ifndef included_MinIntCoarsen
#define included_MinIntCoarsen

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/CoarsenOperator.h"

#include <string>

using namespace SAMRAI;

/**
 * Class MinIntCoarsen implements a simple coarsen operator to
 * take the minimum positive integer of cell-centered integer patch
 * data defined over a Cartesian mesh.  Negative values are ignored, and
 * if all fine values are negative, the coarse value will be set to -1.
 *
 * The findCoarsenOperator() operator function returns true if the input
 * variable is cell-centered int, and the string is "MINIMUM_COARSEN".
 *
 * @see hier::CoarsenOperator
 */

class MinIntCoarsen : public hier::CoarsenOperator
{
 public:
   /**
    * Uninteresting default constructor.
    */
   MinIntCoarsen();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~MinIntCoarsen();

   /**
    * Return true if the variable and name string match cell-centered
    * double weighted averaging; otherwise, return false.
    */
   bool findCoarsenOperator(const std::shared_ptr<hier::Variable>& var,
                            const std::string& op_name) const;

   /**
    * Return name string identifier of this coarsening operation.
    */
   const std::string& getOperatorName() const;

   /**
    * The priority of cell-centered double weighted averaging is 0.
    * It will be performed before any user-defined coarsen operations.
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the weighted averaging operator is the vector of
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector getStencilWidth(const tbox::Dimension& dim) const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the cell-centered double weighted
    * averaging operator.  Coarsening is performed on the intersection of
    * the destination patch and the coarse box.  It is assumed that the
    * fine patch contains sufficient data for the stencil width of the
    * coarsening operator.
    */
   void coarsen(hier::Patch& coarse, const hier::Patch& fine,
                const int dst_component, const int src_component,
                const hier::Box& coarse_box,
                const hier::IntVector& ratio) const;

 private:
   std::string d_name_id;
};

#endif
