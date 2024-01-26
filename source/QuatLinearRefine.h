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
#ifndef included_QuatLinearRefine
#define included_QuatLinearRefine

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/BoxOverlap.h"

#include <string>
using namespace SAMRAI;

/**
 * Class QuatLinearRefine implements linear
 * interpolation for node-centered quaternion patch data defined over a
 * Cartesian mesh.  It is derived from the hier::RefineOperator base class. The
 * numerical operations for interpolation use FORTRAN numerical routines.
 *
 * The findRefineOperator() operator function returns true if the input
 * variable is node-centered double, and the string is "MY_LINEAR_REFINE".
 *
 * @see hier::RefineOperator
 */

class QuatLinearRefine : public hier::RefineOperator
{
 public:
   /**
    * Uninteresting default constructor.
    */
   QuatLinearRefine(const int quat_symm_rotation_id = 0);

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~QuatLinearRefine();

   /**
    * Return true if the variable and name string match node-centered
    * double linear interpolation; otherwise, return false.
    */
   bool findRefineOperator(const std::shared_ptr<hier::Variable>& var,
                           const std::string& op_name) const;

   /**
    * Return name string identifier of this refinement operator.
    */
   const std::string& getOperatorName() const;

   /**
    * The priority of node-centered double linear interpolation is 0.
    * It will be performed before any user-defined interpolation operations.
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the linear interpolation operator is the vector
    * of zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector getStencilWidth(const tbox::Dimension& dim) const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the node-centered double linear
    * interpolation operator.  Interpolation is performed on the intersection
    * of the destination patch and the fine box.   It is assumed that the
    * coarse patch contains sufficient data for the stencil width of the
    * refinement operator.
    */
   void refine(hier::Patch& fine, const hier::Patch& coarse,
               const int dst_component, const int src_component,
               const hier::BoxOverlap& fine_box,
               const hier::IntVector& ratio) const;

   void refine(hier::Patch& fine, const hier::Patch& coarse,
               const int dst_component, const int src_component,
               const hier::Box& fine_box, const hier::IntVector& ratio) const;

 private:
   std::string d_name_id;
   int d_quat_symm_rotation_id;
};

#endif
