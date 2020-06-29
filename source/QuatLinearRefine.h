// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
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
