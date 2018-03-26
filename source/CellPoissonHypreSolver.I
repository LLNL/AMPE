// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
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
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
/*
 * Description:	Hypre solver interface for diffusion-like elliptic problems.
 * adapted from SAMRAI library class
 */
#define SAMRAI_INLINE_KEYWORD inline

SAMRAI_INLINE_KEYWORD
void CellPoissonHypreSolver::setSolnIdDepth( const int depth ) {
   d_soln_depth = depth;
   return;
}

SAMRAI_INLINE_KEYWORD
void CellPoissonHypreSolver::setRhsIdDepth( const int depth ) {
   d_rhs_depth = depth;
   return;
}





/*
********************************************************************
* Specify bc using the default internal bc coefficient object.     *
* Clear up data supporting external bc coefficient setter.         *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void CellPoissonHypreSolver::setBoundaries(
   const std::string& boundary_type,
   const int fluxes,
   const int flags,
   int* bdry_types)
{
   d_physical_bc_simple_case.setBoundaries(boundary_type,
                                           fluxes,
                                           flags,
                                           bdry_types);
   d_physical_bc_coef_strategy = &d_physical_bc_simple_case;
   d_physical_bc_variable.reset();
   return;
}



/*
********************************************************************
* Set the physical boundary condition object.                      *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void CellPoissonHypreSolver::setPhysicalBcCoefObject(
   const solv::RobinBcCoefStrategy *physical_bc_coef_strategy,
   const boost::shared_ptr< hier::Variable > variable )
{
   d_physical_bc_coef_strategy = physical_bc_coef_strategy;
   d_physical_bc_variable = variable;
   return;
}


SAMRAI_INLINE_KEYWORD
void CellPoissonHypreSolver::setStoppingCriteria(
   const double residual_tol)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( residual_tol>=0.0 );
#endif
   d_relative_residual_tol = residual_tol;
}
