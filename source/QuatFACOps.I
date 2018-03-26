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
 * Created from solv::CellPoissonFACOps.I in the SAMRAI library
 */

#include <cassert>
#define SAMRAI_INLINE_KEYWORD inline

/*
********************************************************************
* Set verbosity                                                    *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setVerbose(bool verbose)
{
   d_verbose = verbose;
}

/*
********************************************************************
* Set the FACPreconditioner object.  This is the FAC solver        *
* object that uses this object.  Having a pointer to it            *
* enables callbacks.                                               *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setSolver(
  const solv::FACPreconditioner *solver ) {
  d_solver = solver;
}

/*
********************************************************************
* Set the physical boundary condition object.                      *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setPhysicalBcCoefObject(
   const solv::RobinBcCoefStrategy *physical_bc_coef )
{
   d_physical_bc_coef = physical_bc_coef;
   d_bc_helper.setCoefImplementation(physical_bc_coef);
}


/*
********************************************************************
* Set the patch data id for the flux.                              *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setFluxId( int flux_id ) {
   d_flux_id = flux_id;
}


/*
********************************************************************
* Set the solver tolerance                                         *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setLevelSolverTolerance( double tol )
{
   d_levelsolver_tolerance = tol;
}

/*
********************************************************************
* Set the solver tolerance for the coarse level                    *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setCoarsestLevelSolverTolerance( double tol ) {
   d_coarse_levelsolver_tolerance = tol;
}


/*
********************************************************************
* Set the maximum number of solver iterations                      *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setLevelSolverMaxIterations(int max_iterations)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if( max_iterations < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid number of max iterations\n");
   }
#endif
   d_levelsolver_max_iterations = max_iterations;
}

/*
********************************************************************
* Set the maximum number of iterations for the coarse level        *
* solver.                                                          *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setCoarsestLevelSolverMaxIterations(
   int max_iterations ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( max_iterations < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid number of max iterations\n");
   }
#endif
   d_coarse_levelsolver_max_iterations = max_iterations;
}


/*
********************************************************************
* Set the coarse-fine discretization method.                       *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setCoarseFineDiscretization(
   const string &coarsefine_method ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Cannot change coarse-fine\n"
                 << "discretization method while operator state\n"
                 << "is initialized because that causes a\n"
                 << "corruption in the state.\n");
   }
#endif
   d_cf_discretization = coarsefine_method;
}


/*
********************************************************************
* Set the prolongation method                                      *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD void
QuatFACOps::setProlongationMethod(
   const string &prolongation_method ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Cannot change prolongation method\n"
                 << "while operator state is initialized because that\n"
                 << "causes a corruption in the state.\n");
   }
#endif
   d_prolongation_method = prolongation_method;
}
