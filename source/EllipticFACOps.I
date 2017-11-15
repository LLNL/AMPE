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
#define SAMRAI_INLINE_KEYWORD inline

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setDPatchDataId( int id ) {
   d_poisson_spec.setDPatchDataId(id);
   d_D_is_set = true;
}

SAMRAI_INLINE_KEYWORD
int EllipticFACOps::getDPatchDataId() const{
   return d_poisson_spec.getDPatchDataId();
}

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setDConstant( const double scalar ) {
   d_poisson_spec.setDConstant(scalar);
   d_D_is_set = true;
}

SAMRAI_INLINE_KEYWORD
double EllipticFACOps::getDConstant( ) const{
   return d_poisson_spec.getDConstant();
}

SAMRAI_INLINE_KEYWORD
bool EllipticFACOps::dIsConstant( ) const{
   return d_poisson_spec.dIsConstant();
}

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setCPatchDataId( int id ) {
   d_poisson_spec.setCPatchDataId(id);
   d_C_is_set = true;
}

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setCConstant( const double scalar ) {
   if ( scalar == 0.0 ) {
      d_poisson_spec.setCZero();
   }
   else {
      d_poisson_spec.setCConstant(scalar);
   }
   d_C_is_set = true;
}



/*
********************************************************************
* Set the patch data id for the flux.                              *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setFluxId( int flux_id ) {
   d_flux_id = flux_id;
#ifdef DEBUG_CHECK_ASSERTIONS
   checkInputPatchDataIndices();
#endif
}




/*
********************************************************************
* Set the choice for smoothing algorithm.                          *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setSmoothingChoice(
   const std::string &smoothing_choice )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( smoothing_choice != "redblack" ) {
      TBOX_ERROR(d_object_name << ": Bad smoothing choice '"
                 << smoothing_choice
                 << "' in EllipticFACOps::setSmoothingChoice.");
   }
#endif
   d_smoothing_choice = smoothing_choice;
   return;
}


/*
********************************************************************
* Set the choice for the coarse level solver.                      *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setCoarsestLevelSolverChoice( const std::string &choice ) {
#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef HAVE_HYPRE
   if ( choice == "hypre" ) {
      TBOX_ERROR(d_object_name << ": HYPRe library is not available.\n");
   }
#endif
#endif
   if (  choice == "redblack"
      || choice == "hypre" ) {
      d_coarse_solver_choice = choice;
   }
   else {
      TBOX_ERROR(d_object_name << ": Bad coarse level solver choice '"
                 << choice
                 << "' in scapCellPoissonOpsX::setCoarseLevelSolver.");
   }
   return;
}


/*
********************************************************************
* Set the tolerance for the coarse level solver.                   *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setCoarsestLevelSolverTolerance( double tol ) {
   d_coarse_solver_tolerance = tol;
#ifdef HAVE_HYPRE
   if( d_coarse_solver_choice == "hypre" )
      d_hypre_solver.setStoppingCriteria( tol );
#endif
   return;
}


/*
********************************************************************
* Set the tolerance for the coarse level solver.                   *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setCoarsestLevelSolverMaxIterations(
   int max_iterations ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if( max_iterations < 0 ) {
      TBOX_ERROR(d_object_name << ": Invalid number of max iterations\n");
   }
#endif
   d_coarse_solver_max_iterations = max_iterations;
   return;
}


/*
********************************************************************
* Set the coarse-fine discretization method.                       *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setCoarseFineDiscretization(
   const std::string &coarsefine_method ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Cannot change coarse-fine\n"
                 << "discretization method while operator state\n"
                 << "is initialized because that causes a\n"
                 << "corruption in the state.\n");
   }
#endif
   d_cf_discretization = coarsefine_method;
   return;
}


/*
********************************************************************
* Set the prolongation method                                      *
********************************************************************
*/

SAMRAI_INLINE_KEYWORD
void EllipticFACOps::setProlongationMethod(
   const std::string &prolongation_method ) {
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_hierarchy ) {
      TBOX_ERROR(d_object_name << ": Cannot change prolongation method\n"
                 << "while operator state is initialized because that\n"
                 << "causes a corruption in the state.\n");
   }
#endif
   d_prolongation_method = prolongation_method;
   return;
}

