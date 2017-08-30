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
#include "CPODESSolver.h"

#ifdef HAVE_SUNDIALS

#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"


// CPODES includes
#ifndef included_cpodes_h
#define included_cpodes_h
extern "C" {
#include "cpodes.h"
}
#endif

#ifndef included_cpodes_impl_h
#define included_cpodes_impl_h
extern "C" {
#include "cpodes/cpodes_impl.h"
}
#endif

#ifndef included_cpodes_spils_impl_h
#define included_cpodes_spils_impl_h
extern "C" {
#include "cpodes/cpodes_spils_impl.h"
}
#endif

#ifndef included_cvspgmr_h
#define included_cvspgmr_h
extern "C" {
#include "cpodes/cpodes_spgmr.h"
}
#endif

#include <cassert>
using namespace solv;

const int CPODESSolver::STAT_OUTPUT_BUFFER_SIZE = 256;

/*
*************************************************************************
*                                                                       *
* Static member functions that provide linkage with CPODES package.      *
* See header file for CPODESAbstractFunctions for more information. *
*                                                                       *
*************************************************************************
*/

int CPODESSolver::CPODESRHSFuncEval( realtype t,
                                     N_Vector y,
                                     N_Vector y_dot,
                                     void* my_solver,
                                     int fd_flag)
{
   return ((CPODESSolver*)my_solver)->getCPODESFunctions()->
      evaluateRHSFunction(t, SABSVEC_CAST(y), SABSVEC_CAST(y_dot), fd_flag);
}

int CPODESSolver::CPODESProjEval( realtype t,
                                  N_Vector y,
                                  N_Vector corr,
                                  realtype epsProj,
                                  N_Vector err,
                                  void* my_solver)
{
   return ((CPODESSolver*)my_solver)->getCPODESFunctions()->
      applyProjection(t, SABSVEC_CAST(y), SABSVEC_CAST(corr), epsProj, SABSVEC_CAST(err));
}

/*
*************************************************************************
*                                                                       *
* Static member functions that provide linkage with CPSpgmr package.        *
*                                                                       *
*************************************************************************
*/

int CPODESSolver::CPSpgmrPrecondSet(realtype t,
                                   N_Vector y,
                                   N_Vector fy,
                                   int jok,
                                   booleantype *jcurPtr,
                                   realtype gamma,
                                   void *my_solver,
                                   N_Vector vtemp1,
                                   N_Vector vtemp2,
                                   N_Vector vtemp3)
{
   
   int success;
   success = ((CPODESSolver*)my_solver)->getCPODESFunctions()->
             CPSpgmrPrecondSet(t,
                               SABSVEC_CAST(y),
                               SABSVEC_CAST(fy),
                               jok,
                               jcurPtr,
                               gamma,
                               SABSVEC_CAST(vtemp1),
                               SABSVEC_CAST(vtemp2),
                               SABSVEC_CAST(vtemp3));

   return (success);
}

int CPODESSolver::CPSpgmrPrecondSolve(realtype t,
                                     N_Vector y,
                                     N_Vector fy,
                                     N_Vector r,
                                     N_Vector z,
                                     realtype gamma,
                                     realtype delta,
                                     int lr,
                                     void *my_solver,
                                     N_Vector vtemp)
{
   int success;
   success = ((CPODESSolver*)my_solver)->getCPODESFunctions()->
             CPSpgmrPrecondSolve(t,
                                 SABSVEC_CAST(y),
                                 SABSVEC_CAST(fy),
                                 SABSVEC_CAST(r),
                                 SABSVEC_CAST(z),
                                 gamma,
                                 delta,
                                 lr,
                                 SABSVEC_CAST(vtemp));

   return (success);
}


/*
*************************************************************************
*                                                                       *
* CPODESSolver constructor and destructor.                         *
*                                                                       *
*************************************************************************
*/
CPODESSolver::CPODESSolver(
   const std::string& object_name,
   CPODESAbstractFunctions* my_functions,
   const bool uses_preconditioner)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(!(my_functions == 0));

   d_object_name         = object_name;
   d_cpode_functions      = my_functions;
   d_uses_preconditioner = uses_preconditioner;

   d_solution_vector = 0;
   d_solution_deriv_vector = (solv::SundialsAbstractVector*)NULL;
   d_weight_vector = (solv::SundialsAbstractVector*)NULL;

   d_print_diagnostics = false;
   d_max_precond_steps = 20;   // the CPODES default

   /*
    * Set default parameters to safe values or to CPODES/CPSpgmr defaults.
    */

   /*
    * CPODES memory record and log file.
    */
   d_cpode_mem           = NULL;
   d_cpode_log_file      = NULL;
   d_cpode_log_file_name = "cpodes.log";
 
   /*
    * ODE parameters.
    */ 
   d_t_0 = 0.0;
   d_user_t_f = 0.0;
   d_actual_t_f = 0.0;
   d_ic_vector = ((solv::SundialsAbstractVector*)NULL);
 
   /*
    * ODE integration parameters.
    */

   setLinearMultistepMethod(CP_BDF);
   setIterationType(CP_FUNCTIONAL);
   setToleranceType("scalar");
   setRelativeTolerance(0.0);
   setAbsoluteTolerance(0.0);
   d_absolute_tolerance_vector = (solv::SundialsAbstractVector*)NULL;
   setSteppingMethod(CP_NORMAL);

   d_max_order = -1;
   d_max_num_internal_steps = -1;
   d_max_num_warnings = -1;
   d_init_step_size = -1;
   d_max_step_size = -1;
   d_min_step_size = -1;

   /*
    * Newton solvers parameters.
    */
   d_nonlin_conv_coef = 0.1;

   /* 
    * CPSpgmr parameters.  
    * 
    * Note that when the maximum krylov dimension and CPSpgmr
    * tolerance scale factor are set to 0, CPSpgmr uses its
    * internal default values.  These are described in the header for 
    * this class.
    */
   setPreconditioningType(PREC_NONE);
   setGramSchmidtType(MODIFIED_GS);
   setMaxKrylovDimension(0);
   setCPSpgmrToleranceScaleFactor(0);

   d_CPODE_needs_initialization = true;
}

CPODESSolver::~CPODESSolver()
{
   if (d_weight_vector != NULL) {
     d_weight_vector->freeVector();
   }
   if (d_solution_deriv_vector != NULL) {
     d_solution_deriv_vector->freeVector();
   }
   if (d_cpode_log_file) fclose(d_cpode_log_file);
   if (d_cpode_mem) CPodeFree(&d_cpode_mem);
}

/*
*************************************************************************
*                                                                       *
* Functions to initialize linear solver and reset CPODES structure.      *
*                                                                       *
*************************************************************************
*/

void CPODESSolver::initialize(solv::SundialsAbstractVector* solution)
{
//#ifdef DEBUG_CHECK_ASSERTIONS
//   TBOX_ASSERT(!(solution == (solv::SundialsAbstractVector*)NULL));
//   TBOX_ASSERT(d_solution_vector == (solv::SundialsAbstractVector*)NULL);
//   double l1norm=solution->L1Norm();
//   assert( l1norm==l1norm );
//#endif
   d_solution_vector = solution;
   d_CPODE_needs_initialization = true;
   initializeCPODES();
}

void CPODESSolver::initializeCPODES() 
{
   TBOX_ASSERT(!(d_solution_vector == (solv::SundialsAbstractVector*)NULL));

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   if (d_CPODE_needs_initialization) {
  
      // Allocate solution derivative vector

      if (d_solution_deriv_vector != NULL) {
        d_solution_deriv_vector->freeVector();
      }

      d_solution_deriv_vector = d_solution_vector->makeNewVector();

      if (d_weight_vector != NULL) {
        d_weight_vector->freeVector();
      }

      d_weight_vector = d_solution_vector->makeNewVector();

      /*
       * Set CPODES log file.
       */ 
      if (d_cpode_log_file) {
         fclose(d_cpode_log_file);
      }
      d_cpode_log_file = fopen(d_cpode_log_file_name.c_str(), "w");

      /*
       * Make sure that either the relative tolerance or the
       * absolute tolerance has been set to a nonzero value.
       */
      bool tolerance_error = false;
      if (d_use_scalar_absolute_tolerance) {
         if ( (d_relative_tolerance == 0.0) &&
              (d_absolute_tolerance_scalar == 0.0) ) { 
            tolerance_error = true;
         }
      } else {
         if ( (d_relative_tolerance == 0.0) &&
             (d_absolute_tolerance_vector->maxNorm() == 0.0) ) {
            tolerance_error = true;
         }
      }

      if (tolerance_error && d_cpode_log_file) {
         fprintf(d_cpode_log_file, 
                 "%s: Both relative and absolute tolerance have value 0.0", 
                 d_object_name.c_str());
      }

      /*
       * CPODES function pointer.
       */
      CPRhsFn RHSFunc = CPODESSolver::CPODESRHSFuncEval;

      /*
       * Free previously allocated Cpodes memory.  Note that the
       * CPReInit() function is not used since the d_neq variable
       * might have been changed from the previous initializeCPODES()
       * call.
       */
      if (d_cpode_mem) CPodeFree(&d_cpode_mem);

      /*
       * Allocate main memory for CPODES package.
       */

      d_cpode_mem = CPodeCreate(d_linear_multistep_method, d_iteration_type);

      N_Vector y0 = ((void *)d_ic_vector != NULL)? d_ic_vector -> getNVector(): NULL;

      int ierr = CPodeInitExpl(d_cpode_mem, RHSFunc, d_t_0, y0);
      CPODE_SAMRAI_ERROR(ierr);

      /*
       * Set tolerance parameters based on type
       */
      if(d_tolerance_type == "vector") {
        ierr = CPodeSVtolerances(d_cpode_mem, d_relative_tolerance,
                                 d_absolute_tolerance_vector -> getNVector());
      } else {
        ierr = CPodeSStolerances(d_cpode_mem, d_relative_tolerance, d_absolute_tolerance_scalar);
      }
      CPODE_SAMRAI_ERROR(ierr);

      ierr = CPodeSetUserData(d_cpode_mem, this);
      CPODE_SAMRAI_ERROR(ierr);

      /*
       * Set the constraint function
       */
      CPProjFn ProjFunc = CPODESSolver::CPODESProjEval;

      ierr = CPodeProjDefine(d_cpode_mem, ProjFunc);
      CPODE_SAMRAI_ERROR(ierr);

      if (d_print_diagnostics) {
         const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
         CPodePrintDiagnostics(d_cpode_mem, mpi.getRank()==0);
      }

      CPodeSetMaxPrecondSteps(d_cpode_mem, d_max_precond_steps);

      /*
       * If the iteration type is set to NEWTON, then set the Newton 
       * convergence factor and initialize the CPSpgmr linear solver.
       */
      if (d_iteration_type == CP_NEWTON) {

        ierr = CPodeSetNonlinConvCoef(d_cpode_mem, d_nonlin_conv_coef);
        CPODE_SAMRAI_ERROR(ierr);

        ierr = CPSpgmr(d_cpode_mem,
                       d_precondition_type, 
                       d_max_krylov_dim);
        CPODE_SAMRAI_ERROR(ierr);
        
        if( !(d_max_order < 1) ) {
          ierr = CPodeSetMaxOrd(d_cpode_mem, d_max_order);
          CPODE_SAMRAI_ERROR(ierr);
        }

        /*
         * Setup CPSpgmr function pointers.
         */
        CPSpilsPrecSetupExplFn precond_set   = NULL;
        CPSpilsPrecSolveExplFn precond_solve = NULL;
        
        if (d_uses_preconditioner) {
          precond_set          = CPODESSolver::CPSpgmrPrecondSet;
          precond_solve = CPODESSolver::CPSpgmrPrecondSolve;
          CPSpilsSetPrecFnExpl(d_cpode_mem, precond_set, precond_solve);
        } 
        
        // Set convergence tolerance
        CPSpilsSetDelt(d_cpode_mem, d_tol_scale_factor);

        // Set GS type
        CPSpilsSetGSType(d_cpode_mem, d_gram_schmidt_type);

        if( !(d_max_num_internal_steps < 0) ) {
          ierr = CPodeSetMaxNumSteps(d_cpode_mem, d_max_num_internal_steps);
          CPODE_SAMRAI_ERROR(ierr);
        }

        if( !(d_max_num_warnings < 0) ) {
          ierr = CPodeSetMaxHnilWarns(d_cpode_mem, d_max_num_warnings);
          CPODE_SAMRAI_ERROR(ierr);
        }

        if( !(d_init_step_size < 0) ) {
          ierr = CPodeSetInitStep(d_cpode_mem, d_init_step_size);
          CPODE_SAMRAI_ERROR(ierr);
        }

        if( !(d_max_step_size < 0) ) {
          ierr = CPodeSetMaxStep(d_cpode_mem, d_max_step_size);
          CPODE_SAMRAI_ERROR(ierr);
        }

        if( !(d_min_step_size < 0) ) {
          ierr = CPodeSetMinStep(d_cpode_mem, d_min_step_size);
          CPODE_SAMRAI_ERROR(ierr);
        }
      }

      d_CPODE_needs_initialization = (y0 == NULL);

   } // if no need to initialize CPODES, function does nothing

}


void CPODESSolver::reinitializeAfterRegrid() 
{
   TBOX_ASSERT(!(d_solution_vector == (solv::SundialsAbstractVector*)NULL));

// Disable Intel warning on real comparison
#ifdef __INTEL_COMPILER
#pragma warning (disable:1572)
#endif

   // Allocate solution derivative vector

   if (d_solution_deriv_vector != NULL) {
      d_solution_deriv_vector->freeVector();
   }

   d_solution_deriv_vector = d_solution_vector->makeNewVector();

   if (d_weight_vector != NULL) {
      d_weight_vector->freeVector();
   }

   d_weight_vector = d_solution_vector->makeNewVector();
}


/*
*************************************************************************
*                                                                       *
* Integrate system of ODEs to d_t_f.  If necessary, re-initialize       *
* CPODES.                                                                *
*                                                                       *
*************************************************************************
*/
int CPODESSolver::solve()
{
   int retval = CP_SUCCESS;

   initializeCPODES();

   /* 
    * Check to make sure that user specified final value for t
    * is greater than initial value for t.
    */
   TBOX_ASSERT( d_user_t_f > d_t_0);

   /*
    * See cpodes.h header file for definition of return types.
    */
   retval = CPode(d_cpode_mem,
                  d_user_t_f,
                  &d_actual_t_f,
                  d_solution_vector -> getNVector(),
                  d_solution_deriv_vector -> getNVector(),
                  d_stepping_method);                  

   return( retval );
}

/*
*************************************************************************
*                                                                       *
* Setting CPODES log file name and print flag for CPODES statistics.    *
*                                                                       *
*************************************************************************
*/

void CPODESSolver::setLogFileData(
   const std::string& log_fname)
{
   TBOX_ASSERT(!log_fname.empty());

   if ( !(log_fname == d_cpode_log_file_name) ) {
      d_cpode_log_file_name = log_fname;
      d_CPODE_needs_initialization = true;
   }
}

/*
*************************************************************************
*                                                                       *
* Accessor functions for user-defined function and linear solver.       *
*                                                                       *
*************************************************************************
*/

void CPODESSolver::setCPODESFunctions(
   CPODESAbstractFunctions* my_functions,
   const bool uses_preconditioner)
{
   TBOX_ASSERT(!(my_functions == (CPODESAbstractFunctions*)NULL));

   d_cpode_functions = my_functions;
   d_uses_preconditioner = uses_preconditioner;
   d_CPODE_needs_initialization = true;
}

CPODESAbstractFunctions* CPODESSolver::getCPODESFunctions() const
{
   return(d_cpode_functions);
}
/*
*************************************************************************
*                                                                       *
* Accessor functions for CPODES integration parameters.                  *
*                                                                       *
*************************************************************************
*/
void CPODESSolver::setLinearMultistepMethod(int linear_multistep_method)
{
   TBOX_ASSERT( (linear_multistep_method == CP_ADAMS) || 
           (linear_multistep_method == CP_BDF) );
   d_linear_multistep_method = linear_multistep_method;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setIterationType(int iteration_type)
{
   TBOX_ASSERT( (iteration_type == CP_FUNCTIONAL) || 
           (iteration_type == CP_NEWTON) );
   d_iteration_type = iteration_type;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setToleranceType(std::string tolerance_type)
{
   TBOX_ASSERT( (tolerance_type == "scalar") || 
           (tolerance_type == "vector") );

   d_tolerance_type = tolerance_type;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setRelativeTolerance(double relative_tolerance)
{
   TBOX_ASSERT(relative_tolerance >= 0.0);

   d_relative_tolerance = relative_tolerance;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setAbsoluteTolerance(double absolute_tolerance)
{
   TBOX_ASSERT(absolute_tolerance >= 0.0);

   d_absolute_tolerance_scalar = absolute_tolerance;
   d_use_scalar_absolute_tolerance = true;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setAbsoluteTolerance(
   solv::SundialsAbstractVector* absolute_tolerance)
{
   TBOX_ASSERT( !(absolute_tolerance == (solv::SundialsAbstractVector*)NULL) );
   TBOX_ASSERT( absolute_tolerance->vecMin() >= 0.0 ); 

   d_absolute_tolerance_vector = absolute_tolerance;
   d_use_scalar_absolute_tolerance = false;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setSteppingMethod(int stepping_method)
{
   TBOX_ASSERT( (stepping_method == CP_NORMAL) || 
           (stepping_method == CP_ONE_STEP) );
#
   d_stepping_method = stepping_method;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setInitialValueOfIndependentVariable(double t_0)
{
   d_t_0 = t_0;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setFinalValueOfIndependentVariable(double t_f,
   bool cpode_needs_initialization)
{
   d_user_t_f = t_f;
   d_CPODE_needs_initialization = cpode_needs_initialization;
}

void CPODESSolver::setInitialConditionVector(
   solv::SundialsAbstractVector* ic_vector)
{
   TBOX_ASSERT( !(ic_vector == (solv::SundialsAbstractVector*)NULL) );

   d_ic_vector = ic_vector;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setMaximumLinearMultistepMethodOrder(
   int max_order)
{
   TBOX_ASSERT( max_order >= 0);

   d_max_order = max_order;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setMaximumNumberOfInternalSteps(
   int max_num_internal_steps)
{
   TBOX_ASSERT( max_num_internal_steps >= 0 );

   d_max_num_internal_steps = max_num_internal_steps;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setMaximumNumberOfNilStepWarnings(
   int max_num_warnings)
{
   TBOX_ASSERT( max_num_warnings >= 0 );

   d_max_num_warnings = max_num_warnings;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setInitialStepSize(double init_step_size)
{
   TBOX_ASSERT( init_step_size >= 0.0 );

   d_init_step_size = init_step_size;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setMaximumAbsoluteStepSize(double max_step_size)
{
   TBOX_ASSERT( max_step_size >= 0.0 );

   d_max_step_size = max_step_size;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setMinimumAbsoluteStepSize(double min_step_size)
{
   TBOX_ASSERT( min_step_size >= 0.0 );

   d_min_step_size = min_step_size;
   d_CPODE_needs_initialization = true;
}

/*
*************************************************************************
*                                                                       *
* Accessor functions for CPSpgmr parameters.                            *
*                                                                       *
*************************************************************************
*/
void CPODESSolver::setPreconditioningType(int precondition_type)
{
   TBOX_ASSERT( (precondition_type == PREC_NONE) ||
           (precondition_type == PREC_LEFT) ||
           (precondition_type == PREC_RIGHT) ||
           (precondition_type == PREC_BOTH) );

   d_precondition_type = precondition_type;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setGramSchmidtType(int gs_type)
{
   TBOX_ASSERT( (gs_type == CLASSICAL_GS) ||
           (gs_type == MODIFIED_GS) );

   d_gram_schmidt_type = gs_type;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setMaxKrylovDimension(int max_krylov_dim)
{
   TBOX_ASSERT( max_krylov_dim >= 0 );

   d_max_krylov_dim = max_krylov_dim;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setCPSpgmrToleranceScaleFactor(double tol_scale_factor)
{
   TBOX_ASSERT( tol_scale_factor >= 0 );

   d_tol_scale_factor = tol_scale_factor;
   d_CPODE_needs_initialization = true;
}

void CPODESSolver::setCPNonlinConvCoef(double coef)
{
   TBOX_ASSERT( coef > 0. );

   d_nonlin_conv_coef = coef;
   d_CPODE_needs_initialization = true;
}

/*
*************************************************************************
*                                                                       *
* Accessor functions for results of CPODES integration step.             *
*                                                                       *
*************************************************************************
*/
solv::SundialsAbstractVector* CPODESSolver::getSolutionVector() const
{
   return (d_solution_vector);
}

solv::SundialsAbstractVector* CPODESSolver::getWeightVector()
{
   int ierr = CPodeGetErrWeights(d_cpode_mem, d_weight_vector->getNVector());
   CPODE_SAMRAI_ERROR(ierr);

   return (d_weight_vector);
}

int CPODESSolver::getDkyVector(double t, int k,
   solv::SundialsAbstractVector* dky) const
{
   int return_code;

   return_code = CPodeGetDky(d_cpode_mem, t, k, dky->getNVector());

   return (return_code);
}

double CPODESSolver::getActualFinalValueOfIndependentVariable() const
{
   return (d_actual_t_f);
}


/*
*************************************************************************
*                                                                       *
* Access methods for CPODES statistics.                                  *
*                                                                       *
*************************************************************************
*/

void CPODESSolver::printStatistics(std::ostream& os) const
{
   printCPODESStatistics(os);
   printCPSpgmrStatistics(os);
}

void CPODESSolver::printCPODESStatistics(std::ostream& os) const
{

   char buf[STAT_OUTPUT_BUFFER_SIZE];
   
   os << "\nCPODESSolver: CPODES statistics... " << std::endl;
   
   sprintf(buf, "lenrw           = %5d     leniw            = %5d\n", 
           getCPODESMemoryUsageForDoubles(), 
           getCPODESMemoryUsageForIntegers());
   os << buf;
   sprintf(buf, "nst             = %5d     nfe              = %5d\n", 
           getNumberOfInternalStepsTaken(),
           getNumberOfRHSFunctionCalls());
   os << buf;
   sprintf(buf, "nni             = %5d     nsetups          = %5d\n", 
           getNumberOfNewtonIterations(),
           getNumberOfLinearSolverSetupCalls());
   os << buf;
   sprintf(buf, "netf            = %5d     ncfn             = %5d\n", 
           getNumberOfLocalErrorTestFailures(),
           getNumberOfNonlinearConvergenceFailures());
   os << buf;
   sprintf(buf, "qu              = %5d     qcur             = %5d\n", 
           getOrderUsedDuringLastInternalStep(),
           getOrderToBeUsedDuringNextInternalStep());
   os << buf;
   sprintf(buf, "\nhu              = %e      hcur             = %e\n", 
           getStepSizeForLastInternalStep(),
           getStepSizeForNextInternalStep());
   os << buf;
   sprintf(buf, "tcur            = %e      tolsf            = %e\n", 
           getCurrentInternalValueOfIndependentVariable(),
           getCPODESSuggestedToleranceScalingFactor());
   os << buf;
}

int CPODESSolver::getNumberOfInternalStepsTaken() const
{
   long int r;
   int ierr = CPodeGetNumSteps(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfRHSFunctionCalls() const
{
   long int r;
   int ierr = CPodeGetNumFctEvals(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfLinearSolverSetupCalls() const
{
   long int r;
   int ierr = CPodeGetNumLinSolvSetups(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfNewtonIterations() const
{
   long int r;
   int ierr = CPodeGetNumNonlinSolvIters(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfNonlinearConvergenceFailures() const
{
   long int r;
   int ierr = CPodeGetNumNonlinSolvConvFails(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfLocalErrorTestFailures() const
{
   long int r;
   int ierr = CPodeGetNumErrTestFails(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getOrderUsedDuringLastInternalStep() const
{
   int r;
   int ierr = CPodeGetLastOrder(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getOrderToBeUsedDuringNextInternalStep() const
{
   int r;
   int ierr = CPodeGetCurrentOrder(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getCPODESMemoryUsageForDoubles() const
{
   long int r1;
   long int r2;
   int ierr = CPodeGetWorkSpace(d_cpode_mem, &r1, &r2);
   CPODE_SAMRAI_ERROR(ierr);
   return r1;
}

int CPODESSolver::getCPODESMemoryUsageForIntegers() const
{
   long int r1;
   long int r2;
   int ierr = CPodeGetWorkSpace(d_cpode_mem, &r1, &r2);
   CPODE_SAMRAI_ERROR(ierr);
   return r2;
}

double CPODESSolver::getStepSizeForLastInternalStep() const
{
   realtype r;
   int ierr = CPodeGetLastStep(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

double CPODESSolver::getStepSizeForNextInternalStep() const
{
   realtype r;
   int ierr = CPodeGetCurrentStep(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

double CPODESSolver::getCurrentInternalValueOfIndependentVariable() const
{
   realtype r;
   int ierr = CPodeGetCurrentStep(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

double CPODESSolver::getCPODESSuggestedToleranceScalingFactor() const
{
   realtype r;
   int ierr = CPodeGetTolScaleFactor(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

/*
*************************************************************************
*                                                                       *
* Access methods for CPSpgmr statistics.                                *
*                                                                       *
*************************************************************************
*/

void CPODESSolver::printCPSpgmrStatistics(std::ostream& os) const
{
   if (d_iteration_type == CP_NEWTON) {
      char buf[STAT_OUTPUT_BUFFER_SIZE];
      
      os << "CPODESSolver: CPSpgmr statistics... " << std::endl;
      
      sprintf(buf, "spgmr_lrw       = %5d     spgmr_liw        = %5d\n", 
              getCPSpgmrMemoryUsageForDoubles(),
              getCPSpgmrMemoryUsageForIntegers());
      os << buf;
      sprintf(buf, "nli             = %5d     ncfl             = %5d\n", 
              getNumberOfLinearIterations(),
              getNumberOfLinearConvergenceFailures());
      os << buf;
      sprintf(buf, "npe             = %5d     nps              = %5d\n", 
              getNumberOfPreconditionerEvaluations(),
              getNumberOfPrecondSolveCalls());
      os << buf;
   } else {
      
      os << "\nCPODESSolver not set to use NEWTON iteration . . . \n"
         << std::endl;
      
   }
}
   
int CPODESSolver::getNumberOfPreconditionerEvaluations() const
{
   long int r;
   int ierr = CPSpilsGetNumPrecEvals(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfLinearIterations() const
{
   long int r;
   int ierr = CPSpilsGetNumLinIters(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfPrecondSolveCalls() const
{
   long int r;
   int ierr = CPSpilsGetNumPrecSolves(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getNumberOfLinearConvergenceFailures() const
{
   long int r;
   int ierr = CPSpilsGetNumConvFails(d_cpode_mem, &r);
   CPODE_SAMRAI_ERROR(ierr);
   return r;
}

int CPODESSolver::getCPSpgmrMemoryUsageForDoubles() const
{
   long int r1;
   long int r2;
   int ierr = CPodeGetWorkSpace(d_cpode_mem, &r1, &r2);
   CPODE_SAMRAI_ERROR(ierr);
   return r1;
}

int CPODESSolver::getCPSpgmrMemoryUsageForIntegers() const
{
   long int r1;
   long int r2;
   int ierr = CPodeGetWorkSpace(d_cpode_mem, &r1, &r2);
   CPODE_SAMRAI_ERROR(ierr);
   return r2;
}

/*
*************************************************************************
*                                                                       *
* Print CPODESSolver object data to given output stream.            *
*                                                                       *
*************************************************************************
*/
void CPODESSolver::printClassData(std::ostream& os) const
{
   os << "\nCPODESSolver object data members..." << std::endl;
   os << "Object name = "
      << d_object_name << std::endl;

   os << "this = " << (CPODESSolver*)this << std::endl;
   os << "d_solution_vector = " 
      << (solv::SundialsAbstractVector*)d_solution_vector << std::endl;

   os << "d_CPODE_functions = " 
      << (CPODESAbstractFunctions*)d_cpode_functions << std::endl;

   os << "&d_cpode_mem = " << d_cpode_mem << std::endl;
   os << "d_cpode_log_file = " << (FILE*) d_cpode_log_file << std::endl;
   os << "d_cpode_log_file_name = " << d_cpode_log_file_name << std::endl;

   os << std::endl;
   os << "CPODES parameters..." << std::endl;
   os << "d_t_0 = "
      << d_t_0 << std::endl;
   os << "d_ic_vector = "
      << (solv::SundialsAbstractVector*)d_ic_vector << std::endl;

   os << "d_linear_multistep_method = "
      << d_linear_multistep_method << std::endl;
   os << "d_iteration_type = "
      << d_iteration_type << std::endl;
   os << "d_tolerance_type = "
      << d_tolerance_type << std::endl;
   os << "d_relative_tolerance = "
      << d_relative_tolerance << std::endl;
   os << "d_use_scalar_absolute_tolerance = ";
   if (d_use_scalar_absolute_tolerance) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   } 
   os << "d_absolute_tolerance_scalar = "
      << d_absolute_tolerance_scalar << std::endl;
   os << "d_absolute_tolerance_vector= " << std::endl;
   d_absolute_tolerance_vector->printVector(); 

   os << "Optional CPODES inputs (see CPODES docs for details):" 
      << std::endl;

   os << "maximum linear multistep method order = "
      << d_max_order << std::endl;
   os << "maximum number of internal steps = "
      << d_max_num_internal_steps << std::endl;
   os << "maximum number of nil internal step warnings = "
      << d_max_num_warnings << std::endl;

   os << "initial step size = "
      << d_init_step_size << std::endl;
   os << "maximum absolute value of step size = "
      << d_max_step_size << std::endl;
   os << "minimum absolute value of step size = "
      << d_min_step_size << std::endl;
   os << "last step size = "
      << getStepSizeForLastInternalStep() << std::endl;
   os << "...end of CPODES parameters\n" << std::endl;

   os << std::endl;
   os << "d_nonlin_conv_coef = "
        << d_nonlin_conv_coef << std::endl;

   os << std::endl;
   os << "CPSpgmr parameters..." << std::endl;
   os << "d_precondition_type = "
        << d_precondition_type << std::endl;
   os << "d_gram_schmidt_type = "
        << d_gram_schmidt_type << std::endl;
   os << "d_max_krylov_dim = "
        << d_max_krylov_dim << std::endl;
   os << "d_tol_scale_factor = "
        << d_tol_scale_factor << std::endl;
   os << "...end of CPSpgmr parameters\n" << std::endl;

   os << "d_CPODE_needs_initialization = ";
   if (d_CPODE_needs_initialization) {
      os << "true" << std::endl;
   } else {
      os << "false" << std::endl;
   } 

   os << "...end of CPODESSolver object data members\n" << std::endl;
}

std::vector< solv::SundialsAbstractVector* >*
CPODESSolver::getVectorsRequiringRegrid( void )const
{
   assert( d_cpode_mem!=NULL );

   std::vector< solv::SundialsAbstractVector* >* sundials_vec =
      new std::vector< solv::SundialsAbstractVector* >;

   CPodeMem mem = (CPodeMem) d_cpode_mem;

//   addVectorToList( sundials_vec, mem->cp_Vabstol );
   addVectorToList( sundials_vec, mem->cp_ewt );
   addVectorToList( sundials_vec, mem->cp_y );
   addVectorToList( sundials_vec, mem->cp_yp );
   addVectorToList( sundials_vec, mem->cp_acor );
   addVectorToList( sundials_vec, mem->cp_tempv );
   addVectorToList( sundials_vec, mem->cp_ftemp );
   addVectorToList( sundials_vec, mem->cp_acorP );
   addVectorToList( sundials_vec, mem->cp_errP );
//   addVectorToList( sundials_vec, mem->cp_yC );
//   addVectorToList( sundials_vec, mem->cp_ctol );
//   addVectorToList( sundials_vec, mem->cp_ctemp );
//   addVectorToList( sundials_vec, mem->cp_tempvP1 );
//   addVectorToList( sundials_vec, mem->cp_tempvP2 );
//   addVectorToList( sundials_vec, mem->cp_VabstolQ );
//   addVectorToList( sundials_vec, mem->cp_ewtQ );
//   addVectorToList( sundials_vec, mem->cp_yQ );
//   addVectorToList( sundials_vec, mem->cp_acorQ );
//   addVectorToList( sundials_vec, mem->cp_tempvQ );
//   addVectorToList( sundials_vec, mem->cp_yy0 );
//   addVectorToList( sundials_vec, mem->cp_yp0 );

   // The arrays are size qmax+1, so the <= is intentional
   for ( int i = 0; i <= mem->cp_qmax_alloc; i++ ) {
      addVectorToList( sundials_vec, mem->cp_zn[i] );
//      addVectorToList( sundials_vec, mem->cp_znQ[i] );
   }
   
   CPSpilsMem cpspils_mem = (CPSpilsMem) mem->cp_lmem;

   addVectorToList( sundials_vec, cpspils_mem->s_ytemp );
   addVectorToList( sundials_vec, cpspils_mem->s_yptemp );
   addVectorToList( sundials_vec, cpspils_mem->s_x );
   addVectorToList( sundials_vec, cpspils_mem->s_ycur );
   addVectorToList( sundials_vec, cpspils_mem->s_ypcur );
   addVectorToList( sundials_vec, cpspils_mem->s_fcur );

   SpgmrMem spgmr_mem = (SpgmrMem) cpspils_mem->s_spils_mem;

   // The V array is size l_max+1, so the <= is intentional
   for ( int i = 0; i <= spgmr_mem->l_max; i++ ) {
      addVectorToList( sundials_vec, spgmr_mem->V[i] );
   }

   addVectorToList( sundials_vec, spgmr_mem->xcor );
   addVectorToList( sundials_vec, spgmr_mem->vtemp );

   return sundials_vec;
}

void CPODESSolver::addVectorToList(
   std::vector< solv::SundialsAbstractVector* >* sundials_vec,
   N_Vector& n ) const
{
   if ( n != NULL ) {

      solv::SundialsAbstractVector* y = SABSVEC_CAST( n );
      sundials_vec->push_back( y );

   }
}

#endif
