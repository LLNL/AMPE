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
 * This class is a high-level wrapper for a quaternion system
 * solver.  See the header file QuatSysSolver.h for additional
 * documentation of the class member functions and data.
 *
 * This file was adapted from solv::CellPoissonFACSolver.C
 * in the SAMRAI library.
 */

#include "QuatSysSolver.h"
#include <cassert>

#include "SAMRAI/tbox/TimerManager.h"
using namespace std;

/*
*************************************************************************
*                                                                       *
* Constructor sets uninitialized solver state.                          *
* Set default iteration and convergence parameters.                     *
*                                                                       *
* By default settings:                                                  *
*   - State is uninitialized                                            *
*   - Logging is disabled                                               *
*   - Context for internal data is set based on object name.            *
*                                                                       *
*************************************************************************
*/


QuatSysSolver::QuatSysSolver(
   const int ql,
   const string &             object_name,
   boost::shared_ptr<tbox::Database> database )
  : d_object_name(object_name),
    d_fac_ops(new QuatFACOps(ql, object_name+"::fac_ops",database)),
    d_fac_solver(object_name+"::fac_precond",d_fac_ops),
    d_bc_object(NULL),
    d_simple_bc(tbox::Dimension(NDIM),object_name+"::bc"),
    d_ln_min(-1),
    d_ln_max(-1),
    d_context(hier::VariableDatabase::getDatabase()
	      ->getContext(object_name+"::CONTEXT")) ,
    d_solver_is_initialized(false),
    d_verbose(false),
    d_precond_maxiters(10)
{
  // Get timers
  t_set_op_coef = tbox::TimerManager::getManager()->getTimer("QuatSysSolver::setOperatorCoefficients");
  t_solve_system = tbox::TimerManager::getManager()->getTimer("QuatSysSolver::solveSystem");

  /*
    Set some default parameters.  These may be overridden by the respective
    values in the input database (if supplied).
  */
  setLevelSolverTolerance(1.e-2);
  setLevelSolverMaxIterations(20);
  setCoarsestLevelSolverTolerance(1.e-2);
  setCoarsestLevelSolverMaxIterations(20);

  /*
   * The default RobinBcCoefStrategy used,
   * SimpleCellRobinBcCoefs only works with constant refine
   * for prolongation.  So we use constant refinement
   * for prolongation by default.
   */
  setProlongationMethod("CONSTANT_REFINE");

  /*
   * The FAC operator optionally uses (via callback through the pointer
   * being supplied here) the FAC solver object to get data for logging.
   */
  d_fac_ops->setSolver((const FACPreconditioner*)(&d_fac_solver));

  if ( database ) {
     getFromInput(database);
  }
}

/*************************************************************************
*                                                                       *
* Deallocate internal data.                                             *
*                                                                       *
*************************************************************************/
QuatSysSolver::~QuatSysSolver()
{
  deallocateSolverState();
}


/*
********************************************************************
* Set state from database                                          *
*                                                                  *
* In constructing the FAC operator and solver objects owned by     *
* this class, the input database (assuming one was provided) is    *
* deliberately NOT passed into the respective constructors to      *
* preclude these objects from directly accessing the database      *
* themselves.  Instead, all input parameters controlling the FAC   *
* iteration are read by this wrapper object and passed into the    *
* FAC operator and solver objects through their public accessor    *
* functions.  This is done so that all of the FAC solver input     *
* parameters are set in one place (rather than being spread over   *
* a few databases and sub-databases) and to reduce the opportunity *
* for some inconsistent or unintended parameter specification.     *
********************************************************************
*/
void
QuatSysSolver::getFromInput(const boost::shared_ptr<tbox::Database>& input_db )
{
   if ( input_db->isBool("verbose") ) {
     bool verbose = input_db->getBool("verbose");
     /*
       Controls printing of solver convergence information
       to the screen
     */
     setVerbose(verbose);
   }
   if ( input_db->isString("coarse_fine_discretization") ) {
     string s = input_db->getString("coarse_fine_discretization");
     // Sets the coarse-fine discretization to be used by the FAC operator
     setCoarseFineDiscretization(s);
   }
   if ( input_db->isDouble("levelsolver_tolerance") ) {
     double tol = input_db->getDouble("levelsolver_tolerance");
     // Sets the level solver tolerance to be used by the FAC operator
     setLevelSolverTolerance(tol);
     if ( input_db->isDouble("coarse_levelsolver_tolerance") ) {
       tol = input_db->getDouble("coarse_levelsolver_tolerance");
     }
     // Sets the coarsest level solver tolerance to be used by the FAC operator
     setCoarsestLevelSolverTolerance(tol);
   }
   if ( input_db->isInteger("levelsolver_max_iterations") ) {
     int itr = input_db->getInteger("levelsolver_max_iterations");
     // Sets the level solver max iteration to be used by the FAC operator
     setLevelSolverMaxIterations(itr);
     if ( input_db->isInteger("coarse_levelsolver_max_iterations") ) {
       itr = input_db->getInteger("coarse_levelsolver_max_iterations");
     }
     // Sets the coarsest level solver max iteration to be used by the FAC operator
     setCoarsestLevelSolverMaxIterations(itr);
   }
   if ( input_db->isString("prolongation_method") ) {
     string s = input_db->getString("prolongation_method");
     // Sets the prolongation method to be used by the FAC operator
     setProlongationMethod(s);
   }
}


/*************************************************************************
*                                                                       *
* Prepare internal data for solve.                                      *
* Allocate scratch data.  Create vectors for u and f                    *
* required by the FACPreconditioner interface.                     *
* Set up internal boundary condition object.                            *
* Share data to coordinate with FAC preconditioner and                  *
* QuatFAC operator.                                                     *
*                                                                       *
*************************************************************************/
void
QuatSysSolver::initializeSolverState(const int q_soln_id,
   const int q_rhs_id,
   const int weight_id,
   boost::shared_ptr< hier::PatchHierarchy > hierarchy)
{
   TBOX_ASSERT(hierarchy);
   assert( q_soln_id>=0 );

   if (d_bc_object == NULL) {
      TBOX_ERROR(
         d_object_name << ": No BC coefficient strategy object!\n"
                       << "Use either setBoundaries or setPhysicalBcCoefObject\n"
                       << "to specify the boundary conidition.\n");
   }

   if(!d_solver_is_initialized){
#ifdef DEBUG_CHECK_ASSERTIONS
     if ( q_soln_id < 0 || q_rhs_id < 0 || weight_id < 0 ) {
       TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
     }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
     if (!hierarchy) {
        TBOX_ERROR(d_object_name << ": NULL hierarchy pointer not allowed\n"
                               << "in inititialization.");
     }
#endif
     d_hierarchy = hierarchy;

     d_ln_min = 0;
     d_ln_max = hierarchy->getFinestLevelNumber();
     if ( d_ln_min == -1 ) {
       d_ln_min = 0;
     }
     if ( d_ln_max == -1 ) {
       d_ln_min = d_hierarchy->getFinestLevelNumber();
     }

#ifdef DEBUG_CHECK_ASSERTIONS
     if ( d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max ) {
       TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
		  << "inititialization.\n");
     }
#endif

     // Set the boundary condition object

     if ( d_bc_object == &d_simple_bc ) {
       d_simple_bc.setHierarchy(d_hierarchy, d_ln_min, d_ln_max);
     }

     d_weight_id = weight_id;

     // Wrap the temporary vectors
     createVectorWrappers(q_soln_id, q_rhs_id, d_uv, d_fv);

     // Initialize the FAC solver used for the preconditioner
     d_fac_solver.initializeSolverState(*d_uv, *d_fv);

     d_solver_is_initialized = true;
   }
}



void
QuatSysSolver::resetSolverState(
   const int q_soln_id,
   const int  q_rhs_id,
   const int weight_id,
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy)
{
   assert( q_soln_id!=-1 );
   assert( weight_id!=-1 );
   
   if( d_solver_is_initialized ){
      deallocateSolverState();
      initializeSolverState(q_soln_id, q_rhs_id, weight_id,
                            hierarchy);
   }
}



void
QuatSysSolver::deallocateSolverState()
{
  if ( d_hierarchy ) {

    // Deallocate FAC preconditioner
    d_fac_solver.deallocateSolverState();

    /*
     * Delete internally managed data.
     */
    //    int ln;
    //    for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
    //      d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_weight_id);
    //    }

    d_hierarchy.reset();
    d_ln_min = -1;
    d_ln_max = -1;
    d_solver_is_initialized = false;

    destroyVectorWrappers(d_uv, d_fv);
  }
}



void
QuatSysSolver::setBoundaries(const string& boundary_type,
            const            int fluxes,
            const             int flags,
            int *            bdry_types)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  if ( d_bc_object != NULL && d_bc_object != &d_simple_bc ) {
    TBOX_ERROR(d_object_name << ": Bad attempt to set boundary condition\n"
        << "by using default bc object after it has been overriden.\n");
  }
#endif
  d_simple_bc.setBoundaries(boundary_type, fluxes, flags, bdry_types);
  d_bc_object = &d_simple_bc;
  d_fac_ops->setPhysicalBcCoefObject(d_bc_object);
}


void
QuatSysSolver::setOperatorCoefficients(
   const double        time_step,
   const double        epsilon_q,
   const double  quat_grad_floor,
   const string quat_smooth_floor_type,
   const int         mobility_id,
   const int   mobility_deriv_id,
   const int        diff_coef_id,
   const int  diff_coef_deriv_id,
   const int           grad_q_id,
   const int                q_id )
{
  t_set_op_coef->start();

#ifdef DEBUG_CHECK_ASSERTIONS
  if ( !d_solver_is_initialized ) {
    TBOX_ERROR(d_object_name << ".setOperatorCoefficients(): uninitialized\n"
	       << "solver state.  You must call initializeSolverState()\n"
	       << "before using this function.\n");
  }
  if ( mobility_id < 0 || diff_coef_id < 0 ||
       grad_q_id < 0 || q_id < 0) {
    TBOX_ERROR(d_object_name << "::setOperatorCoefficients: Bad patch data id.\n");
  }
#endif

  d_fac_ops->setOperatorCoefficients(
     time_step, epsilon_q, mobility_id, mobility_deriv_id, diff_coef_id, diff_coef_deriv_id,
     grad_q_id, q_id, quat_grad_floor, quat_smooth_floor_type );

  t_set_op_coef->stop();
}


bool
QuatSysSolver::solveSystem(const int q_soln_id,
                           const int q_rhs_id,
                           const int q_ewt_id)
{
  t_solve_system->start();

  assert( d_weight_id!=-1 );
  assert( q_ewt_id!=-1 );

  boost::shared_ptr<solv::SAMRAIVectorReal<double> > solution;
  solution.reset();
  boost::shared_ptr<solv::SAMRAIVectorReal<double> > rhs;
  rhs.reset();

  createVectorWrappers(q_soln_id, q_rhs_id, solution, rhs);

  d_fac_ops->setWeightIds(q_ewt_id, d_weight_id);

  // Divide by the square root of the mobility
  d_fac_ops->divideMobilitySqrt(q_rhs_id);

  //  setResidualTolerance(d_precond_tol);

  solution->setToScalar(0., false);

  if ( d_bc_object == &d_simple_bc ) {
    /*
     * Knowing that we are using the SimpleCellRobinBcCoefsX
     * implementation of RobinBcCoefStrategy, we must save
     * the ghost data in u before solving.
     * The solver overwrites it, but SimpleCellRobinBcCoefs
     * needs to get to access it repeatedly.
     */
    d_simple_bc.cacheDirichletData(q_soln_id);
  }

  bool solver_rval = d_fac_solver.solveSystem( *solution, *rhs );

  // Multiply by the square root of the mobility
  d_fac_ops->multiplyMobilitySqrt(q_soln_id);

  if ( d_bc_object == &d_simple_bc ) {
    /*
     * Restore the Dirichlet cell data that were overwritten by the
     * solve process.  We do this to be backward compatible with the
     * user code.
     */
    d_simple_bc.restoreDirichletData(q_soln_id);
  }

  if (d_verbose) {
    printFACConvergenceFactors(solver_rval);
  }

  d_fac_ops->setWeightIds(-1, -1);

  destroyVectorWrappers(solution, rhs);

  t_solve_system->stop();

  return solver_rval;
}



void
QuatSysSolver::evaluateRHS(
   const double      epsilon_q,
   const double quat_grad_floor,
   const string quat_floor_type,
   const int  diffusion_coef_id, // D_q(phi)
   const int         grad_q_id,
   const int    grad_q_copy_id,// for computation of diffusion coefficient
   const int        rotations_id,
   const int       mobility_id,
   const int       solution_id,
   int                  rhs_id,
   const bool use_gradq_for_flux )
{
   assert( diffusion_coef_id>=0 );
   assert( grad_q_id>=0 );
   assert( mobility_id>=0 );
   assert( solution_id>=0 );
   assert( rhs_id>=0 );

   d_fac_ops->evaluateRHS(epsilon_q, diffusion_coef_id, grad_q_id, grad_q_copy_id, 
                          quat_grad_floor, quat_floor_type, mobility_id, 
                          rotations_id, solution_id, rhs_id,
                          use_gradq_for_flux);
}



void
QuatSysSolver::multiplyDQuatDPhiBlock(
   const int          q_id,
   const int operator_q_id)
{
   assert( q_id>=0 );
   assert( operator_q_id>=0 );

   d_fac_ops->multiplyDQuatDPhiBlock(q_id, operator_q_id);
}



void
QuatSysSolver::applyProjection(const int      q_id,
				    const int   corr_id,
				    const int    err_id)
{
   assert( q_id>=0 );
   assert( corr_id>=0 );
   assert( err_id>=0 );

   for (int ln=d_ln_max; ln>=d_ln_min; ln--) {
      d_fac_ops->applyProjectionOnLevel(q_id, corr_id, err_id, ln);
   }
}


void
QuatSysSolver::printFACConvergenceFactors(const int solver_ret)
{
  double avg_factor, final_factor;
  getFACConvergenceFactors(avg_factor, final_factor);
  tbox::pout << "Preconditioner FAC iteration ";
  tbox::pout << (solver_ret?"":"NOT ") << "converged " << "\n"
	     << "     iterations: " << getNumberOfIterations() << "\n"
	     << "     residual: " << getResidualNorm() << "\n"
	     << "     average convergence: " << avg_factor << "\n"
	     << "     final convergence: " << final_factor << "\n"
	     << flush;
}


void
QuatSysSolver::createVectorWrappers(int q_u,
			 int q_f,
			 boost::shared_ptr<solv::SAMRAIVectorReal<double> > & uv,
			 boost::shared_ptr<solv::SAMRAIVectorReal<double> > & fv)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert( d_weight_id!=-1 );
#endif

  hier::VariableDatabase &vdb(*hier::VariableDatabase::getDatabase());
  boost::shared_ptr< hier::Variable > variable;

  /*
   *  If the vector uv has not yet been allocated, or if either of
   *  the component ids have changed, then reconstruct the vector
   */
  if ( !uv
       || uv->getComponentDescriptorIndex(0) != q_u ) {
    uv.reset( new solv::SAMRAIVectorReal<double>(d_object_name+"::uv",
						  d_hierarchy,
						  d_ln_min,
						  d_ln_max) );
    vdb.mapIndexToVariable(q_u, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!variable) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                                  << q_u << "\n");
      }
      boost::shared_ptr<pdat::CellVariable<double> > cell_variable(
         BOOST_CAST<pdat::CellVariable<double>, hier::Variable>(variable));
      if (!cell_variable) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << q_u
                                  << " is not a cell-double variable.\n");
      }
#endif
    uv->addComponent( variable, q_u, d_weight_id );
  }


  /*
   *  If the vector d_fv has not yet been allocated, or if either of
   *  the component ids have changed, then reconstruct the vector
   */
  if ( !fv
       || fv->getComponentDescriptorIndex(0) != q_f ) {
    fv.reset( new solv::SAMRAIVectorReal<double>(d_object_name+"::fv",
						  d_hierarchy,
						  d_ln_min,
						  d_ln_max) );
    vdb.mapIndexToVariable(q_f, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if (!variable) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                                  << q_f << "\n");
      }
      boost::shared_ptr<pdat::CellVariable<double> > cell_variable(
         BOOST_CAST<pdat::CellVariable<double>, hier::Variable>(variable));
      if (!cell_variable) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << q_f
                                  << " is not a cell-double variable.\n");
      }
#endif
    fv->addComponent( variable, q_f, d_weight_id );

  }
}


/*
***********************************************************************
* Delete the vector wrappers.  Do not freeVectorComponents because    *
* we do not control their data allocation.  The user does that.       *
***********************************************************************
*/
void
QuatSysSolver::destroyVectorWrappers(
   boost::shared_ptr<solv::SAMRAIVectorReal<double> > & uv,
   boost::shared_ptr<solv::SAMRAIVectorReal<double> > & fv)
{
  uv.reset();
  fv.reset();
}
