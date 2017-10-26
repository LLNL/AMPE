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
#include "EllipticFACSolver.h"
#include "EllipticFACOps.h"

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>

using namespace std;

/*
*************************************************************************
*                                                                       *
* Constructor sets uninitialized solver state.                          *
* Set default iteration and convergence parameters.                     *
*                                                                       *
* By default settings:                                                  *
*   - Poisson equation specified has D=1, C=0.                          *
*   - State is uninitialized                                            *
*   - Logging is disabled                                               *
*   - Context for internal data is set based on object name.            *
*                                                                       *
*************************************************************************
*/


EllipticFACSolver::EllipticFACSolver (
   const std::string &object_name,
   boost::shared_ptr<EllipticFACOps> fac_ops,
   const boost::shared_ptr<tbox::Database>& database )
   :
   d_object_name(object_name),
   d_fac_ops(fac_ops),
   d_fac_precond(object_name+"::fac_precond",d_fac_ops),
   d_bc_object(NULL),
   d_simple_bc(tbox::Dimension(NDIM),object_name+"::bc"),
   d_ln_min(-1),
   d_ln_max(-1),
   d_context(hier::VariableDatabase::getDatabase()
             ->getContext(object_name+"::CONTEXT")) ,
   d_solver_is_initialized(false),
   d_enable_logging(false),
   d_verbose(false)
{
   boost::shared_ptr<pdat::CellVariable<double> > vol_var(
      new pdat::CellVariable<double>(tbox::Dimension(NDIM),object_name+"::weight"));
   d_vol_id= hier::VariableDatabase::getDatabase()->
      registerVariableAndContext(vol_var, d_context, hier::IntVector(tbox::Dimension(NDIM),0) );

   setMaxCycles(10);
   setResidualTolerance(1e-6);
   setCoarseFineDiscretization("Ewing");
#ifdef HAVE_HYPRE
   setCoarsestLevelSolverChoice("hypre");
   setCoarsestLevelSolverTolerance(1e-2);
#else
   setCoarsestLevelSolverChoice("redblack");
   setCoarsestLevelSolverTolerance(1e-8);
   setCoarsestLevelSolverMaxIterations(500);
#endif

   /*
    * The default RobinBcCoefStrategy used,
    * SimpleCellRobinBcCoefs only works with constant refine
    * for prolongation.  So we use constant refinement
    * for prolongation by default.
    */
   setProlongationMethod("CONSTANT_REFINE");

   /*
    * The FAC operator optionally uses the preconditioner
    * to get data for logging.
    */
   d_fac_ops->setPreconditioner( 
      (const FACPreconditioner*)(&d_fac_precond) );

   if ( database ) {
      getFromInput(database);
   }

   return;
}

/*
*************************************************************************
*                                                                       *
* Destructor for EllipticFACSolver.                                *
* Deallocate internal data.                                             *
*                                                                       *
*************************************************************************
*/


EllipticFACSolver::~EllipticFACSolver()
{
   deallocateSolverState();
   return;
}


/*
********************************************************************
* Set state from database                                          *
*                                                                  *
* Do not allow FAC preconditioner and Poisson FAC operators to be  *
* set from database, as that may cause them to be inconsistent     *
* with this object if user does not coordinate the inputs          *
* correctly.  This is also why we don't allow direct access to     *
* those objects.  The responsibility for maintaining consistency   *
* lies in the public functions to set parameters, so use them      *
* instead of setting the parameters directly in this function.     *
********************************************************************
*/

void
EllipticFACSolver::getFromInput(const boost::shared_ptr<tbox::Database>& database )
{
   if ( database->isBool("enable_logging") ) {
      bool logging = database->getBool("enable_logging");
      enableLogging(logging);
   }
   if ( database->isBool("verbose") ) {
      bool verbose = database->getBool("verbose");
      setVerbose(verbose);
   }
   if ( database->isInteger("max_cycles") ) {
      int max_cycles = database->getInteger("max_cycles");
      setMaxCycles(max_cycles);
   }
   
   double residual_tol = 1.e-6;
   if (database->isDouble("residual_tol")) {
      residual_tol = database->getDouble("residual_tol");
      setResidualTolerance(residual_tol, -1.);
   }
   if (database->isDouble("relative_residual_tol")) {
      double relative_residual_tol = database->getDouble(
            "relative_residual_tol");
      setResidualTolerance(residual_tol, relative_residual_tol);
   }
   if ( database->isString("coarse_fine_discretization") ) {
      std::string s = database->getString("coarse_fine_discretization");
      setCoarseFineDiscretization(s);
   }
   if ( database->isString("prolongation_method") ) {
      std::string s = database->getString("prolongation_method");
      setProlongationMethod(s);
   }
   if ( database->isString("coarse_solver_choice") ) {
      std::string s = database->getString("coarse_solver_choice");
      setCoarsestLevelSolverChoice(s);
   }
   if ( database->isDouble("coarse_solver_tolerance") ) {
      double tol = database->getDouble("coarse_solver_tolerance");
      setCoarsestLevelSolverTolerance(tol);
   }
   if ( database->isInteger("coarse_solver_max_iterations") ) {
      int itr = database->getInteger("coarse_solver_max_iterations");
      setCoarsestLevelSolverMaxIterations(itr);
   }

   return;
}


/*
*************************************************************************
*                                                                       *
* Prepare internal data for solve.                                      *
* Allocate scratch data.  Create vectors for u and f                    *
* required by the FACPreconditioner interface.                    *
* Set up internal boundary condition object.                            *
* Share data to coordinate with FAC preconditioner and                  *
* Poisson FAC operator.                                                 *
*                                                                       *
*************************************************************************
*/
void
EllipticFACSolver::initializeSolverState(
   const int solution ,
   const int rhs ,
   const boost::shared_ptr< hier::PatchHierarchy >& hierarchy,
   const int coarse_level,
   const int fine_level )
{
   TBOX_ASSERT(hierarchy);
   
   if ( d_bc_object == NULL ) {
      TBOX_ERROR(d_object_name << ": No BC coefficient strategy object!\n"
                 << "Use either setBoundaries or setPhysicalBcCoefObject\n"
                 << "to specify the boundary conidition.\n");
   }

   if ( ! d_solver_is_initialized ) {
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( solution < 0 || rhs < 0 ) {
         TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
      }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
      if ( !hierarchy ) {
         TBOX_ERROR(d_object_name << ": NULL hierarchy pointer not allowed\n"
                    << "in inititialization.");
      }
#endif
      d_hierarchy = hierarchy;

      d_ln_min = coarse_level;
      d_ln_max = fine_level;
      if ( d_ln_min == -1 ) {
         d_ln_min = 0;
      }
      if ( d_ln_max == -1 ) {
         d_ln_max = d_hierarchy->getFinestLevelNumber();
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      if ( d_ln_min < 0 || d_ln_max < 0 || d_ln_min > d_ln_max ) {
         TBOX_ERROR(d_object_name << ": Bad range of levels in\n"
                    << "inititialization.\n");
      }
#endif

      for (int ln=d_ln_min; ln<=d_ln_max; ++ln ) {
         d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_vol_id);
      }

      d_fac_ops->computeVectorWeights(
         d_hierarchy,
         d_vol_id,
         d_ln_min,
         d_ln_max );

      createVectorWrappers(solution,rhs);

      d_fac_precond.initializeSolverState(*d_uv, *d_fv);
      //d_hopscell.reset(new math::HierarchyCellDataOpsReal<double>(hierarchy,
      //                                                         d_ln_min,
      //                                                         d_ln_max));

      d_solver_is_initialized = true;
   }

   return;
}



void
EllipticFACSolver::finalizeCoefficients()
{

  if ( d_bc_object == &d_simple_bc ) {
    d_simple_bc.setHierarchy(d_hierarchy,
                             d_ln_min,
                             d_ln_max);
    if ( d_fac_ops->dIsConstant() ) {
      d_simple_bc.setDiffusionCoefConstant( d_fac_ops->getDConstant() );
    }
    else {
      d_simple_bc.setDiffusionCoefId( d_fac_ops->getDPatchDataId() );
    }
  }

   d_fac_ops->finalizeCoefficients();

   return;
}

void
EllipticFACSolver::deallocateSolverState()
{
   if ( d_hierarchy ) {

      d_fac_precond.deallocateSolverState();

      /*
       * Delete internally managed data.
       */
      int ln;
      for ( ln=d_ln_min; ln<=d_ln_max; ++ln ) {
         d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_vol_id);
      }

      d_hierarchy.reset();
      d_ln_min = -1;
      d_ln_max = -1;
      d_solver_is_initialized = false;

      destroyVectorWrappers();
   }
   return;
}



void
EllipticFACSolver::resetSolverState(
   const int soln_id,
   const int  rhs_id,
   const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   if( d_solver_is_initialized ){
      assert( hierarchy );
      assert( soln_id>=0 );
      assert( rhs_id>=0 );
      deallocateSolverState();
      initializeSolverState(soln_id, rhs_id, hierarchy);
   }
}


void
EllipticFACSolver::setBoundaries(const std::string& boundary_type,
                                 const int                 fluxes,
                                 const int                  flags,
                                 int *                 bdry_types)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_bc_object != NULL && d_bc_object != &d_simple_bc ) {
      TBOX_ERROR(d_object_name << ": Bad attempt to set boundary condition\n"
                 << "by using default bc object after it has been overriden.\n");
   }
#endif
   d_simple_bc.setBoundaries(boundary_type,
                             fluxes,
                             flags,
                             bdry_types);
   d_bc_object = &d_simple_bc;
   d_fac_ops->setPhysicalBcCoefObject(d_bc_object);
}


/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an initialized solver state.                      *
* Before solving, set the final piece of the boundary condition,        *
* which is not known until now, and initialize some internal            *
* solver quantities.                                                    *
*                                                                       *
*************************************************************************
*/
bool
EllipticFACSolver::solveSystem(const int u_id,
                               const int f_id,
                               const int ew_id)
{
   //tbox::pout<<"EllipticFACSolver::solveSystem() for object "<<d_object_name<<endl;
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( u_id!=-1 );
   assert( f_id!=-1 );
   if ( !d_solver_is_initialized ) {
      TBOX_ERROR(d_object_name << ".solveSystem(int,int): uninitialized\n"
                 << "solver state.  You must call initializeSolverState()\n"
                 << "before using this function.  Or you can use\n"
                 << "solveSystem(int,int,...) to initialize the solver,\n"
                 << "solve and deallocate the solver.\n");
   }
   if ( u_id < 0 || f_id < 0 ) {
      TBOX_ERROR(d_object_name << ": Bad patch data id.\n");
   }
#endif
   if ( d_bc_object == &d_simple_bc ) {
      /*
       * Knowing that we are using the SimpelCellRobinBcCoefsX
       * implementation of RobinBcCoefStrategy, we must save
       * the ghost data in u before solving.
       * The solver overwrites it, but SimpleCellRobinBcCoefs
       * needs to get to access it repeatedly.
       */
      d_simple_bc.cacheDirichletData(u_id);
   }

   createVectorWrappers(u_id,f_id);

   d_fac_ops->setWeightIds(ew_id, d_vol_id);

#if 0
   d_fac_precond.printClassData(tbox::pout);

   math::HierarchyCellDataOpsReal<double> hopscell(d_hierarchy);
   double normf = hopscell.weightedRMSNorm(f_id, ew_id, d_vol_id);
   tbox::pout << "EllipticFACSolver ("<<d_object_name<<"), f: Weighted RMS norm on composite grid = " << normf << endl;
   double normu = hopscell.weightedRMSNorm(u_id, ew_id, d_vol_id);
   tbox::pout << "EllipticFACSolver ("<<d_object_name<<"), u: Weighted RMS norm on composite grid = " << normu << endl;
#endif

   bool solver_rval = d_fac_precond.solveSystem( *d_uv, *d_fv );

   if ( d_bc_object == &d_simple_bc ) {
      /*
       * Restore the Dirichlet cell data that were overwritten by the
       * solve process.  We do this to be backward compatible with the
       * user code.
       */
      d_simple_bc.restoreDirichletData(u_id);
   }

   if (d_verbose) printFACConvergenceFactors(solver_rval);

   d_fac_ops->setWeightIds(-1,-1);
   
   return solver_rval;
}



/*
*************************************************************************
*                                                                       *
* Solve the linear system and report whether iteration converged.       *
*                                                                       *
* This version is for an uninitialized solver state.                    *
* 1. Initialize the (currently uninitialized) solver state.             *
* 2. Solve.                                                             *
* 3. Deallocate the solver state.                                       *
*                                                                       *
*************************************************************************
*/

 bool
EllipticFACSolver::solveSystem(const int u_id,
                               const int f_id,
                               const int w_id,
                               const boost::shared_ptr< hier::PatchHierarchy >& hierarchy,
                               int coarse_ln,
                               int fine_ln)
{
   TBOX_ASSERT(hierarchy);

   if ( d_enable_logging ) {
      tbox::plog << "EllipticFACSolver::solveSystem (" << d_object_name
           << ")\n";
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   if ( d_solver_is_initialized ) {
      TBOX_ERROR(d_object_name << ".solveSystem(int,int,...): initialized\n"
                 << "solver state.  This function can only used when the\n"
                 << "solver state is uninitialized.  You should deallocate\n"
                 << "the solver state or use solveSystem(int,int).\n");
   }
   if ( !hierarchy ) {
      TBOX_ERROR(d_object_name << ".solveSystem(): Null hierarchy\n"
                 << "specified.\n");
   }
#endif
   initializeSolverState( u_id, f_id, hierarchy, coarse_ln, fine_ln );

   bool solver_rval= solveSystem( u_id, f_id, w_id);

   deallocateSolverState();

   return solver_rval;
}




void
EllipticFACSolver::createVectorWrappers(int u,
                                        int f)
{

   hier::VariableDatabase &vdb(*hier::VariableDatabase::getDatabase());
   boost::shared_ptr< hier::Variable > variable;

   if ( !d_uv || d_uv->getComponentDescriptorIndex(0) != u ) {
      d_uv.reset( new solv::SAMRAIVectorReal<double>(d_object_name+"::uv",
                                                d_hierarchy,
                                                d_ln_min,
                                                d_ln_max) );
      vdb.mapIndexToVariable(u, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( !variable ) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                    << u << "\n");
      }
      boost::shared_ptr<pdat::CellVariable<double> > cell_variable ( 
         BOOST_CAST<pdat::CellVariable<double>,hier::Variable>(variable) );
      if ( !cell_variable ) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << u
                    << " is not a cell-double variable.\n");
      }
#endif
      d_uv->addComponent( variable, u, d_vol_id );
   }

   if ( !d_fv || d_fv->getComponentDescriptorIndex(0) != f ) {
      d_fv.reset( new solv::SAMRAIVectorReal<double>(d_object_name+"::fv",
                                                d_hierarchy,
                                                d_ln_min,
                                                d_ln_max));
      vdb.mapIndexToVariable(f, variable);
#ifdef DEBUG_CHECK_ASSERTIONS
      if ( !variable ) {
         TBOX_ERROR(d_object_name << ": No variable for patch data index "
                    << f << "\n");
      }
      boost::shared_ptr<pdat::CellVariable<double> > cell_variable ( 
         BOOST_CAST<pdat::CellVariable<double>,hier::Variable>(variable) );
      if ( !cell_variable ) {
         TBOX_ERROR(d_object_name << ": hier::Patch data index " << f
                    << " is not a cell-double variable.\n");
      }
#endif
      d_fv->addComponent( variable, f, d_vol_id );
   }

   return;
}


void
EllipticFACSolver::printFACConvergenceFactors(const int solver_ret)
{
  d_fac_precond.printClassData(tbox::pout);
  double avg_factor, final_factor;
  getConvergenceFactors(avg_factor, final_factor);
  tbox::pout << "  EllipticFACSolver "<<d_object_name<<", iteration ";
  tbox::pout << (solver_ret?"":"NOT ") << "converged " << "\n"
             << "     iterations: " << getNumberOfIterations() << "\n"
             << "     residual: " << getResidualNorm() << "\n"
             << "     average convergence: " << avg_factor << "\n"
             << "     final convergence: " << final_factor << "\n"
             << std::flush;
}
