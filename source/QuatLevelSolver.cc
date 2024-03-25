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
/*
 * This class implements a level solver for a quaternion system
 * FAC solver.  See the header file QuatLevelSolver.h for
 * additional documentation of the class member functions and
 * data.
 *
 * This file was adapted from solv::CellPoissonHypreSolver.C in
 * the SAMRAI library.
 */

#define SYMMETRIC_STENCIL

#include "QuatLevelSolver.h"
#include "krylov.h"
#include "fc_internal_mangle.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"

#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/math/ArrayDataNormOpsReal.h"


#include <cassert>


#define NPARTS 1
#define PART 0

extern "C" {

#if NDIM == 2

// 2d prototypes for the Fortran subroutines defined in
// quaternionlevelsolver2d.f

void SET_J_IJ2D(const int &, const int &, const int &, const int &,
                const double *, const int &, const int &, const int &,
                const int &, const double *, const int &, const int &,
                const int &, const int &, const double *, const double *,
                const int &, const int &, const int &, const int &,
                const double *, const int &, const int &, const int &,
                const int &, const double *, const int &, const int &,
                const int &, const int &);
void SET_STENCIL2D(const int &, const int &, const int &, const int &,
                   const double &, const double *, const int &, const int &,
                   const int &, const int &, const double *, const int &,
                   const int &, const int &, const int &, const double *,
                   const int &, const int &, const int &, const int &,
                   const double *, const int &, const int &, const int &,
                   const int &, double *);
void SET_SYMMETRIC_STENCIL2D(const int &, const int &, const int &, const int &,
                             const double &, const double *, const int &,
                             const int &, const int &, const int &,
                             const double *, const int &, const int &,
                             const int &, const int &, const double *,
                             const int &, const int &, const int &, const int &,
                             const double *, const int &, const int &,
                             const int &, const int &, double *);
void COPY2D(const int &, const int &, const int &, const int &, const double *,
            const int &, const int &, const int &, const int &, const double *,
            const int &, const int &, const int &, const int &);
void ADJUST_BDRY2D(const double &, const double *, const int &, const int &,
                   const int &, const int &, double *diag,
                   const double *offdiagi, const double *offdiagj,
                   const int *pifirst, const int *pilast, const int *pjfirst,
                   const int *pjlast, const double *acoef, const double *bcoef,
                   const int *aifirst, const int *ailast, const int *ajfirst,
                   const int *ajlast, const double *Ak0, const int *kifirst,
                   const int *kilast, const int *kjfirst, const int *kjlast,
                   const int *lower, const int *upper, const int *location,
                   const double *h);
void ADJUST_QRHS2D(double *rhs, const int *rifirst, const int *rilast,
                   const int *rjfirst, const int *rjlast, const double *Ak0,
                   const int *kifirst, const int *kilast, const int *kjfirst,
                   const int *kjlast, const double *gcoef, const int *aifirst,
                   const int *ailast, const int *ajfirst, const int *ajlast,
                   const int *lower, const int *upper, const int *location);

#endif
#if NDIM == 3

// 3d prototypes for the Fortran subroutines defined in
// quaternionlevelsolver3d.f

void SET_J_IJ3D(const int &, const int &, const int &, const int &, const int &,
                const int &, const double *, const int &, const int &,
                const int &, const int &, const int &, const int &,
                const double *, const int &, const int &, const int &,
                const int &, const int &, const int &, const double *,
                const int &, const int &, const int &, const int &, const int &,
                const int &, const double *, const double *, const int &,
                const int &, const int &, const int &, const int &, const int &,
                const double *, const int &, const int &, const int &,
                const int &, const int &, const int &, const double *,
                const int &, const int &, const int &, const int &, const int &,
                const int &, const double *, const int &, const int &,
                const int &, const int &, const int &, const int &);
void SET_STENCIL3D(const int &, const int &, const int &, const int &,
                   const int &, const int &, const double &, const double *,
                   const int &, const int &, const int &, const int &,
                   const int &, const int &, const double *, const int &,
                   const int &, const int &, const int &, const int &,
                   const int &, const double *, const int &, const int &,
                   const int &, const int &, const int &, const int &,
                   const double *, const int &, const int &, const int &,
                   const int &, const int &, const int &, const double *,
                   const int &, const int &, const int &, const int &,
                   const int &, const int &, double *);
void SET_SYMMETRIC_STENCIL3D(
    const int &, const int &, const int &, const int &, const int &,
    const int &, const double &, const double *, const int &, const int &,
    const int &, const int &, const int &, const int &, const double *,
    const int &, const int &, const int &, const int &, const int &,
    const int &, const double *, const int &, const int &, const int &,
    const int &, const int &, const int &, const double *, const int &,
    const int &, const int &, const int &, const int &, const int &,
    const double *, const int &, const int &, const int &, const int &,
    const int &, const int &, double *);
void COPY3D(const int &, const int &, const int &, const int &, const int &,
            const int &, const double *, const int &, const int &, const int &,
            const int &, const int &, const int &, const double *, const int &,
            const int &, const int &, const int &, const int &, const int &);
void ADJUST_BDRY3D(const double &, const double *, const int &, const int &,
                   const int &, const int &, const int &, const int &,
                   double *diag, const double *offdiagi, const double *offdiagj,
                   const double *offdiagk, const int *pifirst,
                   const int *pilast, const int *pjfirst, const int *pjlast,
                   const int *pkfirst, const int *pklast, const double *acoef,
                   const double *bcoef, const int *aifirst, const int *ailast,
                   const int *ajfirst, const int *ajlast, const int *akfirst,
                   const int *aklast, const double *Ak0, const int *kifirst,
                   const int *kilast, const int *kjfirst, const int *kjlast,
                   const int *kkfirst, const int *kklast, const int *lower,
                   const int *upper, const int *location, const double *h);
void ADJUST_QRHS3D(double *rhs, const int *rifirst, const int *rilast,
                   const int *rjfirst, const int *rjlast, const int *rkfirst,
                   const int *rklast, const double *Ak0, const int *kifirst,
                   const int *kilast, const int *kjfirst, const int *kjlast,
                   const int *kkfirst, const int *kklast, const double *gcoef,
                   const int *aifirst, const int *ailast, const int *ajfirst,
                   const int *ajlast, const int *akfirst, const int *aklast,
                   const int *lower, const int *upper, const int *location);
#endif
}


std::shared_ptr<pdat::OutersideVariable<double> >
    QuatLevelSolver::s_Ak0_var[NDIM];

/*
*************************************************************************
* Constructor                                                           *
*************************************************************************
*/


QuatLevelSolver::QuatLevelSolver(const int ql, const std::string &object_name,
                                 std::shared_ptr<tbox::Database> database)
    : d_dim(tbox::Dimension(NDIM)),
      d_qlen(ql),
      d_object_name(object_name),
      d_ln(-1),
      d_context(hier::VariableDatabase::getDatabase()->getContext(object_name +
                                                                  "::context")),
      d_physical_bc_coef_strategy(&d_physical_bc_simple_case),
      d_physical_bc_simple_case(tbox::Dimension(NDIM),
                                d_object_name + "::simple bc"),
      d_cf_bc_coef(tbox::Dimension(NDIM),
                   object_name + "::coarse-fine bc coefs"),
      d_Ak0_id(-1),
      d_max_iterations(10),
      d_relative_residual_tol(1e-10),
      d_num_pre_relax_steps(1),
      d_num_post_relax_steps(1),
      d_number_iterations(-1),
      d_relative_residual_norm(-1.0),
      d_verbose(false),
      d_solver_type(0),
      d_hypre_object_type(-1),
      d_hypre_data_allocated(false),
      d_print_solver_info(false)
{
   if (database) {
      getFromInput(database);
   }

   t_set_matrix_coefficients = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatLevelSolver::setMatrixCoefficients");
   t_solve_system = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatLevelSolver::solveSystem");
   t_create_hypre_solver = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatLevelSolver::createHypreSolver");
   t_hypre_solve = tbox::TimerManager::getManager()->getTimer(
       "AMPE::QuatLevelSolver::hypreSolve");

   hier::VariableDatabase *vdb = hier::VariableDatabase::getDatabase();
   if (!s_Ak0_var[d_dim.getValue() - 1]) {
      s_Ak0_var[d_dim.getValue() - 1].reset(
          new pdat::OutersideVariable<double>(d_dim, d_object_name + "::Ak0",
                                              d_qlen));
   }
   d_Ak0_id = vdb->registerVariableAndContext(s_Ak0_var[d_dim.getValue() - 1],
                                              d_context,
                                              hier::IntVector::getZero(d_dim));
}


// Set solver from database
// Note that tolerance and max iterations are set from QuatFACOps,
// and depend if the solve is a coarse or fine level solve.

void QuatLevelSolver::getFromInput(std::shared_ptr<tbox::Database> database)
{
   d_print_solver_info =
       database->getBoolWithDefault("print_solver_info", d_print_solver_info);
   d_solver_type = database->getIntegerWithDefault("solver_id", d_solver_type);
   d_num_pre_relax_steps =
       database->getIntegerWithDefault("num_pre_relax_steps",
                                       d_num_pre_relax_steps);
   if (d_num_pre_relax_steps < 0) {
      TBOX_ERROR(d_object_name << ": Number of relaxation steps must be\n"
                               << "non-negative.\n");
   }
   d_num_post_relax_steps =
       database->getIntegerWithDefault("num_post_relax_steps",
                                       d_num_post_relax_steps);
   if (d_num_post_relax_steps < 0) {
      TBOX_ERROR(d_object_name << ": Number of relaxation steps must be\n"
                               << "non-negative.\n");
   }

   // print database just read
   tbox::plog << "QuatLevelSolver database..." << std::endl;
   database->printClassData(tbox::plog);
}


/*
********************************************************************
* Initialize internal data for a given hierarchy level             *
* After setting internal data, propagate the information           *
* to the major algorithm objects.  Allocate data for               *
* storing boundary condition-dependent quantities for              *
* adding to source term before solving.                            *
********************************************************************
*/


void QuatLevelSolver::initializeSolverState(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, int ln)
{
   TBOX_ASSERT(hierarchy);

   deallocateSolverState();

   d_hierarchy = hierarchy;
   d_ln = ln;

   hier::IntVector max_gcw(tbox::Dimension(NDIM), 1);
   d_cf_boundary.reset(
       new hier::CoarseFineBoundary(*d_hierarchy, d_ln, max_gcw));

   d_physical_bc_simple_case.setHierarchy(d_hierarchy, d_ln, d_ln);

   d_number_iterations = -1;
   d_relative_residual_norm = -1.0;

   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(d_ln);
   level->allocatePatchData(d_Ak0_id);
   allocateHypreData();
}


/*
********************************************************************
* Deallocate data initialized by initializeSolverState             *
********************************************************************
*/


void QuatLevelSolver::deallocateSolverState()
{
   if (!d_hierarchy) {
      return;
   }

   d_cf_boundary->clear();
   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(d_ln);
   level->deallocatePatchData(d_Ak0_id);
   deallocateHypreData();
   d_hierarchy.reset();
   d_ln = -1;
}


/*
*************************************************************************
*                                                                       *
* Allocate the HYPRE data structures that depend only on the level      *
* and will not change (grid, stencil, matrix, and std::vectors).             *
*                                                                       *
*************************************************************************
*/

void QuatLevelSolver::allocateHypreData()
{
   if (d_hypre_data_allocated) {
      deallocateHypreData();
   }
   tbox::SAMRAI_MPI::Comm communicator =
       d_hierarchy->getMPI().getCommunicator();

   /*
    * Set up the grid data - only set grid data for local boxes
    */

   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(d_ln);
   std::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry,
                              hier::BaseGridGeometry>(
           d_hierarchy->getGridGeometry()));
   const hier::IntVector ratio = level->getRatioToLevelZero();
   hier::IntVector periodic_shift = grid_geometry->getPeriodicShift(ratio);

   int periodic_flag[NDIM];
   bool is_periodic = false;
   for (int d = 0; d < NDIM; ++d) {
      periodic_flag[d] = periodic_shift[d] != 0;
      is_periodic = is_periodic || periodic_flag[d];
   }


   d_matrix = new HYPRE_SStructMatrix[d_qlen];
   d_linear_rhs = new HYPRE_SStructVector[d_qlen];
   d_linear_sol = new HYPRE_SStructVector[d_qlen];
   if (d_matrix == NULL || d_linear_rhs == NULL || d_linear_sol == NULL) {
      TBOX_ERROR(d_object_name << ": Could not allocate Hypre data\n");
   }

   /* Set the object type (by default HYPRE_SSTRUCT). This determines the
      data structure used to store the matrix.  If you want to use
      unstructured solvers, e.g. BoomerAMG, the object type should be
      HYPRE_PARCSR. If the problem is purely structured (with one part), you
      may want to use HYPRE_STRUCT to access the structured solvers.  */
   if (d_solver_type > 1 && d_solver_type < 4) {
      d_hypre_object_type = HYPRE_PARCSR;

      d_parcsr_solver = new HYPRE_Solver[d_qlen];
      d_parcsr_precond = new HYPRE_Solver[d_qlen];
      d_parcsr_A = new HYPRE_ParCSRMatrix[d_qlen];
      d_parcsr_b = new HYPRE_ParVector[d_qlen];
      d_parcsr_x = new HYPRE_ParVector[d_qlen];

      for (int var = 0; var < d_qlen; var++) {
         d_parcsr_solver[var] = d_parcsr_precond[var] = (HYPRE_Solver)NULL;
         d_parcsr_A[var] = (HYPRE_ParCSRMatrix)NULL;
         d_parcsr_b[var] = d_parcsr_x[var] = (HYPRE_ParVector)NULL;
      }
   } else {
      d_hypre_object_type = HYPRE_STRUCT;

      d_struct_solver = new HYPRE_StructSolver[d_qlen];
      d_struct_precond = new HYPRE_StructSolver[d_qlen];
      d_struct_A = new HYPRE_StructMatrix[d_qlen];
      d_struct_b = new HYPRE_StructVector[d_qlen];
      d_struct_x = new HYPRE_StructVector[d_qlen];

      for (int var = 0; var < d_qlen; var++) {
         d_struct_solver[var] = (HYPRE_StructSolver)NULL;
         d_struct_precond[var] = (HYPRE_StructSolver)NULL;
         d_struct_A[var] = (HYPRE_StructMatrix)NULL;
         d_struct_b[var] = (HYPRE_StructVector)NULL;
         d_struct_x[var] = (HYPRE_StructVector)NULL;
      }
   }

   for (int var = 0; var < d_qlen; var++) {
      d_matrix[var] = (HYPRE_SStructMatrix)NULL;
      d_linear_rhs[var] = (HYPRE_SStructVector)NULL;
      d_linear_sol[var] = (HYPRE_SStructVector)NULL;
   }

   // Create an empty grid object
   HYPRE_SStructGridCreate(communicator, NDIM, NPARTS, &d_grid);

   // Add boxes to the grid
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); p++) {
      const hier::Box &box = p->getBox();
      hier::Index lower = box.lower();
      hier::Index upper = box.upper();
      HYPRE_SStructGridSetExtents(d_grid, PART, &lower[0], &upper[0]);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   tbox::Dimension::dir_t d;
   if (is_periodic) {
      const hier::BoxContainer &level_domain =
          level->getPhysicalDomain(hier::BlockId::zero());
      hier::Box domain_bound(level_domain.front());
      for (hier::BoxContainer::const_iterator i = level_domain.begin();
           i != level_domain.end(); ++i) {
         domain_bound.setLower(
             hier::Index::min(domain_bound.lower(), i->lower()));
         domain_bound.setUpper(
             hier::Index::min(domain_bound.upper(), i->upper()));
      }
      for (d = 0; d < d_dim.getValue(); ++d) {
         if (periodic_flag[d] == true) {
            int tmpi = 1;
            unsigned int p_of_two;
            for (p_of_two = 0; p_of_two < 8 * sizeof(p_of_two) - 1;
                 ++p_of_two) {
               if (tmpi == domain_bound.numberCells(d)) {
                  break;
               }
               if (tmpi > domain_bound.numberCells(d)) {
                  TBOX_ERROR(d_object_name
                             << ": Hypre currently requires\n"
                             << "that grid size in periodic directions be\n"
                             << "powers of two.  (This requirement may go\n"
                             << "away in future versions of hypre.)\n"
                             << "Size problem in direction " << d << "\n"
                             << "Domain bound is " << domain_bound << ",\n"
                             << "Size of " << domain_bound.numberCells()
                             << "\n");
               }
               tmpi = tmpi ? tmpi << 1 : 1;
            }
         }
      }
   }
#endif

   HYPRE_SStructGridSetPeriodic(d_grid, PART, &periodic_shift[0]);

   // Set the variable type and number of variables on each part.
   HYPRE_SStructVariable *vartypes = new HYPRE_SStructVariable[1];
   if (vartypes == NULL) {
      TBOX_ERROR(d_object_name << ": Could not allocate vartypes\n");
   }

   vartypes[0] = HYPRE_SSTRUCT_VARIABLE_CELL;

   HYPRE_SStructGridSetVariables(d_grid, PART, 1, vartypes);

   delete[] vartypes;

   /*
      This is a collective call finalizing the grid assembly.
      The grid is now ``ready to be used''
   */
   HYPRE_SStructGridAssemble(d_grid);

   // Create the graph object
   HYPRE_SStructGraphCreate(communicator, d_grid, &d_graph);
   HYPRE_SStructGraphSetObjectType(d_graph, d_hypre_object_type);

   int i, j;

#ifdef SYMMETRIC_STENCIL
   /*
     if (NDIM==2) {
     int laplacian_offsets[3][NDIM] = {{0,0}, {-1,0}, {0,-1}};
     }
     if (NDIM==3) {
     int laplacian_offsets[4][NDIM] = {{0,0,0}, {-1,0,0}, {0,-1,0}, {0,0,-1}};
     }
   */
   int laplacian_offsets[NDIM + 1][NDIM];
   for (j = 0; j < NDIM; j++) {
      laplacian_offsets[0][j] = 0;  // center
   }
   for (i = 1; i < 2 * NDIM + 1; i += 2) {
      for (j = 0; j < NDIM; j++) {
         laplacian_offsets[(i + 1) / 2][j] =
             (j == (i + 1) / 2 - 1) ? 2 * ((i + 1) % 2) - 1 : 0;
      }
   }

   int stencil_size = NDIM + 1;
#else
   /*
     if (NDIM==2) {
     int laplacian_offsets[5][NDIM] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
     }
     if (NDIM==3) {
     int laplacian_offsets[7][NDIM] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0},
     {0,1,0}, {0,0,-1}, {0,0,1}};
     }
   */
   int laplacian_offsets[2 * NDIM + 1][NDIM];
   for (j = 0; j < NDIM; j++) {
      laplacian_offsets[0][j] = 0;  // center
   }
   for (i = 1; i < 2 * NDIM + 1; i++) {
      for (j = 0; j < NDIM; j++) {
         laplacian_offsets[i][j] =
             (j == (i + 1) / 2 - 1) ? 2 * ((i + 1) % 2) - 1 : 0;
      }
   }

   int stencil_size = 2 * NDIM + 1;
#endif

   HYPRE_SStructStencilCreate(NDIM, stencil_size, &d_stencil);

   for (int entry = 0; entry < stencil_size; entry++) {
      HYPRE_SStructStencilSetEntry(d_stencil, entry, laplacian_offsets[entry],
                                   0);
   }

   // Assign the stencils we created to variables
   HYPRE_SStructGraphSetStencil(d_graph, PART, 0, d_stencil);

   // Assemble the graph
   HYPRE_SStructGraphAssemble(d_graph);

   for (int depth = 0; depth < d_qlen; depth++) {

      // Create empty matrix objects

      HYPRE_SStructMatrixCreate(communicator, d_graph, &d_matrix[depth]);

#ifdef SYMMETRIC_STENCIL
      HYPRE_SStructMatrixSetSymmetric(d_matrix[depth], PART, 0, 0, 1);
#endif

      HYPRE_SStructMatrixSetObjectType(d_matrix[depth], d_hypre_object_type);

      // Indicate that the matrix coefficients are ready to be set
      HYPRE_SStructMatrixInitialize(d_matrix[depth]);

      // Create the rhs and solution std::vectors

      // Create empty std::vector objects
      HYPRE_SStructVectorCreate(communicator, d_grid, &d_linear_rhs[depth]);
      HYPRE_SStructVectorCreate(communicator, d_grid, &d_linear_sol[depth]);

      /* Set the object type for the std::vectors
         to be the same as was already set for the matrix */
      HYPRE_SStructVectorSetObjectType(d_linear_rhs[depth],
                                       d_hypre_object_type);
      HYPRE_SStructVectorSetObjectType(d_linear_sol[depth],
                                       d_hypre_object_type);

      // Indicate that the std::vector coefficients are ready to be set
      HYPRE_SStructVectorInitialize(d_linear_rhs[depth]);
      HYPRE_SStructVectorInitialize(d_linear_sol[depth]);
   }

   d_hypre_data_allocated = true;
}


/*
*************************************************************************
*                                                                       *
* The destructor deallocates solver data.                               *
*                                                                       *
*************************************************************************
*/


QuatLevelSolver::~QuatLevelSolver()
{
   deallocateHypreData();

   if (d_hierarchy) {
      std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(0);
      level->deallocatePatchData(d_Ak0_id);
   }
   hier::VariableDatabase *vdb = hier::VariableDatabase::getDatabase();
   vdb->removePatchDataIndex(d_Ak0_id);
}


/*
*************************************************************************
*                                                                       *
* Deallocate HYPRE solver and data.  HYPRE requires that we             *
* check whether HYPRE has already deallocated this data.                *
* Note that the HYPRE solver, d_solver, was created at                  *
* the end of setMatrixCoefficients.                                     *
*                                                                       *
*************************************************************************
*/


void QuatLevelSolver::deallocateHypreData()
{
   if (d_hypre_data_allocated) {

      // Destroy the Hypre objects

      int depth;
      for (depth = 0; depth < d_qlen; depth++) {
         destroyHypreSolver(depth);
         HYPRE_SStructVectorDestroy(d_linear_sol[depth]);
         HYPRE_SStructVectorDestroy(d_linear_rhs[depth]);
         HYPRE_SStructMatrixDestroy(d_matrix[depth]);
      }

      HYPRE_SStructStencilDestroy(d_stencil);

      delete[] d_matrix;
      delete[] d_linear_rhs;
      delete[] d_linear_sol;

      if (d_hypre_object_type == HYPRE_PARCSR) {
         delete[] d_parcsr_solver;
         delete[] d_parcsr_precond;
         delete[] d_parcsr_A;
         delete[] d_parcsr_b;
         delete[] d_parcsr_x;
      } else if (d_hypre_object_type == HYPRE_STRUCT) {
         delete[] d_struct_solver;
         delete[] d_struct_precond;
         delete[] d_struct_A;
         delete[] d_struct_b;
         delete[] d_struct_x;
      }

      HYPRE_SStructGraphDestroy(d_graph);
      HYPRE_SStructGridDestroy(d_grid);

      d_hypre_data_allocated = false;
   }
}


/*
*************************************************************************
*                                                                       *
* Copy data into the HYPRE std::vector structures.                           *
*                                                                       *
*************************************************************************
*/

void QuatLevelSolver::copyToHypre(const int component,
                                  const pdat::CellData<double> &src_data,
                                  HYPRE_SStructVector dst_vector)
{
   // Tracer t("QuatLevelSolver::copyToHypre");
   assert(src_data.getDepth() == d_qlen);


   hier::Box box(src_data.getBox());
   const hier::Box &gbox = src_data.getGhostBox();
   hier::Index lo = box.lower();
   hier::Index up = box.upper();

   double *values = new double[box.size()];
   if (values == NULL) {
      TBOX_ERROR(d_object_name << ": Could not allocate values\n");
   }

#if NDIM == 2
   COPY2D(box.lower(0), box.upper(0), box.lower(1), box.upper(1),
          src_data.getPointer(component), gbox.lower(0), gbox.upper(0),
          gbox.lower(1), gbox.upper(1), values, box.lower(0), box.upper(0),
          box.lower(1), box.upper(1));
#endif
#if NDIM == 3
   COPY3D(box.lower(0), box.upper(0), box.lower(1), box.upper(1), box.lower(2),
          box.upper(2), src_data.getPointer(component), gbox.lower(0),
          gbox.upper(0), gbox.lower(1), gbox.upper(1), gbox.lower(2),
          gbox.upper(2), values, box.lower(0), box.upper(0), box.lower(1),
          box.upper(1), box.lower(2), box.upper(2));
#endif

   HYPRE_SStructVectorSetBoxValues(dst_vector, PART, &lo[0], &up[0], 0, values);

   delete[] values;
}

/*
*************************************************************************
*                                                                       *
* Copy data out of the HYPRE std::vector structures.                         *
*                                                                       *
*************************************************************************
*/

void QuatLevelSolver::copyFromHypre(const int depth,
                                    const HYPRE_SStructVector src_vector,
                                    const int component,
                                    pdat::CellData<double> &dst_data)
{
   // Tracer t("QuatLevelSolver::copyFromHypre");
   assert(dst_data.getDepth() == d_qlen);

   hier::Box box(dst_data.getBox());
   const hier::Box &gbox = dst_data.getGhostBox();

   double *values = new double[box.size()];
   if (values == NULL) {
      TBOX_ERROR(d_object_name << ": Could not allocate values\n");
   }

   hier::Index lower = box.lower();
   hier::Index upper = box.upper();
   HYPRE_SStructVectorGetBoxValues(src_vector, PART, &lower[0], &upper[0],
                                   depth, values);

   // copy from values into dst_data
#if NDIM == 2
   COPY2D(box.lower(0), box.upper(0), box.lower(1), box.upper(1), values,
          box.lower(0), box.upper(0), box.lower(1), box.upper(1),
          dst_data.getPointer(component), gbox.lower(0), gbox.upper(0),
          gbox.lower(1), gbox.upper(1));
#endif
#if NDIM == 3
   COPY3D(box.lower(0), box.upper(0), box.lower(1), box.upper(1), box.lower(2),
          box.upper(2), values, box.lower(0), box.upper(0), box.lower(1),
          box.upper(1), box.lower(2), box.upper(2),
          dst_data.getPointer(component), gbox.lower(0), gbox.upper(0),
          gbox.lower(1), gbox.upper(1), gbox.lower(2), gbox.upper(2));
#endif

   delete[] values;
}


/*
*************************************************************************
*                                                                       *
* Set the matrix coefficients for the linear system.                    *
* The matrix coefficients are dependent on the problem                  *
* specification described by the PoissonSpecificiations                 *
* object and by the boundary condition.                                 *
*                                                                       *
*************************************************************************
*/

void QuatLevelSolver::setMatrixCoefficients(const double gamma,
                                            const int sqrt_mobility_id,
                                            const int face_coef_id)
{
   t_set_matrix_coefficients->start();

   if (d_physical_bc_coef_strategy == NULL) {
      TBOX_ERROR(d_object_name << ": No BC coefficient strategy object!\n"
                               << "Use either setBoundaries or "
                                  "setPhysicalBcCoefObject\n"
                               << "to specify the boundary conidition.  Do it "
                                  "before\n"
                               << "calling setMatrixCoefficients.");
   }

   /* The value of the ghost cell based on the Robin boundary condition
    * can be written as the sum of a constant, k0, plus a multiple of the
    * internal cell value, k1*ui.  k1*ui depends on the value of u so it
    * contributes to the product Au,
    * while the constant k0 contributes the right hand side f.
    * We save Ak0 = A*k0(a) to add to f when solving.
    * We assume unit g here because we will multiply it in just before
    * solving, thus allowing everything that does not affect A to change
    * from solve to solve.
    */
   std::shared_ptr<pdat::OutersideData<double> > Ak0;

   /*
    * Loop over patches and set matrix entries for each patch.
    */

   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(d_ln);

   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        ++pi) {
      hier::Patch &patch = **pi;

      // Get the cell width array
      std::shared_ptr<geom::CartesianPatchGeometry> pg(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                 hier::PatchGeometry>(
              patch.getPatchGeometry()));

      const double *h = pg->getDx();

      /*
        We make a local copy of the patch box so that we can get non-const
        integer pointers to the lower and upper box corners that can be passed
        to Hypre (whose prototypes don't accept const modifiers)
      */
      hier::Box patch_box(patch.getBox());
      hier::Index lower = patch_box.lower();
      hier::Index upper = patch_box.upper();
      int num_cells = patch_box.size();

      Ak0 = std::dynamic_pointer_cast<pdat::OutersideData<double>,
                                      hier::PatchData>(
          patch.getPatchData(d_Ak0_id));
      assert(Ak0->getDepth() == d_qlen);

      Ak0->fillAll(0.0);

      // Storage for diagonal and off-diagonal entries,
      pdat::CellData<double> diagonal(patch_box, d_qlen,
                                      hier::IntVector(tbox::Dimension(NDIM),
                                                      0));
      pdat::SideData<double> off_diagonal(patch_box, d_qlen,
                                          hier::IntVector(tbox::Dimension(NDIM),
                                                          0));

      std::shared_ptr<pdat::SideData<double> > face_coef_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch.getPatchData(face_coef_id)));
      //  Set the J_ij blocks

      for (int depth = 0; depth < d_qlen; depth++) {

#if NDIM == 2
         SET_J_IJ2D(lower[0], upper[0], lower[1], upper[1],
                    face_coef_data->getPointer(0), lower[0], upper[0] + 1,
                    lower[1], upper[1], face_coef_data->getPointer(1), lower[0],
                    upper[0], lower[1], upper[1] + 1, h,
                    diagonal.getPointer(depth), lower[0], upper[0], lower[1],
                    upper[1], off_diagonal.getPointer(0, depth), lower[0],
                    upper[0] + 1, lower[1], upper[1],
                    off_diagonal.getPointer(1, depth), lower[0], upper[0],
                    lower[1], upper[1] + 1);
#endif
#if NDIM == 3
         SET_J_IJ3D(lower[0], upper[0], lower[1], upper[1], lower[2], upper[2],
                    face_coef_data->getPointer(0), lower[0], upper[0] + 1,
                    lower[1], upper[1], lower[2], upper[2],
                    face_coef_data->getPointer(1), lower[0], upper[0], lower[1],
                    upper[1] + 1, lower[2], upper[2],
                    face_coef_data->getPointer(2), lower[0], upper[0], lower[1],
                    upper[1], lower[2], upper[2] + 1, h,
                    diagonal.getPointer(depth), lower[0], upper[0], lower[1],
                    upper[1], lower[2], upper[2],
                    off_diagonal.getPointer(0, depth), lower[0], upper[0] + 1,
                    lower[1], upper[1], lower[2], upper[2],
                    off_diagonal.getPointer(1, depth), lower[0], upper[0],
                    lower[1], upper[1] + 1, lower[2], upper[2],
                    off_diagonal.getPointer(2, depth), lower[0], upper[0],
                    lower[1], upper[1], lower[2], upper[2] + 1);
#endif
      }

      std::shared_ptr<pdat::CellData<double> > sqrt_mobility_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(sqrt_mobility_id)));

      const hier::Box &m_gbox = sqrt_mobility_data->getGhostBox();
      const hier::Index &mlower = m_gbox.lower();
      const hier::Index &mupper = m_gbox.upper();

      /*
       * Walk physical domain boundaries and adjust off-diagonals
       * before computation of diagonal entries.
       * The exterior cell's value is
       * uo = ( h*gamma + ui*(beta-h*alpha/2) )/( beta+h*alpha/2 )
       *   = k0 + k1*ui
       * where k0 = h*gamma/( beta+h*alpha/2 )
       * k1 = ( beta-h*alpha/2 )/( beta+h*alpha/2 )
       * Split coupling between interior-exterior cells
       * into two parts: interior-interior coupling (k1)
       * and rhs contribution (k0).
       */
      {
         const std::vector<hier::BoundaryBox> &surface_boxes =
             pg->getCodimensionBoundaries(1);
         const int n_bdry_boxes = static_cast<int>(surface_boxes.size());
         for (int n = 0; n < n_bdry_boxes; ++n) {

            const hier::BoundaryBox &boundary_box = surface_boxes[n];
            if (boundary_box.getBoundaryType() != 1) {
               TBOX_ERROR(d_object_name << ": Illegal boundary type in "
                                        << "QuatLevelSolver::"
                                           "setMatrixCoefficients\n");
            }
            const hier::BoundaryBoxUtils bbu(boundary_box);
            const hier::BoundaryBox trimmed_boundary_box =
                bbu.trimBoundaryBox(patch.getBox());
            const hier::Box bccoef_box = bbu.getSurfaceBoxFromBoundaryBox();
            std::shared_ptr<pdat::ArrayData<double> > acoef_data(
                std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1));
            std::shared_ptr<pdat::ArrayData<double> > bcoef_data(
                std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1));
            std::shared_ptr<pdat::ArrayData<double> > gcoef_data;
            d_physical_bc_coef_strategy->setBcCoefs(acoef_data, bcoef_data,
                                                    gcoef_data,
                                                    d_physical_bc_variable,
                                                    patch, boundary_box);
            adjustBoundaryEntries(gamma, *sqrt_mobility_data, diagonal,
                                  off_diagonal, patch_box, *acoef_data,
                                  *bcoef_data, bccoef_box, *Ak0,
                                  trimmed_boundary_box, h);
         }
      }

      /*
       * Walk coarse-fine boundaries and adjust off-diagonals
       * according data in ghost cells.
       */
      if (d_ln > 0) {
         /*
          * There are potentially coarse-fine boundaries to deal with.
          */

         std::vector<hier::BoundaryBox> surface_boxes;

         if (NDIM == 2) {
            surface_boxes = d_cf_boundary->getEdgeBoundaries(pi->getGlobalId());
         } else if (NDIM == 3) {
            surface_boxes = d_cf_boundary->getFaceBoundaries(pi->getGlobalId());
         }

         const int n_bdry_boxes = surface_boxes.size();
         for (int n = 0; n < n_bdry_boxes; ++n) {

            const hier::BoundaryBox &boundary_box = surface_boxes[n];
            if (boundary_box.getBoundaryType() != 1) {
               TBOX_ERROR(d_object_name << ": Illegal boundary type in "
                                        << "QuatLevelSolver::"
                                           "setMatrixCoefficients\n");
            }
            const hier::BoundaryBoxUtils bbu(boundary_box);
            const hier::BoundaryBox trimmed_boundary_box =
                bbu.trimBoundaryBox(patch.getBox());
            const hier::Box bccoef_box = bbu.getSurfaceBoxFromBoundaryBox();
            std::shared_ptr<pdat::ArrayData<double> > acoef_data(
                std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1));
            std::shared_ptr<pdat::ArrayData<double> > bcoef_data(
                std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1));
            std::shared_ptr<pdat::ArrayData<double> > gcoef_data;

            /*
           Reset invalid ghost data id to help detect use in setBcCoefs.
         */
            d_cf_bc_coef.setGhostDataId(-1,
                                        hier::IntVector(tbox::Dimension(NDIM),
                                                        0));
            d_cf_bc_coef.setBcCoefs(acoef_data, bcoef_data, gcoef_data,
                                    d_coarsefine_bc_variable, patch,
                                    boundary_box);

            adjustBoundaryEntries(gamma, *sqrt_mobility_data, diagonal,
                                  off_diagonal, patch_box, *acoef_data,
                                  *bcoef_data, bccoef_box, *Ak0,
                                  trimmed_boundary_box, h);
         }
      }

#ifdef SYMMETRIC_STENCIL
      int nentries = NDIM + 1;
#else
      int nentries = 2 * NDIM + 1;
#endif

      int nvalues = nentries * num_cells;
      double *values = new double[nvalues];
      if (values == NULL) {
         TBOX_ERROR(d_object_name << ": Could not allocate values\n");
      }

      for (int depth = 0; depth < d_qlen; depth++) {

#ifdef SYMMETRIC_STENCIL

#if NDIM == 2
         SET_SYMMETRIC_STENCIL2D(lower[0], upper[0], lower[1], upper[1], gamma,
                                 sqrt_mobility_data->getPointer(), mlower[0],
                                 mupper[0], mlower[1], mupper[1],
                                 diagonal.getPointer(depth), lower[0], upper[0],
                                 lower[1], upper[1],
                                 off_diagonal.getPointer(0, depth), lower[0],
                                 upper[0] + 1, lower[1], upper[1],
                                 off_diagonal.getPointer(1, depth), lower[0],
                                 upper[0], lower[1], upper[1] + 1, values);
#endif
#if NDIM == 3
         SET_SYMMETRIC_STENCIL3D(
             lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], gamma,
             sqrt_mobility_data->getPointer(), mlower[0], mupper[0], mlower[1],
             mupper[1], mlower[2], mupper[2], diagonal.getPointer(depth),
             lower[0], upper[0], lower[1], upper[1], lower[2], upper[2],
             off_diagonal.getPointer(0, depth), lower[0], upper[0] + 1,
             lower[1], upper[1], lower[2], upper[2],
             off_diagonal.getPointer(1, depth), lower[0], upper[0], lower[1],
             upper[1] + 1, lower[2], upper[2],
             off_diagonal.getPointer(2, depth), lower[0], upper[0], lower[1],
             upper[1], lower[2], upper[2] + 1, values);
#endif

         int intra_q_indices[NDIM + 1];

#else

#if NDIM == 2
         SET_STENCIL2D(lower[0], upper[0], lower[1], upper[1], gamma,
                       sqrt_mobility_data->getPointer(), mlower[0], mupper[0],
                       mlower[1], mupper[1], diagonal.getPointer(depth),
                       lower[0], upper[0], lower[1], upper[1],
                       off_diagonal.getPointer(0, depth), lower[0],
                       upper[0] + 1, lower[1], upper[1],
                       off_diagonal.getPointer(1, depth), lower[0], upper[0],
                       lower[1], upper[1] + 1, values);
#endif
#if NDIM == 3
         SET_STENCIL3D(lower[0], upper[0], lower[1], upper[1], lower[2],
                       upper[2], gamma, sqrt_mobility_data->getPointer(),
                       mlower[0], mupper[0], mlower[1], mupper[1], mlower[2],
                       mupper[2], diagonal.getPointer(depth), lower[0],
                       upper[0], lower[1], upper[1], lower[2], upper[2],
                       off_diagonal.getPointer(0, depth), lower[0],
                       upper[0] + 1, lower[1], upper[1], lower[2], upper[2],
                       off_diagonal.getPointer(1, depth), lower[0], upper[0],
                       lower[1], upper[1] + 1, lower[2], upper[2],
                       off_diagonal.getPointer(2, depth), lower[0], upper[0],
                       lower[1], upper[1], lower[2], upper[2] + 1, values);
#endif

         int intra_q_indices[2 * NDIM + 1];

#endif

         for (int i = 0; i < nentries; i++)
            intra_q_indices[i] = i;

         HYPRE_SStructMatrixSetBoxValues(d_matrix[depth], PART, &lower[0],
                                         &upper[0], 0, nentries,
                                         intra_q_indices, values);
      }

      delete[] values;
      values = NULL;

   }  // end patch loop

   if (d_print_solver_info)  // dump before assembly
      for (int depth = 0; depth < d_qlen; depth++) {
         std::string name("mat_bA.out_q" + std::to_string(depth));
         HYPRE_SStructMatrixPrint(name.c_str(), d_matrix[depth], 1);
      }

   /* This is a collective call finalizing the matrix assembly.
      The matrix is now ``ready to be used'' */
   for (int depth = 0; depth < d_qlen; depth++) {
      HYPRE_SStructMatrixAssemble(d_matrix[depth]);

      if (d_print_solver_info) {  // dump after assembly
         std::string name("mat_aA.out_q" + std::to_string(depth));
         HYPRE_SStructMatrixPrint(name.c_str(), d_matrix[depth], 1);
      }

      createHypreSolver(depth);
   }

   t_set_matrix_coefficients->stop();
}


/*
**********************************************************************
* Add g*A*k0(a) from physical boundaries to rhs.                     *
* This operation is done for physical as well as cf boundaries,      *
* so it is placed in a function.                                     *
**********************************************************************
*/

void QuatLevelSolver::add_gAk0_toRhs(
    const hier::Patch &patch, const std::vector<hier::BoundaryBox> &bdry_boxes,
    const solv::RobinBcCoefStrategy *robin_bc_coef, pdat::CellData<double> &rhs)
{
   assert(rhs.getDepth() == d_qlen);
   const hier::Box &rhsbox = rhs.getArrayData().getBox();

#ifdef DEBUG_CHECK_ASSERTIONS
   math::ArrayDataNormOpsReal<double> ops;
   for (int depth = 0; depth < d_qlen; depth++) {
      double nb = ops.maxNorm(rhs.getArrayData(), rhsbox);
      assert(nb == nb);
   }
#endif

   /*
    * g*A*k0(a) is the storage for adjustments to be made to the rhs
    * when we solve. This is the value of the weight of the ghost cell
    * value for the interior cell, times k0.  It is independent of u,
    * and so is moved to the rhs.  Before solving, g*A*k0(a) is added
    * to rhs.
    */
   std::shared_ptr<pdat::OutersideData<double> > Ak0(
       SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>, hier::PatchData>(
           patch.getPatchData(d_Ak0_id)));
   assert(Ak0->getDepth() == d_qlen);

   const int n_bdry_boxes = bdry_boxes.size();
   for (int n = 0; n < n_bdry_boxes; ++n) {

      const hier::BoundaryBox &boundary_box = bdry_boxes[n];
#ifdef DEBUG_CHECK_ASSERTIONS
      if (boundary_box.getBoundaryType() != 1) {
         TBOX_ERROR(d_object_name << ": Illegal boundary type in "
                                  << "QuatLevelSolver::add_gAk0_toRhs\n");
      }
#endif
      const int location_index = boundary_box.getLocationIndex();
      const hier::BoundaryBoxUtils bbu(boundary_box);
      const hier::BoundaryBox trimmed_boundary_box =
          bbu.trimBoundaryBox(patch.getBox());
      const hier::Index &lower = trimmed_boundary_box.getBox().lower();
      const hier::Index &upper = trimmed_boundary_box.getBox().upper();
      const hier::Box &Ak0box =
          Ak0->getArrayData(location_index / 2, location_index % 2).getBox();
      const hier::Box bccoef_box = bbu.getSurfaceBoxFromBoundaryBox();
      std::shared_ptr<pdat::ArrayData<double> > acoef_data;
      std::shared_ptr<pdat::ArrayData<double> > bcoef_data;
      std::shared_ptr<pdat::ArrayData<double> > gcoef_data(
          std::make_shared<pdat::ArrayData<double> >(bccoef_box, 1));
      static const double fill_time = 0.0;
      // Set bc coefficients
      robin_bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data,
                                d_physical_bc_variable, patch, boundary_box,
                                fill_time);

#ifdef DEBUG_CHECK_ASSERTIONS
      pdat::ArrayData<double> &Ak0_data =
          Ak0->getArrayData(location_index / 2, location_index % 2);
      double nak0 = ops.maxNorm(Ak0_data, Ak0box);
      if (nak0 != nak0) {
         Ak0->print(Ak0box, std::cerr);
         std::cerr << "Ak0: Norm =" << nak0 << std::endl;
         TBOX_ERROR("Nan in Ak0 before adjusting rhs\n");
      }
      if (nak0 > 1.e25) {
         Ak0->print(Ak0box, std::cerr);
         std::cerr << "Ak0: Norm =" << nak0 << std::endl;
         TBOX_ERROR("Inf in Ak0 before adjusting rhs\n");
      }
#endif

      for (int depth = 0; depth < d_qlen; depth++) {
#ifdef DEBUG_CHECK_ASSERTIONS
         double nb = ops.maxNorm(rhs.getArrayData(), rhsbox);
         assert(nb == nb);
#endif

#if NDIM == 2
         ADJUST_QRHS2D(
             rhs.getPointer(depth), &rhsbox.lower()[0], &rhsbox.upper()[0],
             &rhsbox.lower()[1], &rhsbox.upper()[1],
             Ak0->getPointer(location_index / 2, location_index % 2, depth),
             &Ak0box.lower()[0], &Ak0box.upper()[0], &Ak0box.lower()[1],
             &Ak0box.upper()[1], gcoef_data->getPointer(),
             &bccoef_box.lower()[0], &bccoef_box.upper()[0],
             &bccoef_box.lower()[1], &bccoef_box.upper()[1], &lower[0],
             &upper[0], &location_index);
#endif
#if NDIM == 3
         ADJUST_QRHS3D(
             rhs.getPointer(depth), &rhsbox.lower()[0], &rhsbox.upper()[0],
             &rhsbox.lower()[1], &rhsbox.upper()[1], &rhsbox.lower()[2],
             &rhsbox.upper()[2],
             Ak0->getPointer(location_index / 2, location_index % 2, depth),
             &Ak0box.lower()[0], &Ak0box.upper()[0], &Ak0box.lower()[1],
             &Ak0box.upper()[1], &Ak0box.lower()[2], &Ak0box.upper()[2],
             gcoef_data->getPointer(), &bccoef_box.lower()[0],
             &bccoef_box.upper()[0], &bccoef_box.lower()[1],
             &bccoef_box.upper()[1], &bccoef_box.lower()[2],
             &bccoef_box.upper()[2], &lower[0], &upper[0], &location_index);
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
         nak0 = ops.maxNorm(Ak0_data, Ak0box);
         if (nak0 != nak0) {
            Ak0->print(Ak0box, std::cerr);
            std::cerr << "Ak0: Norm =" << nak0 << std::endl;
            TBOX_ERROR("Nan in Ak0 after adjusting rhs\n");
         }
         double na = ops.maxNorm(rhs.getArrayData(), rhsbox);
         if (na != na) {
            Ak0->print(Ak0box, std::cerr);
            std::cerr << "Norm before=" << nb << std::endl;
            std::cerr << "Norm after=" << na << std::endl;
            TBOX_ERROR("Nan in rhs\n");
         }
#endif
      }
   }
}


/*
*************************************************************************
* Create the hypre solver and set it according to the current state.    *
*************************************************************************
*/
void QuatLevelSolver::createHypreSolver(int depth)
{
   t_create_hypre_solver->start();

#ifdef HAVE_MPI
   const tbox::SAMRAI_MPI &mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   MPI_Comm communicator = mpi.getCommunicator();
#else
   MPI_Comm communicator;
#endif

   // If we are using the parCSR or structured solvers, get the appropriate
   // matrix, solution and rhs objects
   if (d_hypre_object_type == HYPRE_PARCSR) {
      HYPRE_SStructMatrixGetObject(d_matrix[depth], (void **)&d_parcsr_A);
      HYPRE_SStructVectorGetObject(d_linear_rhs[depth], (void **)&d_parcsr_b);
      HYPRE_SStructVectorGetObject(d_linear_sol[depth], (void **)&d_parcsr_x);
   } else if (d_hypre_object_type == HYPRE_STRUCT) {
      HYPRE_SStructMatrixGetObject(d_matrix[depth],
                                   (void **)&d_struct_A[depth]);
      HYPRE_SStructVectorGetObject(d_linear_rhs[depth],
                                   (void **)&d_struct_b[depth]);
      HYPRE_SStructVectorGetObject(d_linear_sol[depth],
                                   (void **)&d_struct_x[depth]);
   }

   switch (d_solver_type) {

      case 0:  // multigrid preconditioned PCG

         if (d_struct_solver[depth] != NULL) destroyHypreSolver(depth);

         HYPRE_StructPCGCreate(communicator, &d_struct_solver[depth]);

         // PCG parameters
         HYPRE_StructPCGSetPrintLevel(d_struct_solver[depth],
                                      0);  // print each PCG iteration
         HYPRE_StructPCGSetLogging(d_struct_solver[depth], 1);

         // Use PFMG as precondititioner
         HYPRE_StructPFMGCreate(communicator, &d_struct_precond[depth]);

         // Set PFMG parameters
         HYPRE_StructPFMGSetTol(d_struct_precond[depth], 0.);
         HYPRE_StructPFMGSetMaxIter(d_struct_precond[depth], 1);
         HYPRE_StructPFMGSetNumPreRelax(d_struct_precond[depth], 1);
         HYPRE_StructPFMGSetNumPostRelax(d_struct_precond[depth], 1);
         HYPRE_StructPFMGSetPrintLevel(d_struct_precond[depth], 0);
         HYPRE_StructPFMGSetZeroGuess(d_struct_precond[depth]);

         // Set the preconditioner
         HYPRE_StructPCGSetPrecond(d_struct_solver[depth],
                                   HYPRE_StructPFMGSolve, HYPRE_StructPFMGSetup,
                                   d_struct_precond[depth]);

         // Finish the setup
         HYPRE_StructPCGSetup(d_struct_solver[depth], d_struct_A[depth],
                              d_struct_b[depth], d_struct_x[depth]);

         break;
      case 1:  // multigrid only

         if (d_struct_solver[depth] != NULL) destroyHypreSolver(depth);

         HYPRE_StructPFMGCreate(communicator, &d_struct_solver[depth]);

         // Set PFMG parameters
         HYPRE_StructPFMGSetNumPreRelax(d_struct_solver[depth], 1);
         HYPRE_StructPFMGSetNumPostRelax(d_struct_solver[depth], 1);
         HYPRE_StructPFMGSetPrintLevel(d_struct_solver[depth], 1);
         HYPRE_StructPFMGSetLogging(d_struct_solver[depth], 1);

         // Finish the setup
         HYPRE_StructPFMGSetup(d_struct_solver[depth], d_struct_A[depth],
                               d_struct_b[depth], d_struct_x[depth]);

         break;
      case 2:  // AMG preconditioned PCG

         HYPRE_ParCSRPCGCreate(communicator, &d_parcsr_solver[depth]);

         // Set the PCG paramaters
         HYPRE_PCGSetMaxIter(d_parcsr_solver[depth], 100);
         HYPRE_PCGSetTol(d_parcsr_solver[depth], 1.0e-06);
         HYPRE_PCGSetPrintLevel(d_parcsr_solver[depth], 2);
         HYPRE_PCGSetLogging(d_parcsr_solver[depth], 1);

         // Use BoomerAMG as preconditioner
         HYPRE_BoomerAMGCreate(&d_parcsr_precond[depth]);
         HYPRE_BoomerAMGSetCoarsenType(d_parcsr_precond[depth], 6);
         HYPRE_BoomerAMGSetStrongThreshold(d_parcsr_precond[depth], 0.25);
         HYPRE_BoomerAMGSetTol(d_parcsr_precond[depth], 0.0);
         HYPRE_BoomerAMGSetPrintLevel(d_parcsr_precond[depth], 1);
         HYPRE_BoomerAMGSetPrintFileName(d_parcsr_precond[depth],
                                         "ex9.out."
                                         "log");
         HYPRE_BoomerAMGSetMaxIter(d_parcsr_precond[depth], 1);

         // Set the preconditioner
         HYPRE_ParCSRPCGSetPrecond(d_parcsr_solver[depth], HYPRE_BoomerAMGSolve,
                                   HYPRE_BoomerAMGSetup,
                                   d_parcsr_precond[depth]);

         // Finish the setup
         HYPRE_ParCSRPCGSetup(d_parcsr_solver[depth], d_parcsr_A[depth],
                              d_parcsr_b[depth], d_parcsr_x[depth]);

         break;
      case 3:  // AMG only

         HYPRE_BoomerAMGCreate(&d_parcsr_solver[depth]);
         HYPRE_BoomerAMGSetCoarsenType(d_parcsr_solver[depth], 6);
         HYPRE_BoomerAMGSetStrongThreshold(d_parcsr_solver[depth], 0.25);
         HYPRE_BoomerAMGSetTol(d_parcsr_solver[depth], 1.9e-6);
         HYPRE_BoomerAMGSetPrintLevel(d_parcsr_solver[depth], 1);
         HYPRE_BoomerAMGSetPrintFileName(d_parcsr_solver[depth], "ex9.out.log");
         HYPRE_BoomerAMGSetMaxIter(d_parcsr_solver[depth], 50);

         // Finish the setup
         HYPRE_BoomerAMGSetup(d_parcsr_solver[depth], d_parcsr_A[depth],
                              d_parcsr_b[depth], d_parcsr_x[depth]);

         break;
      default:
         TBOX_ERROR(d_object_name << ": ERROR: Invalid solver id specified.\n");
   }

   t_create_hypre_solver->stop();
}


void QuatLevelSolver::destroyHypreSolver(int depth)
{

   switch (d_solver_type) {

      case 0:  // multigrid preconditioned PCG
         if (d_struct_solver[depth])
            HYPRE_StructPCGDestroy(d_struct_solver[depth]);
         if (d_struct_precond[depth])
            HYPRE_StructPFMGDestroy(d_struct_precond[depth]);
         d_struct_solver[depth] = d_struct_precond[depth] =
             (HYPRE_StructSolver)NULL;
         break;
      case 1:  // multigrid only
         if (d_struct_solver[depth])
            HYPRE_StructPFMGDestroy(d_struct_solver[depth]);
         d_struct_solver[depth] = (HYPRE_StructSolver)NULL;
         break;
      case 2:  // AMR preconditioned PCG
         if (d_parcsr_solver[depth])
            HYPRE_ParCSRPCGDestroy(d_parcsr_solver[depth]);
         if (d_parcsr_precond[depth])
            HYPRE_BoomerAMGDestroy(d_parcsr_precond[depth]);
         d_parcsr_solver[depth] = d_parcsr_precond[depth] = (HYPRE_Solver)NULL;
         break;
      case 3:  // AMG only
         if (d_parcsr_solver[depth])
            HYPRE_BoomerAMGDestroy(d_parcsr_solver[depth]);
         d_parcsr_solver[depth] = (HYPRE_Solver)NULL;
         break;
      default:
         TBOX_ERROR(d_object_name << ": ERROR: Invalid solver id specified.\n");
   }
}


/*
*************************************************************************
*                                                                       *
* Solve the linear system.  This routine assumes that the boundary      *
* conditions and the matrix coefficients have been specified.           *
*                                                                       *
*************************************************************************
*/

int QuatLevelSolver::solveSystem(const int q_rhs_id, const int q_solution_id,
                                 bool homogeneous_bc)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(q_rhs_id >= 0);
   assert(q_solution_id >= 0);
   for (int depth = 0; depth < d_qlen; depth++)
      assert(d_struct_solver[depth] != NULL);
#endif

   // Tracer t("QuatLevelSolver::solveSystem");

   t_solve_system->start();

   if (d_physical_bc_coef_strategy == NULL) {
      TBOX_ERROR(d_object_name << ": No BC coefficient strategy object!\n"
                               << "Use either setBoundaries or "
                                  "setPhysicalBcCoefObject\n"
                               << "to specify the boundary conidition.  Do it "
                                  "before\n"
                               << "calling solveSystem.");
   }

   if (d_physical_bc_coef_strategy == &d_physical_bc_simple_case) {
      /*
       * If we are using the simple bc implementation, the final piece
       * of information it requires is the Dirichlet boundary value
       * set in the ghost cells.  Now that we have the ghost cell data,
       * we can complete the boundary condition setup.
       */
      d_physical_bc_simple_case.cacheDirichletData(q_solution_id);
   }

   /*
    * Modify right-hand-side to account for boundary conditions and
    * copy solution and right-hand-side to HYPRE structures.
    */

   /*
    * At coarse-fine boundaries, we expect ghost cells to have correct
    * values to be used in our bc, so u provides the ghost cell data.
    * Assume that the user only provided data for the immediate first
    * ghost cell, so pass zero for the number of extensions fillable.
    */
   d_cf_bc_coef.setGhostDataId(q_solution_id,
                               hier::IntVector(tbox::Dimension(NDIM), 0));

   math::ArrayDataNormOpsReal<double> ops;

   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(d_ln);
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); p++) {
      std::shared_ptr<hier::Patch> patch = *p;

      const hier::Box box = patch->getBox();

      // Get references to the hierarchy solution data
      std::shared_ptr<pdat::CellData<double> > q_solution_data_(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(q_solution_id)));
      TBOX_ASSERT(q_solution_data_);
      TBOX_ASSERT(q_solution_data_->getDepth() == d_qlen);

      pdat::CellData<double> &q_solution_data = *q_solution_data_;
      pdat::CellData<double> q_rhs_data(box, d_qlen,
                                        hier::IntVector(tbox::Dimension(NDIM),
                                                        0));

      // Copy solution from the hierarchy into the HYPRE std::vector
      for (int depth = 0; depth < d_qlen; depth++) {
         copyToHypre(depth, q_solution_data, d_linear_sol[depth]);
         // newcopyToHypre(d_linear_sol[depth], q_solution_data, depth, box );
      }
      q_rhs_data.copy(*(patch->getPatchData(q_rhs_id)));
#ifdef DEBUG_CHECK_ASSERTIONS
      double nb = ops.maxNorm(q_rhs_data.getArrayData(), box);
      assert(nb == nb);
#endif

      if (!homogeneous_bc) {
         // Add g*A*k0(a) from physical and coarse-fine boundaries to rhs.
         add_gAk0_toRhs(*patch,
                        patch->getPatchGeometry()->getCodimensionBoundaries(1),
                        d_physical_bc_coef_strategy, q_rhs_data);
         add_gAk0_toRhs(*patch,
                        d_cf_boundary->getBoundaries(patch->getGlobalId(), 1),
                        &d_cf_bc_coef, q_rhs_data);
      }

      // Copy rhs from the hierarchy into the HYPRE std::vector
      for (int depth = 0; depth < d_qlen; depth++) {
         copyToHypre(depth, q_rhs_data, d_linear_rhs[depth]);
         // newcopyToHypre(d_linear_rhs[depth], q_rhs_data, depth, box );
      }

   }  // end patch loop

   // Reset invalid ghost data id to help detect erroneous further use.
   d_cf_bc_coef.setGhostDataId(-1, hier::IntVector(tbox::Dimension(NDIM), 0));

   // Finish assembly of the Hypre solution and rhs std::vectors
   for (int depth = 0; depth < d_qlen; depth++) {
      HYPRE_SStructVectorAssemble(d_linear_sol[depth]);
      HYPRE_SStructVectorAssemble(d_linear_rhs[depth]);
   }

   // Dump state of solve data before the solve
   if (d_print_solver_info) {
      for (int depth = 0; depth < d_qlen; depth++) {
         char fname[15];
         snprintf(fname, 15, "sol0.out_q%d", depth);
         HYPRE_SStructVectorPrint(fname, d_linear_sol[depth], 1);
         snprintf(fname, 15, "mat0.out_q%d", depth);
         HYPRE_SStructMatrixPrint(fname, d_matrix[depth], 1);
         snprintf(fname, 15, "rhs.out_q%d", depth);
         HYPRE_SStructVectorPrint(fname, d_linear_rhs[depth], 1);
      }
   }

   t_hypre_solve->start();

   // Solve the system
   switch (d_solver_type) {

      case 0:  // PCG with SysPFMG
         d_relative_residual_norm = 0.;
         d_number_iterations = 0;
         for (int depth = 0; depth < d_qlen; depth++) {

            HYPRE_StructPCGSetTol(d_struct_solver[depth],
                                  d_relative_residual_tol);
            HYPRE_StructPCGSetMaxIter(d_struct_solver[depth], d_max_iterations);
            HYPRE_StructPCGSolve(d_struct_solver[depth], d_struct_A[depth],
                                 d_struct_b[depth], d_struct_x[depth]);
            double relative_residual_norm;
            HYPRE_StructPCGGetFinalRelativeResidualNorm(
                d_struct_solver[depth], &relative_residual_norm);
            int number_iterations;
            HYPRE_StructPCGGetNumIterations(d_struct_solver[depth],
                                            &number_iterations);

            if (d_verbose) {
               tbox::pout << "PFMG preconditioned PCG on level " << d_ln
                          << ": Residual norm = " << relative_residual_norm
                          << "  "
                          << " Iterations = " << number_iterations << std::endl;
            }

            if (relative_residual_norm > d_relative_residual_norm)
               d_relative_residual_norm = relative_residual_norm;
            if (number_iterations > d_number_iterations)
               d_number_iterations = number_iterations;
         }

         break;
      case 1:  // SysPFMG
         d_relative_residual_norm = 0.;
         d_number_iterations = 0;
         for (int depth = 0; depth < d_qlen; depth++) {
            double relative_residual_norm;
            int number_iterations;

            HYPRE_StructPFMGSetTol(d_struct_solver[depth],
                                   d_relative_residual_tol);
            HYPRE_StructPFMGSetMaxIter(d_struct_solver[depth],
                                       d_max_iterations);
            HYPRE_StructPFMGSolve(d_struct_solver[depth], d_struct_A[depth],
                                  d_struct_b[depth], d_struct_x[depth]);
            HYPRE_StructPFMGGetFinalRelativeResidualNorm(
                d_struct_solver[depth], &relative_residual_norm);
            HYPRE_StructPFMGGetNumIterations(d_struct_solver[depth],
                                             &number_iterations);

            if (d_verbose) {
               tbox::pout << "PFMG on level " << d_ln
                          << ": Residual norm = " << relative_residual_norm
                          << "  "
                          << " Iterations = " << number_iterations << std::endl;
            }

            if (relative_residual_norm > d_relative_residual_norm)
               d_relative_residual_norm = relative_residual_norm;
            if (number_iterations > d_number_iterations)
               d_number_iterations = number_iterations;
         }

         break;
      case 2:  // PCG with AMG
         HYPRE_ParCSRPCGSolve(d_parcsr_solver[0], d_parcsr_A[0], d_parcsr_b[0],
                              d_parcsr_x[0]);
         HYPRE_PCGGetNumIterations(d_parcsr_solver[0], &d_number_iterations);
         HYPRE_PCGGetFinalRelativeResidualNorm(d_parcsr_solver[0],
                                               &d_relative_residual_norm);
         break;
      case 3:  // AMG
         HYPRE_BoomerAMGSolve(d_parcsr_solver[0], d_parcsr_A[0], d_parcsr_b[0],
                              d_parcsr_x[0]);
         HYPRE_BoomerAMGGetNumIterations(d_parcsr_solver[0],
                                         &d_number_iterations);
         HYPRE_BoomerAMGGetFinalRelativeResidualNorm(d_parcsr_solver[0],
                                                     &d_relative_residual_norm);
         break;
      default:
         TBOX_ERROR(d_object_name << ": ERROR: Invalid solver id specified.\n");
   }

   t_hypre_solve->stop();

   // Gather the solution std::vector if the  object  type is parcsr

   if (d_hypre_object_type == HYPRE_PARCSR) {
      for (int depth = 0; depth < d_qlen; depth++) {
         HYPRE_SStructVectorGather(d_linear_sol[depth]);
      }
   }

   // Print the solution and other info
   if (d_print_solver_info) {
      // Dump solution std::vector after the solve
      for (int depth = 0; depth < d_qlen; depth++) {
         char fname[20];
         snprintf(fname, 20, "sstruct_%d.out", depth);
         HYPRE_SStructVectorPrint(fname, d_linear_sol[depth], 0);
      }

      tbox::pout << std::endl
                 << "Iterations = " << d_number_iterations << std::endl
                 << "Final Relative Residual Norm = "
                 << d_relative_residual_norm << std::endl;
   }

   /*
    * Pull the solution std::vector out of the HYPRE structures
    */
   for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
        ++ip) {
      const std::shared_ptr<hier::Patch> &patch = *ip;
      std::shared_ptr<pdat::CellData<double> > q_solution_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(q_solution_id)));

      for (int depth = 0; depth < d_qlen; depth++) {
         copyFromHypre(0, d_linear_sol[depth], depth, *q_solution_data);
      }
   }

   t_solve_system->stop();

   return (d_relative_residual_norm <= d_relative_residual_tol);
}


void QuatLevelSolver::adjustBoundaryEntries(
    double gamma, pdat::CellData<double> &sqrt_mobility_data,
    pdat::CellData<double> &diagonal,
    const pdat::SideData<double> &off_diagonal, const hier::Box &patch_box,
    const pdat::ArrayData<double> &acoef_data,
    const pdat::ArrayData<double> &bcoef_data, const hier::Box bccoef_box,
    pdat::OutersideData<double> &Ak0,
    const hier::BoundaryBox &trimmed_boundary_box, const double h[NDIM])
{
   assert(Ak0.getDepth() == d_qlen);

   const hier::Index patch_lo = patch_box.lower();
   const hier::Index patch_up = patch_box.upper();
   const int location_index = trimmed_boundary_box.getLocationIndex();
   const hier::Index &lower = trimmed_boundary_box.getBox().lower();
   const hier::Index &upper = trimmed_boundary_box.getBox().upper();

   const hier::Box &m_gbox = sqrt_mobility_data.getGhostBox();
   const hier::Index &mlower = m_gbox.lower();
   const hier::Index &mupper = m_gbox.upper();

   pdat::ArrayData<double> &Ak0_data =
       Ak0.getArrayData(location_index / 2, location_index % 2);
   const hier::Box &Ak0_box = Ak0_data.getBox();

   for (int depth = 0; depth < d_qlen; depth++) {

#if NDIM == 2
      ADJUST_BDRY2D(
          gamma, sqrt_mobility_data.getPointer(), mlower[0], mupper[0],
          mlower[1], mupper[1], diagonal.getPointer(depth),
          off_diagonal.getPointer(0, depth), off_diagonal.getPointer(1, depth),
          &patch_lo[0], &patch_up[0], &patch_lo[1], &patch_up[1],
          acoef_data.getPointer(), bcoef_data.getPointer(),
          &bccoef_box.lower()[0], &bccoef_box.upper()[0],
          &bccoef_box.lower()[1], &bccoef_box.upper()[1],
          Ak0.getPointer(location_index / 2, location_index % 2, depth),
          &Ak0_box.lower()[0], &Ak0_box.upper()[0], &Ak0_box.lower()[1],
          &Ak0_box.upper()[1], &lower[0], &upper[0], &location_index, h);
#endif
#if NDIM == 3
      ADJUST_BDRY3D(
          gamma, sqrt_mobility_data.getPointer(), mlower[0], mupper[0],
          mlower[1], mupper[1], mlower[2], mupper[2],
          diagonal.getPointer(depth), off_diagonal.getPointer(0, depth),
          off_diagonal.getPointer(1, depth), off_diagonal.getPointer(2, depth),
          &patch_lo[0], &patch_up[0], &patch_lo[1], &patch_up[1], &patch_lo[2],
          &patch_up[2], acoef_data.getPointer(), bcoef_data.getPointer(),
          &bccoef_box.lower()[0], &bccoef_box.upper()[0],
          &bccoef_box.lower()[1], &bccoef_box.upper()[1],
          &bccoef_box.lower()[2], &bccoef_box.upper()[2],
          Ak0.getPointer(location_index / 2, location_index % 2, depth),
          &Ak0_box.lower()[0], &Ak0_box.upper()[0], &Ak0_box.lower()[1],
          &Ak0_box.upper()[1], &Ak0_box.lower()[2], &Ak0_box.upper()[2],
          &lower[0], &upper[0], &location_index, h);
#endif
#ifdef DEBUG_CHECK_ASSERTIONS
      math::ArrayDataNormOpsReal<double> ops;
      double na = ops.maxNorm(Ak0_data, Ak0_box);
      if (na != na) {
         Ak0.print(Ak0_box, std::cerr);
         std::cerr << "depth=" << depth << ", Ak0: Norm after=" << na
                   << std::endl;
         TBOX_ERROR("Nan in Ak0!!!\n");
      }
      if (na > 1.e25) {
         Ak0.print(Ak0_box, std::cerr);
         std::cerr << "depth=" << depth << ", Ak0: Norm =" << na << std::endl;
         TBOX_ERROR("Inf in Ak0!!!\n");
      }
#endif
   }
}


void QuatLevelSolver::freeVariables()
{
   for (int d = 0; d < NDIM; ++d) {
      s_Ak0_var[d].reset();
   }
}
