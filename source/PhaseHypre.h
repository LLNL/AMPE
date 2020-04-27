/*************************************************************************
 * Adapted from SAMRAI test suite
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 ************************************************************************/
#ifndef included_PhaseHypre
#define included_PhaseHypre

#include "InterpolationType.h"
#include "CellPoissonHypreSolver.h"
#include "PoissonSpecifications.h"

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"


using namespace SAMRAI;

/*!
 * @brief Class to solve a sample Poisson equation on a SAMR grid.
 *
 * This class demonstrates how use the FAC Poisson solver
 * class to solve Poisson's equation on a single level
 * within a hierarchy.
 *
 * We set up and solve the following problem:
 *
 *   2d: div(grad(u)) = -2 (pi^2) sin(pi x) sin(pi y)
 *
 *   3d: div(grad(u)) = -3 (pi^2) sin(pi x) sin(pi y) sin(pi z)
 *
 * which has the exact solution
 *
 *   2d: u = sin(pi x) sin(pi y)
 *
 *   3d: u = sin(pi x) sin(pi y) sin(pi z)
 *
 * This class inherits and implements virtual functions from
 * - mesh::StandardTagAndInitStrategy to initialize data
 *   on the SAMR grid.
 * - appu::VisDerivedDataStrategy to write out certain data
 *   in a vis file, such as the error of the solution.
 *
 * Inputs:  The only input parameter for this class is
 * "fac_poisson", the input database for the solv::PhaseHypreSolver
 * object.  See the documentation for solv::PhaseHypreSolver
 * for its input parameters.
 */
class PhaseHypre : public mesh::StandardTagAndInitStrategy,
                   public appu::VisDerivedDataStrategy
{

 public:
   /*!
    * @brief Constructor.
    *
    * If you want standard output and logging,
    * pass in valid pointers for those streams.
    *
    * @param object_name Ojbect name
    * @param dim
    * @param fac_solver
    * @param bc_coefs
    */
   PhaseHypre(
       const std::string& object_name, const tbox::Dimension& dim,
       const std::shared_ptr<solv::LocationIndexRobinBcCoefs>& bc_coefs,
       std::shared_ptr<tbox::Database> database, const double epsilon,
       const double omega, const double delta, const double mobility,
       const double gamma);

   virtual ~PhaseHypre();

   //@{ @name mesh::StandardTagAndInitStrategy virtuals

   /*!
    * @brief Allocate and initialize data for a new level
    * in the patch hierarchy.
    *
    * This is where you implement the code for initialize data on
    * the grid.  All the information needed to initialize the grid
    * are in the arguments.
    *
    * @see mesh::StandardTagAndInitStrategy::initializeLevelData()
    */
   virtual void initializeLevelData(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int level_number, const double init_data_time,
       const bool can_be_refined, const bool initial_time,
       const std::shared_ptr<hier::PatchLevel>& old_level,
       const bool allocate_data);

   /*!
    * @brief Reset any internal hierarchy-dependent information.
    */
   virtual void resetHierarchyConfiguration(
       const std::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
       int coarsest_level, int finest_level);

   //@}

   //@{ @name appu::VisDerivedDataStrategy virtuals

   virtual bool packDerivedDataIntoDoubleBuffer(
       double* buffer, const hier::Patch& patch, const hier::Box& region,
       const std::string& variable_name, int depth_id,
       double simulation_time) const;

   //@}

   /*!
    * @brief Solve using HYPRE Poisson solver
    *
    * Set up the linear algebra problem and use a
    * solv::PhaseHypreSolver object to solve it.
    * -# Set initial guess
    * -# Set boundary conditions
    * -# Specify Poisson equation parameters
    * -# Call solver
    */
   int solve(EnergyInterpolationType phase_interp_func_type, double,
             std::string phase_well_func_type);

#ifdef HAVE_HDF5
   /*!
    * @brief Set up external plotter to plot internal
    * data from this class.
    *
    * After calling this function, the external
    * data writer may be used to write the
    * viz file for this object.
    *
    * The internal hierarchy is used and must be
    * established before calling this function.
    * (This is commonly done by building a hierarchy
    * with the mesh::StandardTagAndInitStrategy virtual
    * functions implemented by this class.)
    *
    * @param viz_writer VisIt writer
    */
   int setupPlotter(appu::VisItDataWriter& plotter) const;
#endif

   double compareSolutionWithExact();

   void setCtoZero();

   void setC(const int phi_id, const double gamma,
             const EnergyInterpolationType phi_interp_func_type,
             const double phi_well_scale, const std::string phi_well_func_type);

 private:
   void setCOnPatchPrivate(std::shared_ptr<pdat::CellData<double> > cd_phi,
                           std::shared_ptr<pdat::CellData<double> > cd_m,
                           std::shared_ptr<pdat::CellData<double> > cd_c,
                           const double gamma, const char* phi_interp_func_type,
                           const double phi_well_scale,
                           const char* phi_well_func_type,
                           const hier::Box& pbox);

   std::string d_object_name;

   const tbox::Dimension d_dim;

   std::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   //@{
   /*!
    * @name Major algorithm objects.
    */

   /*!
    * @brief poisson solver.
    */
   CellPoissonHypreSolver d_poisson_solver;

   //@}

   //@{

   /*!
    * @name Private state variables for solution.
    */

   /*!
    * @brief Context owned by this object.
    */
   std::shared_ptr<hier::VariableContext> d_context;

   /*!
    * @brief Descriptor indices of internal data.
    *
    * These are initialized in the constructor and never change.
    */
   int d_comp_soln_id, d_exact_id, d_rhs_id;
   int d_m_id;
   int d_c_id;

   //@}

   /*!
    * PFM parameters
    */
   double d_epsilon;
   double d_omega;
   double d_delta;
   double d_mobility;
   double d_gamma;

   PoissonSpecifications d_poisson_spec;

   bool d_C_is_set;
   bool d_D_is_set;
   bool d_M_is_set;

   std::shared_ptr<solv::RobinBcCoefStrategy> d_physical_bc_coef;
};

#endif  // included_PhaseHypre
