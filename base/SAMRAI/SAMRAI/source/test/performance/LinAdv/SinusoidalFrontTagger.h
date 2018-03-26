/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   SinusoidalFrontTagger class declaration
 *
 ************************************************************************/
#ifndef included_SinusoidalFrontTagger
#define included_SinusoidalFrontTagger

#include <string>
#include <boost/shared_ptr.hpp>
#include "SAMRAI/tbox/Database.h"

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/Timer.h"

using namespace SAMRAI;

/*!
 * @brief Class to tag a sinusoidal "front" in given domain.
 */
class SinusoidalFrontTagger:
   public mesh::StandardTagAndInitStrategy,
   public appu::VisDerivedDataStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   SinusoidalFrontTagger(
      /*! Ojbect name */
      const std::string& object_name,
      const tbox::Dimension& dim,
      /*! Input database */
      tbox::Database* database = NULL);

   ~SinusoidalFrontTagger();

   //@{ @name SAMRAI::mesh::StandardTagAndInitStrategy virtuals

public:
   /*!
    * @brief Allocate and initialize data for a new level
    * in the patch hierarchy.
    *
    * This is where you implement the code for initialize data on the
    * grid.  Nevermind when it is called or where in the program that
    * happens.  All the information you need to initialize the grid
    * are in the arguments.
    *
    * @see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData()
    */
   virtual void
   initializeLevelData(
      /*! Hierarchy to initialize */
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      /*! Level to initialize */
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      /*! Whether level is being introduced for the first time */
      const bool initial_time,
      /*! Level to copy data from */
      const boost::shared_ptr<hier::PatchLevel>& old_level =
         boost::shared_ptr<hier::PatchLevel>(),
      /*! Whether data on new patch needs to be allocated */
      const bool allocate_data = true);

   virtual void
   resetHierarchyConfiguration(
      /*! New hierarchy */
      const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
      /*! Coarsest level */ int coarsest_level,
      /*! Finest level */ int finest_level);

   virtual void
   applyGradientDetector(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation);

   //@}

   void
   initializePatchData(
      hier::Patch& patch,
      const double init_data_time,
      const bool initial_time,
      const bool allocate_data);

   bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_index) const;

   void
   setTime(
      double time);

public:
   /*!
    * @brief Deallocate internally managed patch data on level.
    */
   void
   deallocatePatchData(
      hier::PatchLevel& level);

   /*!
    * @brief Deallocate internally managed patch data on hierarchy.
    */
   void
   deallocatePatchData(
      hier::PatchHierarchy& hierarchy);

#ifdef HAVE_HDF5
   /*!
    * @brief Tell a VisIt plotter which data to write for this class.
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& writer);
#endif

   /*
    * Compute patch data allocated by this class, on a hierarchy.
    */
   void
   computeHierarchyData(
      hier::PatchHierarchy& hierarchy,
      double time);

   /*!
    * @brief Compute distance and tag data for a level.
    */
   void
   computeLevelData(
      const hier::PatchHierarchy& hierarchy,
      const int ln,
      const double time,
      const int dist_id,
      const int tag_id,
      const boost::shared_ptr<hier::PatchLevel>& old_level =
         boost::shared_ptr<hier::PatchLevel>()) const;

   /*!
    * @brief Compute distance and tag data for a patch.
    */
   void
   computePatchData(
      const hier::Patch& patch,
      const double time,
      pdat::NodeData<double>* dist_data,
      pdat::CellData<int>* tag_data) const;

private:
   std::string d_name;

   const tbox::Dimension d_dim;

   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Period of sinusoid.
    */
   double d_period[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*!
    * @brief Initial displacement.
    */
   double d_init_disp[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*!
    * @brief Front velocity.
    */
   double d_velocity[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*!
    * @brief Amplitude of sinusoid.
    */
   double d_amplitude;

   /*!
    * @brief ghost cell width of internal data.
    *
    * Optional.  Meant to influence gridding parameters.  Defaults to zero.
    */
   hier::IntVector d_ghost_cell_width;

   /*!
    * @brief Number of cells to tag around cells intersecting the front.
    */
   hier::IntVector d_buffer_cells;

   tbox::Array<double> d_buffer_space;

   boost::shared_ptr<hier::VariableContext> d_context;

   /*!
    * @brief Distance from the front in the x direction.
    */
   int d_dist_id;
   /*!
    * @brief Value of tag based on distance from front.
    */
   int d_tag_id;

   /*!
    * @brief Whether to allocate data on the mesh.
    */
   bool d_allocate_data;

   /*!
    * @brief Front time.
    */
   double d_time;

   boost::shared_ptr<tbox::Timer> t_setup;
   boost::shared_ptr<tbox::Timer> t_node_pos;
   boost::shared_ptr<tbox::Timer> t_distance;
   boost::shared_ptr<tbox::Timer> t_tag_cells;
   boost::shared_ptr<tbox::Timer> t_copy;

};

#endif  // included_ssup_SinusoidalFrontTagger
