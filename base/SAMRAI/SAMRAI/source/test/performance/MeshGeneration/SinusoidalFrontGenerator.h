/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   SinusoidalFrontGenerator class declaration
 *
 ************************************************************************/
#ifndef included_SinusoidalFrontGenerator
#define included_SinusoidalFrontGenerator

#include "MeshGenerationStrategy.h"

#include <string>

/*
 * SAMRAI classes
 */
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"

#include <boost/shared_ptr.hpp>

using namespace SAMRAI;

/*!
 * @brief Class to tag a sinusoidal "front" in given domain.
 *
 * Inputs:
 *
 * buffer_distance_0, buffer_distance_1, ...:
 * buffer_distance[ln] is the buffer distance when tagging ON
 * level ln.  We tag the fronts and buffer the tags by this amount.
 * If buffer_distance[ln] is not given, then tag by buffer_distance[ln].
 *
 * shrink_distance_0, shrink_distance_1, ...:
 * shrink_distance[ln] is the shink distance when tagging ON
 * level ln by shrinking the boundaries of level ln.
 * We use this method ONLY when shrink_distance[ln] is not given.
 */
class SinusoidalFrontGenerator:
   public MeshGenerationStrategy
{

public:
   /*!
    * @brief Constructor.
    */
   SinusoidalFrontGenerator(
      /*! Ojbect name */
      const std::string& object_name,
      const tbox::Dimension& dim,
      /*! Input database */
      const boost::shared_ptr<tbox::Database> &database = boost::shared_ptr<tbox::Database>() );

   ~SinusoidalFrontGenerator();

   /*!
    * @brief Set tas on the tag level.
    */
   virtual void setTags(
      bool &exact_tagging,
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      int tag_ln,
      int tag_data_id);

   //@{ @name SAMRAI::mesh::StandardTagAndInitStrategy virtuals

public:

   /*!
    * @brief Set the domain, possibly scaling up the specifications.
    *
    * Take the domain_boxes, xlo and xhi to be the size for the
    * (integer) value of autoscale_base_nprocs.  Scale the problem
    * from there to the number of process running by doubling the
    * dimension starting with the j direction.
    *
    * The number of processes must be a power of 2 times the value
    * of autoscale_base_nprocs.
    */
   void setDomain(
      hier::BoxContainer &domain,
      double xlo[],
      double xhi[],
      int autoscale_base_nprocs,
      const tbox::SAMRAI_MPI &mpi);

   virtual void
   resetHierarchyConfiguration(
      /*! New hierarchy */
      const boost::shared_ptr<hier::PatchHierarchy> &new_hierarchy,
      /*! Coarsest level */ const int coarsest_level,
      /*! Finest level */ const int finest_level);

   //@}

   void
   initializePatchData(
      hier::Patch& patch,
      const double init_data_time,
      const bool initial_time,
      const bool allocate_data)
      {
         NULL_USE(patch);
         NULL_USE(init_data_time);
         NULL_USE(initial_time);
         NULL_USE(allocate_data);
         TBOX_ERROR("Should not be here.");
      }

   bool
   packDerivedDataIntoDoubleBuffer(
      double* buffer,
      const hier::Patch& patch,
      const hier::Box& region,
      const std::string& variable_name,
      int depth_index) const;

public:

#ifdef HAVE_HDF5
   /*!
    * @brief Tell a VisIt plotter which data to write for this class.
    */
   int
   registerVariablesWithPlotter(
      appu::VisItDataWriter& writer);
#endif

   /*!
    * @brief Compute distance and tag data for a patch.
    *
    * This method is not specific to data on the hierarchy,
    * so it is of more general use.  It does not require the
    * hierarchy.
    */
   void computeFrontsData(
      pdat::NodeData<double>* dist_data,
      pdat::CellData<int>* tag_data,
      const std::vector<double> &buffer_distance,
      const double xlo[],
      const double dx[],
      const double time ) const;

private:

   /*!
    * @brief Compute distance and tag data for a patch.
    */
   void
   computePatchData(
      const hier::Patch& patch,
      const double time,
      pdat::NodeData<double>* dist_data,
      pdat::CellData<int>* tag_data) const;

   std::string d_name;

   const tbox::Dimension d_dim;

   /*!
    * @brief PatchHierarchy for use in implementations of some
    * abstract interfaces that do not specify a hierarch.
    */
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

   /*
    * @brief Whether to generate tags by buffering around the fronts
    * or by shrinking the tag level.
    *
    * 'b' means tag around the fronts and buffer the tags.  's' means
    * tag by shrinking the tag level.  The amount of buffer or
    * shinkage is given by d_buffer_shrink_distance.
    */
   std::vector<char> d_buffer_shrink;

   /*!
    * @brief Buffer or shrink distances for generating tags.
    *
    * d_buffer_shrink determines how to generate tags on each level.
    * If buffering on level ln, then d_buffer_shrink_distance[ln]
    * means the buffer distances.  If shrinking on level ln, then it
    * means the shrink distances.
    */
   std::vector<std::vector<double> >  d_buffer_shrink_distance;

#if 0
   /*!
    * @brief Number of cells to tag around cells intersecting the front.
    */
   hier::IntVector d_buffer_cells;

   /*!
    * @brief Buffer size in physical space.
    *
    * d_buffer_space[ln*d_dim + d] is the buffer size for level ln in direction d.
    */
   tbox::Array<double> d_buffer_space;

   /*!
    * @brief Shrink distance, for generating L2 by shrinking.
    */
   tbox::Array<double> d_shrink_distance;
#endif

   boost::shared_ptr<tbox::Timer> t_setup;
   boost::shared_ptr<tbox::Timer> t_node_pos;
   boost::shared_ptr<tbox::Timer> t_distance;
   boost::shared_ptr<tbox::Timer> t_tag_cells;
   boost::shared_ptr<tbox::Timer> t_copy;

};

#endif  // included_ssup_SinusoidalFrontGenerator
