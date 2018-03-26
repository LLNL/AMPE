/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   SinusoidalFrontTagger class implementation
 *
 ************************************************************************/
#include "SinusoidalFrontTagger.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/MDA_Access.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iomanip>

using namespace SAMRAI;



SinusoidalFrontTagger::SinusoidalFrontTagger(
   const std::string& object_name,
   const tbox::Dimension& dim,
   tbox::Database* database):
   d_name(object_name),
   d_dim(dim),
   d_amplitude(0.2),
   d_ghost_cell_width(dim, 0),
   d_buffer_cells(dim, 1),
   d_allocate_data(true),
   d_time(0.5)
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(variable_db != NULL);
#endif

   tbox::Array<double> init_disp;
   tbox::Array<double> velocity;
   tbox::Array<double> period;

   if (database != NULL) {
      d_allocate_data =
         database->getBoolWithDefault("allocate_data",
            d_allocate_data);
      if (database->isInteger("buffer_cells")) {
         database->getIntegerArray("buffer_cells",
            &d_buffer_cells[0], d_dim.getValue());
      }
      for (int ln = 0; true; ++ln) {
         std::string name("buffer_space_");
         name = name + tbox::Utilities::intToString(ln);
         if (database->isDouble(name)) {
            d_buffer_space.resizeArray(d_dim.getValue() * (ln + 1));
            database->getDoubleArray(name, &d_buffer_space[d_dim.getValue() * ln], d_dim.getValue());
         } else {
            break;
         }
      }
      if (database->isDouble("period")) {
         period =
            database->getDoubleArray("period");
      }
      if (database->isDouble("init_disp")) {
         init_disp =
            database->getDoubleArray("init_disp");
      }
      if (database->isDouble("velocity")) {
         velocity =
            database->getDoubleArray("velocity");
      }
      d_amplitude =
         database->getDoubleWithDefault("amplitude",
            d_amplitude);
      d_time =
         database->getDoubleWithDefault("time",
            d_time);

      if (database->isInteger("ghost_cell_width")) {
         database->getIntegerArray("ghost_cell_width",
            &d_ghost_cell_width[0], d_dim.getValue());
      }
   }

   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      d_init_disp[idim] = idim < init_disp.size() ? init_disp[idim] : 0.0;
      d_velocity[idim] = idim < velocity.size() ? velocity[idim] : 0.0;
      d_period[idim] = idim < period.size() ? period[idim] : 1.0e20;
   }

   const std::string context_name = d_name + std::string(":context");
   d_context = variable_db->getContext(context_name);

   boost::shared_ptr<hier::Variable> dist_var(
      new pdat::NodeVariable<double>(dim, d_name + ":dist"));
   d_dist_id = variable_db->registerVariableAndContext(dist_var,
         d_context,
         d_ghost_cell_width);

   boost::shared_ptr<hier::Variable> tag_var(
      new pdat::CellVariable<int>(dim, d_name + ":tag"));
   d_tag_id = variable_db->registerVariableAndContext(tag_var,
         d_context,
         d_ghost_cell_width);

   t_setup = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::setup");
   t_node_pos = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::node_pos");
   t_distance = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::distance");
   t_tag_cells = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::tag_cells");
   t_copy = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontTagger::copy");
}



SinusoidalFrontTagger::~SinusoidalFrontTagger()
{
}



void SinusoidalFrontTagger::initializeLevelData(
   /*! Hierarchy to initialize */
   const boost::shared_ptr<hier::PatchHierarchy>& base_hierarchy,
   /*! Level to initialize */
   const int ln,
   const double init_data_time,
   const bool can_be_refined,
   /*! Whether level is being introduced for the first time */
   const bool initial_time,
   /*! Level to copy data from */
   const boost::shared_ptr<hier::PatchLevel>& old_base_level,
   const bool allocate_data)
{
   NULL_USE(can_be_refined);
   NULL_USE(old_base_level);

   TBOX_ASSERT(base_hierarchy);

   /*
    * Reference the level object with the given index from the hierarchy.
    */
   boost::shared_ptr<hier::PatchLevel> level(
      base_hierarchy->getPatchLevel(ln));

   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      hier::Patch& patch = **pi;
      initializePatchData(patch,
         init_data_time,
         initial_time,
         allocate_data);
   }

#if 0
   if (d_allocate_data) {
      /*
       * If instructed, allocate all patch data on the level.
       * Allocate only persistent data.  Scratch data will
       * generally be allocated and deallocated as needed.
       */
      if (allocate_data) {
         level->allocatePatchData(d_dist_id);
         level->allocatePatchData(d_tag_id);
      }
      computeLevelData(base_hierarchy, ln, d_time /*init_data_time*/,
         d_dist_id, d_tag_id, old_base_level);
   }
#endif
}



void SinusoidalFrontTagger::initializePatchData(
   hier::Patch& patch,
   const double init_data_time,
   const bool initial_time,
   const bool allocate_data)
{
   NULL_USE(initial_time);

   if (d_allocate_data) {
      /*
       * If instructed, allocate all patch data on the level.
       * Allocate only persistent data.  Scratch data will
       * generally be allocated and deallocated as needed.
       */
      if (allocate_data) {
         if (!patch.checkAllocated(d_dist_id)) {
            patch.allocatePatchData(d_dist_id);
         }
         if (!patch.checkAllocated(d_tag_id)) {
            patch.allocatePatchData(d_tag_id);
         }
         boost::shared_ptr<pdat::NodeData<double> > dist_data(
            patch.getPatchData(d_dist_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::CellData<int> > tag_data(
            patch.getPatchData(d_tag_id),
            boost::detail::dynamic_cast_tag());
         TBOX_ASSERT(dist_data);
         TBOX_ASSERT(tag_data);
         computePatchData(patch, init_data_time,
            dist_data.get(), tag_data.get());
      }
   }
}



void SinusoidalFrontTagger::resetHierarchyConfiguration(
   /*! New hierarchy */ const boost::shared_ptr<hier::PatchHierarchy>& new_hierarchy,
   /*! Coarsest level */ int coarsest_level,
   /*! Finest level */ int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
   d_hierarchy = new_hierarchy;
   TBOX_ASSERT(d_hierarchy);
}



void SinusoidalFrontTagger::applyGradientDetector(
   const boost::shared_ptr<hier::PatchHierarchy>& base_hierarchy_,
   const int ln,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation)
{
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation);
   TBOX_ASSERT(base_hierarchy_);
   boost::shared_ptr<hier::PatchLevel> level_(
      base_hierarchy_->getPatchLevel(ln));
   TBOX_ASSERT(level_);

   hier::PatchLevel& level = *level_;

   for (hier::PatchLevel::iterator pi(level.begin());
        pi != level.end(); ++pi) {
      hier::Patch& patch = **pi;

      boost::shared_ptr<hier::PatchData> tag_data(
         patch.getPatchData(tag_index),
         boost::detail::dynamic_cast_tag());
      if (!tag_data) {
         TBOX_ERROR("Data index " << tag_index
                                  << " does not exist for patch.\n");
      }
      boost::shared_ptr<pdat::CellData<int> > tag_cell_data_(
         tag_data,
         boost::detail::dynamic_cast_tag());
      if (!tag_cell_data_) {
         TBOX_ERROR("Data index " << tag_index
                                  << " is not cell int data.\n");
      }

      if (d_allocate_data) {
         // Use internally stored data.
         boost::shared_ptr<hier::PatchData> saved_tag_data(
            patch.getPatchData(d_tag_id),
            boost::detail::dynamic_cast_tag());
         tag_cell_data_->copy(*saved_tag_data);
      } else {
         // Compute tag data for patch.
         computePatchData(patch,
            error_data_time,
            NULL,
            tag_cell_data_.get());
      }

   }
}



/*
 * Deallocate patch data allocated by this class.
 */

void SinusoidalFrontTagger::deallocatePatchData(
   hier::PatchHierarchy& hierarchy)
{
   int ln;
   for (ln = 0; ln < hierarchy.getNumberOfLevels(); ++ln) {
      boost::shared_ptr<hier::PatchLevel> level(hierarchy.getPatchLevel(ln));
      deallocatePatchData(*level);
   }
}



/*
 * Deallocate patch data allocated by this class.
 */

void SinusoidalFrontTagger::deallocatePatchData(
   hier::PatchLevel& level)
{
   level.deallocatePatchData(d_dist_id);
   level.deallocatePatchData(d_tag_id);
}



/*
 * Deallocate patch data allocated by this class.
 */
void SinusoidalFrontTagger::computeHierarchyData(
   hier::PatchHierarchy& hierarchy,
   double time)
{
   d_time = time;
   if (!d_allocate_data) return;

   for (int ln = 0; ln < hierarchy.getNumberOfLevels(); ++ln) {
      computeLevelData(hierarchy, ln, time, d_dist_id, d_tag_id);
   }
}



/*
 * Compute the solution data for a level.
 * Can copy data from old level (if any) to support
 * initializeLevelData().
 */

void SinusoidalFrontTagger::computeLevelData(
   const hier::PatchHierarchy& hierarchy,
   const int ln,
   const double time,
   const int dist_id,
   const int tag_id,
   const boost::shared_ptr<hier::PatchLevel>& old_level) const
{
   NULL_USE(old_level);

   const boost::shared_ptr<hier::PatchLevel> level(
      hierarchy.getPatchLevel(ln));

   /*
    * Initialize data in all patches in the level.
    */
   for (hier::PatchLevel::iterator pi(level->begin());
        pi != level->end(); ++pi) {
      hier::Patch& patch = **pi;
      boost::shared_ptr<pdat::NodeData<double> > dist_data;
      if (dist_id >= 0) {
         dist_data = boost::dynamic_pointer_cast<pdat::NodeData<double>,
                                                 hier::PatchData>(patch.getPatchData(dist_id));
      }
      boost::shared_ptr<pdat::CellData<int> > tag_data;
      if (tag_id >= 0) {
         tag_data = boost::dynamic_pointer_cast<pdat::CellData<int>,
                                                hier::PatchData>(patch.getPatchData(tag_id));
      }
      computePatchData(patch, time,
         dist_data.get(),
         tag_data.get());
   }
}



/*
 * Compute the solution data for a patch.
 */

void SinusoidalFrontTagger::computePatchData(
   const hier::Patch& patch,
   const double time,
   pdat::NodeData<double>* dist_data,
   pdat::CellData<int>* tag_data) const
{

   t_setup->start();

   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT(patch.inHierarchy());

   const int ln = patch.getPatchLevelNumber();
   const boost::shared_ptr<hier::PatchLevel> level(
      d_hierarchy->getPatchLevel(ln));
   const hier::IntVector& ratio(level->getRatioToLevelZero());

   boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());

   const double* xlo = patch_geom->getXLower();

   const double* dx = patch_geom->getDx();

   /*
    * Compute the size of buffer to tag around cells crossing front.
    * They should be at least d_buffer_cells in reference index space
    * and at least getBufferSpace(ln) in physical space.
    */
   hier::IntVector buffer(d_buffer_cells);
   for (int i = 0; i < d_dim.getValue(); ++i) {
      const double *buffer_space = getBufferSpace(ln);
      if ( buffer_space != NULL ) {
         int space_based_buffer =
            int(d_buffer_space[ln * d_dim.getValue() + i] / dx[i] + 0.5);
         if (space_based_buffer > buffer(i)) {
            buffer(i) = space_based_buffer;
         }
      }
   }
   buffer *= ratio;

   computeFrontsData( dist_data, tag_data, buffer, xlo, dx, time );
}



/*
 * Compute the various data due to the fronts.
 */
void SinusoidalFrontTagger::computeFrontsData(
   pdat::NodeData<double>* dist_data,
   pdat::CellData<int>* tag_data,
   const hier::IntVector &tag_buffer,
   const double xlo[],
   const double dx[],
   const double time ) const
{
   const tbox::Dimension &dim(tag_buffer.getDim());

   t_setup->start();

   if ( dist_data != NULL && tag_data != NULL ) {
      TBOX_ASSERT( dist_data->getBox().isSpatiallyEqual(tag_data->getBox()) );
   }

   const hier::Box& pbox = tag_data->getBox();

   /*
    * We need at least tag_buffer ghost cells to compute
    * the tags, but the data may not have as many ghost cells.
    * So we create temporary patch data with the required
    * tag_buffer for computing tag values.  (We could give the real
    * data the required ghost cells, but that may affect the
    * regridding algorithm I'm testing.)
    */
   pdat::CellData<int> tmp_tag(pbox, 1, tag_buffer);

   /*
    * Determine what x-cell-index contains the sinusoidal front.
    */

   double wave_number[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      wave_number[idim] = 2 * 3.141592654 / d_period[idim];
   }

   t_setup->stop();


   /*
    * Initialize node x-distances from front.
    */

   t_node_pos->start();
   hier::Box front_box = pbox;
   front_box.grow(tag_buffer);
   front_box.growUpper(hier::IntVector(d_dim, 1));
   // Squash front_box to a single plane.
   front_box.upper(0) = front_box.lower(0) = pbox.lower(0);
   pdat::ArrayData<double> front_x(front_box, 1);
   pdat::ArrayData<int>::iterator aiend(front_x.getBox(), false);
   for ( pdat::ArrayData<int>::iterator ai(front_x.getBox(), true);
         ai != aiend; ++ai ) {

      const hier::Index &index = *ai;
      double y=0.0, siny=0.0, z=0.0, sinz=1.0;

      y = xlo[1] + dx[1]*(index(1) - pbox.lower()[1]);
      siny = sin(wave_number[1] * (y + d_init_disp[1] - d_velocity[1] * time));

      if ( dim.getValue() > 2 ) {
         z = xlo[2] + dx[2]*(index(2) - front_box.lower()[2]);
         sinz = sin(wave_number[2] * (z + d_init_disp[2] - d_velocity[2] * time));
      }

      front_x(index,0) = d_amplitude * siny * sinz + d_init_disp[0];
      // tbox::plog << "front_x" << index << " = " << front_x(index,0) << std::endl;
   }
   t_node_pos->stop();

   /*
    * Initialize tmp_tag to zero then tag specific cells.
    */
   tmp_tag.fill(0);
   hier::BlockId blk0(0);
   pdat::CellData<int>::iterator ciend(tmp_tag.getGhostBox(), false);
   for ( pdat::CellData<int>::iterator ci(tmp_tag.getGhostBox(), true);
         ci != ciend; ++ci ) {

      const pdat::CellIndex &cell_index = *ci;
      const hier::Box cell_box(cell_index, cell_index, blk0);

#if 0
      int node_orientation = 0;
#endif
      double min_distance_to_front =  tbox::MathUtilities<double>::getMax();
      double max_distance_to_front = -tbox::MathUtilities<double>::getMax();
      // tbox::plog << "initial distances to front: " << min_distance_to_front << " .. " << max_distance_to_front << std::endl;
      pdat::NodeIterator niend(cell_box, false);
      for ( pdat::NodeIterator ni(cell_box, true); ni != niend; ++ni ) {

         const pdat::NodeIndex &node_index = *ni;
         hier::Index front_index = node_index;
         front_index(0) = pbox.lower(0);

         double node_x = xlo[0] + dx[0]*( node_index(0) - pbox.lower()(0) );

         double distance_to_front = node_x - front_x(front_index,0);
         min_distance_to_front = tbox::MathUtilities<double>::Min(min_distance_to_front, distance_to_front);
         max_distance_to_front = tbox::MathUtilities<double>::Max(max_distance_to_front, distance_to_front);
         // tbox::plog << "cell_index = " << cell_index << "   node_index = " << node_index << "   node_x = " << node_x << "   front_index = " << front_index << "   front_x = " << front_x(front_index,0) << "   distance to front(" << node_x << ") = " << distance_to_front << std::endl;

#if 0
         while ( distance_to_front > 0.5*d_period[0] ) distance_to_front -= d_period[0];
         while ( distance_to_front < -0.5*d_period[0] ) distance_to_front += d_period[0];
         // tbox::plog << "distance_to_front adjusted = " << distance_to_front << std::endl;
         int distance_sign = distance_to_front < 0 ? -1 : 1;
         if ( node_orientation == 0 ) {
            node_orientation = distance_sign;
         }
         else {
            if ( distance_sign != node_orientation ) {
               // This cell has nodes on both sides of the front.  Tag it and the tag_buffer around it.
               hier::Box cell_and_buffer(cell_index, cell_index);
               cell_and_buffer.grow(tag_buffer);
               tmp_tag.fill(1,cell_and_buffer);
            }
         }
#endif
      }
      // tbox::plog << "distances to front: " << min_distance_to_front << " .. " << max_distance_to_front << std::endl;
      while ( min_distance_to_front < -0.5*d_period[0] ) {
         min_distance_to_front += d_period[0];
         max_distance_to_front += d_period[0];
      }
      while ( max_distance_to_front >  0.5*d_period[0] ) {
         min_distance_to_front -= d_period[0];
         max_distance_to_front -= d_period[0];
      }
      // tbox::plog << "shifted ..........: " << min_distance_to_front << " .. " << max_distance_to_front << std::endl;
      if ( min_distance_to_front < 0 && max_distance_to_front > 0 ) {
         // This cell has nodes on both sides of the front.  Tag it and the tag_buffer around it.
         hier::Box cell_and_buffer(cell_index, cell_index, blk0);
         cell_and_buffer.grow(tag_buffer);
         tmp_tag.fill(1,cell_and_buffer);
      }

   }


   /*
    * Initialize distance data.
    */
   if (dist_data != NULL) {
      t_distance->start();

      pdat::NodeData<double> &dist_to_front(*dist_data);
      pdat::NodeData<double>::iterator ni(dist_to_front.getGhostBox(), true);
      pdat::NodeData<double>::iterator niend(dist_to_front.getGhostBox(), false);
      for ( ; ni != niend; ++ni) {
         const pdat::NodeIndex& index = *ni;
         pdat::NodeIndex front_index(index);
         front_index(0) = 0;
         dist_to_front(index) = xlo[0] + (index(0) - pbox.lower(0)) * dx[0]
            - front_x(front_index,0);
      }
      // dist_to_front.print(dist_to_front.getBox(),0,plog);

      t_distance->stop();
   }

   t_copy->start();
   tag_data->copy(tmp_tag);
   t_copy->stop();
}



#ifdef HAVE_HDF5
int SinusoidalFrontTagger::registerVariablesWithPlotter(
   appu::VisItDataWriter& writer)
{
   /*
    * Register variables with plotter.
    */
   if (d_allocate_data) {
      writer.registerPlotQuantity("Distance to front", "SCALAR", d_dist_id);
      writer.registerPlotQuantity("Tag value", "SCALAR", d_tag_id);
   } else {
      writer.registerDerivedPlotQuantity("Distance to front", "SCALAR", this,
         // hier::IntVector(0),
         1.0,
         "NODE");
      writer.registerDerivedPlotQuantity("Tag value", "SCALAR", this);
   }
   return 0;
}
#endif



bool SinusoidalFrontTagger::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_index) const
{
   NULL_USE(region);
   NULL_USE(depth_index);

   TBOX_ASSERT(d_allocate_data == false);
   if (variable_name == "Distance to front") {
      pdat::NodeData<double> dist_data(patch.getBox(), 1, hier::IntVector(d_dim,
                                          0));
      computePatchData(patch, d_time, &dist_data, NULL);
      pdat::NodeData<double>::iterator ciend(patch.getBox(), false);
      for (pdat::NodeData<double>::iterator ci(patch.getBox(), true);
           ci != ciend; ++ci) {
         *(buffer++) = dist_data(*ci);
      }
   } else if (variable_name == "Tag value") {
      pdat::CellData<int> tag_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      computePatchData(patch, d_time, NULL, &tag_data);
      pdat::CellData<double>::iterator ciend(patch.getBox(), false);
      for (pdat::CellData<double>::iterator ci(patch.getBox(), true);
           ci != ciend; ++ci) {
         *(buffer++) = tag_data(*ci);
      }
   } else {
      TBOX_ERROR("Unrecognized name " << variable_name);
   }
   return true;
}



void SinusoidalFrontTagger::setTime(
   double time)
{
   d_time = time;
}



const double *SinusoidalFrontTagger::getBufferSpace( int ln ) const
{
   if ( d_buffer_space.size() > ln*d_dim.getValue() ) {
      return &d_buffer_space[ln*d_dim.getValue()];
   }
   return NULL;
}
