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

// using namespace std;

SinusoidalFrontTagger::SinusoidalFrontTagger(
   const std::string& object_name,
   const tbox::Dimension& dim,
   tbox::Database* database):
   d_name(object_name),
   d_dim(dim),
   d_period(1.0),
   d_init_disp(dim.getValue()),
   d_velocity(dim.getValue()),
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
      d_period =
         database->getDoubleWithDefault("period",
            d_period);
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
   }

   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      d_init_disp[idim] = idim < init_disp.size() ? init_disp[idim] : 0.0;
      d_velocity[idim] = idim < velocity.size() ? velocity[idim] : 0.0;
   }

   const std::string context_name = d_name + std::string(":context");
   d_context = variable_db->getContext(context_name);

   if (database->isInteger("ghost_cell_width")) {
      database->getIntegerArray("ghost_cell_width",
         &d_ghost_cell_width[0], d_dim.getValue());
   }

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

   boost::shared_ptr<hier::PatchHierarchy> hierarchy(base_hierarchy);
   boost::shared_ptr<hier::PatchLevel> old_level(old_base_level);
   if (old_base_level) {
      TBOX_ASSERT(old_level);
   }
   TBOX_ASSERT(hierarchy);

   /*
    * Reference the level object with the given index from the hierarchy.
    */
   boost::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

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
      computeLevelData(hierarchy, ln, d_time /*init_data_time*/,
         d_dist_id, d_tag_id, old_level);
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

   boost::shared_ptr<hier::PatchHierarchy> hierarchy_(base_hierarchy_);
   TBOX_ASSERT(hierarchy_);
   boost::shared_ptr<hier::PatchLevel> level_(hierarchy_->getPatchLevel(ln));
   TBOX_ASSERT(level_);

   hier::PatchLevel& level = *level_;

   for (hier::PatchLevel::iterator pi(level.begin());
        pi != level.end(); ++pi) {
      hier::Patch& patch = **pi;

      boost::shared_ptr<hier::PatchData> tag_data(
         patch.getPatchData(tag_index));
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
            patch.getPatchData(d_tag_id));
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
         dist_data =
            boost::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(patch.getPatchData(dist_id));
      }
      boost::shared_ptr<pdat::CellData<int> > tag_data;
      if (tag_id >= 0) {
         tag_data =
            boost::dynamic_pointer_cast<pdat::CellData<int>,
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

   const hier::Box& pbox = patch.getBox();
   boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());

   const double* xlo = patch_geom->getXLower();

   const double* dx = patch_geom->getDx();

   // Compute the size of buffer to tag around cells crossing front.
   hier::IntVector buffer(d_buffer_cells);
   for (int i = 0; i < d_dim.getValue(); ++i) {
      if (d_buffer_space.size() > ln * d_dim.getValue() + i) {
         int space_based_buffer =
            int(d_buffer_space[ln * d_dim.getValue() + i] / dx[i] + 0.5);
         if (space_based_buffer > buffer(i)) buffer(i) = space_based_buffer;
      }
   }
   // std::cout << "buffer for ln of " << ln << " is " << buffer << std::endl;

   /*
    * We need at least buffer ghost cells to compute
    * the tags, but the data does not have as many ghost cells.
    * So we create temporary patch data with the required "ghost"
    * buffer for computing tag values.  (We could give the real
    * data the required ghost cells, but that may affect the
    * regridding algorithm I'm testing.)
    */
   hier::IntVector required_tmp_buffer(buffer);
   required_tmp_buffer *= ratio;
   pdat::NodeData<double> tmp_dist(pbox, 1, required_tmp_buffer);
   pdat::CellData<int> tmp_tag(pbox, 1, required_tmp_buffer);

   /*
    * Determine what x-node-index contains the sinusoidal front.
    */

   const double wave_number = 2 * 3.141592654 / d_period;

   t_setup->stop();

   t_node_pos->start();

   hier::Box front_box = pbox;
   front_box.grow(required_tmp_buffer);
   front_box.growUpper(hier::IntVector(d_dim, 1));
   // Squash front_box to a single plane.
   front_box.upper(0) = front_box.lower(0);
   const int ifront = front_box.lower(0);

   pdat::ArrayData<int> front_i_(front_box, 1);

   MDA_Access<int, 2, MDA_OrderColMajor<2> > front_i2;
   MDA_Access<int, 3, MDA_OrderColMajor<3> > front_i3;
   if (d_dim == tbox::Dimension(2)) {
      front_i2 = MDA_Access<int, 2, MDA_OrderColMajor<2> >(
            front_i_.getPointer(0),
            &front_i_.getBox().lower()[0],
            &front_i_.getBox().upper()[0]);
      for (int j = front_i2.beg(1); j < front_i2.end(1); ++j) {
         double y = xlo[1] + dx[1] * (j - pbox.lower(1));
         double siny =
            sin(wave_number * (y + d_init_disp[1] - d_velocity[1] * time));
         double fx = d_amplitude * siny + d_init_disp[0] + d_velocity[0] * time;
         front_i2(ifront, j) = int((fx - xlo[0]) / dx[0]) + pbox.lower(0);
         // std::cout << i << '\t' << j << '\t' << y << '\t' << front_i(i,j) << std::endl;
      }
   } else if (d_dim == tbox::Dimension(3)) {
      front_i3 = MDA_Access<int, 3, MDA_OrderColMajor<3> >(
            front_i_.getPointer(0),
            &front_i_.getBox().lower()[0],
            &front_i_.getBox().upper()[0]);
      for (int k = front_i3.beg(2); k < front_i3.end(2); ++k) {
         double z = xlo[2] + dx[2] * (k - pbox.lower(2));
         double sinz =
            sin(wave_number * (z + d_init_disp[2] - d_velocity[2] * time));
         for (int j = front_i3.beg(1); j < front_i3.end(1); ++j) {
            double y = xlo[1] + dx[1] * (j - pbox.lower(1));
            double siny =
               sin(wave_number * (y + d_init_disp[1] + d_velocity[1] * time));
            double fx = d_amplitude * siny * sinz + d_init_disp[0]
               + d_velocity[0] * time;
            front_i3(ifront, j, k) = int((fx - xlo[0]) / dx[0]) + pbox.lower(0);
            // std::cout << i << '\t' << j << '\t' << k << '\t' << y << '\t' << z << '\t' << front_i(i,j,k) << std::endl;
         }
      }
   }

   t_node_pos->stop();

   if (dist_data != NULL) {
      t_distance->start();

      pdat::NodeData<double>::iterator ni(tmp_dist.getGhostBox(), true);
      pdat::NodeData<double>::iterator niend(tmp_dist.getGhostBox(), false);
      for ( ; ni != niend; ++ni) {
         const pdat::NodeIndex& index = *ni;
         if (d_dim == tbox::Dimension(2)) {
            tmp_dist(index) = xlo[0] + (index(0) - pbox.lower(0)) * dx[0]
               - front_i2(ifront, index(1)) * dx[0];
         } else if (d_dim == tbox::Dimension(3)) {
            tmp_dist(index) = xlo[0] + (index(0) - pbox.lower(0)) * dx[0]
               - front_i3(ifront, index(1), index(2)) * dx[0];
         }
      }
      // tmp_dist.print(tmp_dist.getBox(),0,plog);

      t_distance->stop();
   }

   if (tag_data != NULL) {

      t_tag_cells->start();

      tag_data->fill(0);

      const hier::IntVector tag_growth(buffer);

      if (d_dim == tbox::Dimension(2)) {
         MDA_Access<int, 2, MDA_OrderColMajor<2> > tag_aa(
            tag_data->getPointer(0),
            &tag_data->getGhostBox().lower()[0],
            &tag_data->getGhostBox().upper()[0]);
         for (int j = pbox.lower(1); j <= pbox.upper(1); ++j) {
            int mini = front_i2(ifront, j) - buffer(0);
            int maxi = front_i2(ifront, j) + buffer(0);
            if (mini < pbox.lower() (0)) mini = pbox.lower() (0);
            if (maxi > pbox.upper() (0)) maxi = pbox.upper() (0);
            for (int i = mini; i <= maxi; ++i) {
               tag_aa(i, j) = 1;
            }
         }
      } else if (d_dim == tbox::Dimension(3)) {
         MDA_Access<int, 3, MDA_OrderColMajor<3> > tag_aa(
            tag_data->getPointer(0),
            &tag_data->getGhostBox().lower()[0],
            &tag_data->getGhostBox().upper()[0]);
         for (int k = pbox.lower(2); k <= pbox.upper(2); ++k) {
            for (int j = pbox.lower(1); j <= pbox.upper(1); ++j) {
               int mini = front_i3(ifront, j, k) - buffer(0);
               int maxi = front_i3(ifront, j, k) + buffer(0);
               if (mini < pbox.lower() (0)) mini = pbox.lower() (0);
               if (maxi > pbox.upper() (0)) maxi = pbox.upper() (0);
               for (int i = mini; i <= maxi; ++i) {
                  tag_aa(i, j, k) = 1;
               }
            }
         }
      }

      t_tag_cells->stop();

   }

   t_copy->start();

   /*
    * Copy computed data to output.  Recall that the convention is
    * to send in a NULL pointer to indicate that data is not wanted.
    */
   if (dist_data != NULL) {
      dist_data->copy(tmp_dist);
   }

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
