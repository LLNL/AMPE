/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   SinusoidalFrontGenerator class implementation
 *
 ************************************************************************/
#include "SinusoidalFrontGenerator.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iomanip>

using namespace SAMRAI;



SinusoidalFrontGenerator::SinusoidalFrontGenerator(
   const std::string& object_name,
   const tbox::Dimension& dim,
   const boost::shared_ptr<tbox::Database> &database):
   d_name(object_name),
   d_dim(dim),
   d_hierarchy(),
   d_amplitude(0.2)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
   TBOX_ASSERT(variable_db != NULL);
#endif

   tbox::Array<double> init_disp;
   tbox::Array<double> velocity;
   tbox::Array<double> period;


   // Parameters set by database, with defaults.
   boost::shared_ptr<tbox::Database> sft_db; // SinusoidalFrontGenerator database.

   if (database) {

      sft_db = database->getDatabaseWithDefault("SinusoidalFrontGenerator", sft_db);

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

      /*
       * Input parameters to determine whether to tag by buffering
       * fronts or shrinking level, and by how much.
       */
      const std::string sname( "shrink_distance_" );
      const std::string bname( "buffer_distance_" );
      std::string aname;
      for ( int ln=0; ; ++ln ) {
         const std::string lnstr( tbox::Utilities::intToString(ln) );

         // Look for buffer input first, then shrink input.
         const std::string bnameln = bname + lnstr;
         const std::string snameln = sname + lnstr;

         tbox::Array<double> tmpa;

         if ( database->isDouble(bnameln) ) {
            tmpa = database->getDoubleArray(bnameln);
            if ( tmpa.getSize() != dim.getValue() ) {
               TBOX_ERROR(bnameln << " input parameter must have " << dim << " values");
            }
            d_buffer_shrink.push_back('b');
         }

         if ( database->isDouble(snameln) ) {
            if ( !tmpa.empty() ) {
               TBOX_ERROR("Cannot specify both " << bnameln << " and " << snameln);
            }
            tmpa = database->getDoubleArray(snameln);
            if ( tmpa.getSize() != dim.getValue() ) {
               TBOX_ERROR(snameln << " input parameter must have " << dim << " values");
            }
            d_buffer_shrink.push_back('s');
         }

         if ( !tmpa.empty() ) {
            d_buffer_shrink_distance.resize(d_buffer_shrink_distance.size() + 1);
            d_buffer_shrink_distance.back().insert( d_buffer_shrink_distance.back().end(),
                                                    tmpa.getPointer(),
                                                    tmpa.getPointer()+tmpa.size() );
         }
         else {
            break;
         }

      }
   }

   if ( d_buffer_shrink.empty() ) {
      TBOX_ERROR("SinusoidalFrontGenerator: You must specify either\n"
                 << "buffer_distance_# or shrink_distance_# for each level\n"
                 << "you plan to use, except the finest level.");
   }

   for (int idim = 0; idim < d_dim.getValue(); ++idim) {
      d_init_disp[idim] = idim < init_disp.size() ? init_disp[idim] : 0.0;
      d_velocity[idim] = idim < velocity.size() ? velocity[idim] : 0.0;
      d_period[idim] = idim < period.size() ? period[idim] : 1.0e20;
   }

   t_setup = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::setup");
   t_node_pos = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::node_pos");
   t_distance = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::distance");
   t_tag_cells = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::tag_cells");
   t_copy = tbox::TimerManager::getManager()->
      getTimer("apps::SinusoidalFrontGenerator::copy");
}



SinusoidalFrontGenerator::~SinusoidalFrontGenerator()
{
}




/*
 ***********************************************************************
 ***********************************************************************
 */
void SinusoidalFrontGenerator::setTags(
   bool &exact_tagging,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int tag_ln,
   int tag_data_id )
{
   if ( d_buffer_shrink[tag_ln] == 's' ) {
      setTagsByShrinkingLevel(
         hierarchy,
         tag_ln,
         tag_data_id,
         hier::IntVector::getZero(d_dim),
         &d_buffer_shrink_distance[1][0]);
      exact_tagging = true;
      return;
   }

   const boost::shared_ptr<hier::PatchLevel> &tag_level(
      hierarchy->getPatchLevel(tag_ln));

   resetHierarchyConfiguration(hierarchy, 0, 1);

   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {

      boost::shared_ptr<hier::Patch> patch = *pi;

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         patch->getPatchGeometry(),
         boost::detail::dynamic_cast_tag());

      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_data_id),
         boost::detail::dynamic_cast_tag());

      computeFrontsData(
         NULL /* distance data */,
         tag_data.get(),
         d_buffer_shrink_distance[tag_ln],
         patch_geom->getXLower(),
         patch_geom->getDx(),
         0.0 );

   }

   exact_tagging = false;

   return;
}



void SinusoidalFrontGenerator::setDomain(
   hier::BoxContainer &domain,
   double xlo[],
   double xhi[],
   int autoscale_base_nprocs,
   const tbox::SAMRAI_MPI &mpi)
{
   TBOX_ASSERT( autoscale_base_nprocs <= mpi.getSize() );
   TBOX_ASSERT( !domain.isEmpty() );

   hier::BoxContainer::const_iterator ii = domain.begin();
   ii->getDim();
   const tbox::Dimension &dim = domain.begin()->getDim();

   int doubling_dir = 1;
   while (autoscale_base_nprocs < mpi.getSize()) {
      for ( hier::BoxContainer::iterator bi=domain.begin();
            bi!=domain.end(); ++bi ) {
         hier::Box &input_box = *bi;
         input_box.upper()(doubling_dir) += input_box.numberCells(doubling_dir);
      }
      xhi[doubling_dir] += xhi[doubling_dir] - xlo[doubling_dir];
      doubling_dir = (doubling_dir + 1)%dim.getValue();
      autoscale_base_nprocs *= 2;
      tbox::plog << "autoscale_base_nprocs = " << autoscale_base_nprocs << std::endl
                 << domain.format("IB: ", 2) << std::endl;
   }

   if ( autoscale_base_nprocs != mpi.getSize() ) {
      TBOX_ERROR("If autoscale_base_nprocs (" << autoscale_base_nprocs << ") is given,\n"
                 <<"number of processes (" << mpi.getSize() << ") must be\n"
                 <<"a power-of-2 times the value of autoscale_base_nprocs.");
   }

}



void SinusoidalFrontGenerator::resetHierarchyConfiguration(
   /*! New hierarchy */ const boost::shared_ptr<hier::PatchHierarchy> &new_hierarchy,
   /*! Coarsest level */ const int coarsest_level,
   /*! Finest level */ const int finest_level)
{
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);
   TBOX_ASSERT(new_hierarchy->getDim() == d_dim);
   d_hierarchy = new_hierarchy;
   TBOX_ASSERT(d_hierarchy);
}



/*
 * Compute the solution data for a patch.
 */

void SinusoidalFrontGenerator::computePatchData(
   const hier::Patch& patch,
   const double time,
   pdat::NodeData<double>* dist_data,
   pdat::CellData<int>* tag_data) const
{

   t_setup->start();

   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT(patch.inHierarchy());

   const int ln = patch.getPatchLevelNumber();
   const boost::shared_ptr<hier::PatchLevel> level =
      d_hierarchy->getPatchLevel(ln);

   boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());

   const double* xlo = patch_geom->getXLower();

   const double* dx = patch_geom->getDx();

   computeFrontsData( dist_data, tag_data, d_buffer_shrink_distance[patch.getPatchLevelNumber()], xlo, dx, time );
}



/*
 * Compute the various data due to the fronts.
 */
void SinusoidalFrontGenerator::computeFrontsData(
   pdat::NodeData<double>* dist_data,
   pdat::CellData<int>* tag_data,
   const std::vector<double> &buffer_distance,
   const double xlo[],
   const double dx[],
   const double time ) const
{
   const tbox::Dimension &dim(tag_data->getDim());

   t_setup->start();

   if ( dist_data != NULL && tag_data != NULL ) {
      TBOX_ASSERT( dist_data->getBox().isSpatiallyEqual(tag_data->getBox()) );
   }

   // Compute the buffer in terms of cells.
   hier::IntVector buffer_cells(dim);
   for ( int i=0; i<dim.getValue(); ++i ) {
      buffer_cells(i) = static_cast<int>(0.5 + buffer_distance[i]/dx[i]);
   }

   const hier::Box& pbox = tag_data->getBox();

   /*
    * We need at least buffer_cells ghost cells to compute
    * the tags, but the data may not have as many ghost cells.
    * So we create temporary patch data with the required
    * buffer_cells for computing tag values.  (We could give the real
    * data the required ghost cells, but that may affect the
    * regridding algorithm I'm testing.)
    */
   pdat::CellData<int> tmp_tag(pbox, 1, buffer_cells);

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
   front_box.grow(buffer_cells);
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
         // This cell has nodes on both sides of the front.  Tag it and the buffer_cells around it.
         hier::Box cell_and_buffer(cell_index, cell_index, blk0);
         cell_and_buffer.grow(buffer_cells);
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



bool SinusoidalFrontGenerator::packDerivedDataIntoDoubleBuffer(
   double* buffer,
   const hier::Patch& patch,
   const hier::Box& region,
   const std::string& variable_name,
   int depth_index) const
{
   (void)region;
   (void)depth_index;

   if (variable_name == "Distance to front") {
      pdat::NodeData<double> dist_data(patch.getBox(), 1, hier::IntVector(d_dim,
                                          0));
      computePatchData(patch, 0.0, &dist_data, NULL);
      pdat::NodeData<double>::iterator ciend(patch.getBox(), false);
      for (pdat::NodeData<double>::iterator ci(patch.getBox(), true);
           ci != ciend; ++ci) {
         *(buffer++) = dist_data(*ci);
      }
   } else if (variable_name == "Tag value") {
      pdat::CellData<int> tag_data(patch.getBox(), 1, hier::IntVector(d_dim, 0));
      computePatchData(patch, 0.0, NULL, &tag_data);
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
