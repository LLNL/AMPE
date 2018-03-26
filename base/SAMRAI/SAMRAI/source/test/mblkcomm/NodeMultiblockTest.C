/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for node-centered patch data
 *
 ************************************************************************/

#include "NodeMultiblockTest.h"

#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/NodeVariable.h"

#include "MultiblockTester.h"

using namespace SAMRAI;

NodeMultiblockTest::NodeMultiblockTest(
   const string& object_name,
   const tbox::Dimension& dim,
   boost::shared_ptr<tbox::Database> main_input_db,
   bool do_refine,
   bool do_coarsen,
   const string& refine_option):
   PatchMultiblockTestStrategy(dim),
   d_dim(dim)
{
   NULL_USE(do_refine);
   NULL_USE(do_coarsen);
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(main_input_db);
   TBOX_ASSERT(!refine_option.empty());
#endif

   d_object_name = object_name;

   d_refine_option = refine_option;

   d_finest_level_number = main_input_db->
      getDatabase("PatchHierarchy")->
      getInteger("max_levels") - 1;

   char geom_name[32];

   sprintf(geom_name, "BlockGridGeometry");

   if (main_input_db->keyExists(geom_name)) {
      getGridGeometry().reset(
         new geom::GridGeometry(
            dim,
            geom_name,
            main_input_db->getDatabase(geom_name)));

   } else {
      TBOX_ERROR("NodeMultiblockTest: could not find entry `"
         << geom_name << "' in input.");
   }

   readTestInput(main_input_db->getDatabase("NodeMultiblockTest"));
}

NodeMultiblockTest::~NodeMultiblockTest()
{
}

void NodeMultiblockTest::readTestInput(
   boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(db);
#endif

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));
}

void NodeMultiblockTest::registerVariables(
   MultiblockTester* commtest)
{
   TBOX_ASSERT(commtest != (MultiblockTester *)NULL);

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i].reset(
         new pdat::NodeVariable<double>(d_dim,
                                        d_variable_src_name[i],
                                        d_variable_depth[i]));

      commtest->registerVariable(d_variables[i],
         d_variables[i],
         d_variable_src_ghosts[i],
         d_variable_dst_ghosts[i],
         getGridGeometry(),
         d_variable_refine_op[i]);

   }

}

void NodeMultiblockTest::initializeDataOnPatch(
   hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   const hier::BlockId& block_id,
   char src_or_dst)
{
   NULL_USE(hierarchy);
   NULL_USE(src_or_dst);

   if ((d_refine_option == "INTERIOR_FROM_SAME_LEVEL")
       || ((d_refine_option == "INTERIOR_FROM_COARSER_LEVEL")
           && (level_number < d_finest_level_number))) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::NodeData<double> > node_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = node_data->getGhostBox();

         node_data->fillAll((double)block_id.getBlockValue());

      }
   }
}

void NodeMultiblockTest::tagCellsToRefine(
   hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   int tag_index)
{
   NULL_USE(hierarchy);

   /*
    * Base class sets tags in box array for each level.
    */
   tagCellsInInputBoxes(patch, level_number, tag_index);

}

void NodeMultiblockTest::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw_to_fill) const
{
   NULL_USE(time);

   boost::shared_ptr<hier::PatchGeometry> pgeom(patch.getPatchGeometry());

   const tbox::Array<hier::BoundaryBox> node_bdry =
      pgeom->getCodimensionBoundaries(d_dim.getValue());
   const int num_node_bdry_boxes = node_bdry.getSize();

   tbox::Array<hier::BoundaryBox> edge_bdry;
   int num_edge_bdry_boxes = 0;
   if (d_dim > tbox::Dimension(1)) {
      edge_bdry = pgeom->getCodimensionBoundaries(d_dim.getValue() - 1);
      num_edge_bdry_boxes = edge_bdry.getSize();
   }

   tbox::Array<hier::BoundaryBox> face_bdry;
   int num_face_bdry_boxes = 0;
   if (d_dim == tbox::Dimension(3)) {
      face_bdry = pgeom->getCodimensionBoundaries(d_dim.getValue() - 2);
      num_face_bdry_boxes = face_bdry.getSize();
   }

   for (int i = 0; i < d_variables.getSize(); i++) {

      boost::shared_ptr<pdat::NodeData<double> > node_data(
         patch.getPatchData(d_variables[i], getDataContext()),
         boost::detail::dynamic_cast_tag());

      /*
       * Set node boundary data.
       */
      for (int nb = 0; nb < num_node_bdry_boxes; nb++) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[nb],
               patch.getBox(),
               gcw_to_fill);

         hier::Box patch_node_box =
            pdat::NodeGeometry::toNodeBox(patch.getBox());
         if (!node_bdry[nb].getIsMultiblockSingularity()) {
            pdat::NodeIterator niend(fill_box, false);
            for (pdat::NodeIterator ni(fill_box, true); ni != niend; ++ni) {
               if (!patch_node_box.contains(*ni)) {
                  for (int d = 0; d < node_data->getDepth(); d++) {
                     (*node_data)(*ni, d) =
                        (double)(node_bdry[nb].getLocationIndex() + 100);
                  }
               }
            }
         }
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int eb = 0; eb < num_edge_bdry_boxes; eb++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[eb],
                  patch.getBox(),
                  gcw_to_fill);

            hier::Box patch_node_box =
               pdat::NodeGeometry::toNodeBox(patch.getBox());
            hier::Index plower(patch_node_box.lower());
            hier::Index pupper(patch_node_box.upper());

            if (!edge_bdry[eb].getIsMultiblockSingularity()) {
               pdat::NodeIterator niend(fill_box, false);
               for (pdat::NodeIterator ni(fill_box, true); ni != niend; ++ni) {
                  if (!patch_node_box.contains(*ni)) {
                     bool use_index = true;
                     for (int n = 0; n < d_dim.getValue(); n++) {
                        if (edge_bdry[eb].getBox().numberCells(n) == 1) {
                           if ((*ni)(n) == plower(n) || (*ni)(n) ==
                               pupper(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }

                     if (use_index) {
                        for (int d = 0; d < node_data->getDepth(); d++) {
                           (*node_data)(*ni, d) =
                              (double)(edge_bdry[eb].getLocationIndex() + 100);
                        }
                     }
                  }
               }
            }
         }
      }

      if (d_dim == tbox::Dimension(3)) {
         /*
          * Set face boundary data.
          */
         for (int fb = 0; fb < num_face_bdry_boxes; fb++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(face_bdry[fb],
                  patch.getBox(),
                  gcw_to_fill);

            hier::Box patch_node_box =
               pdat::NodeGeometry::toNodeBox(patch.getBox());
            hier::Index plower(patch_node_box.lower());
            hier::Index pupper(patch_node_box.upper());

            if (!face_bdry[fb].getIsMultiblockSingularity()) {
               pdat::NodeIterator niend(fill_box, false);
               for (pdat::NodeIterator ni(fill_box, true); ni != niend; ++ni) {
                  if (!patch_node_box.contains(*ni)) {
                     bool use_index = true;
                     for (int n = 0; n < d_dim.getValue(); n++) {
                        if (face_bdry[fb].getBox().numberCells(n) == 1) {
                           if ((*ni)(n) == plower(n) || (*ni)(n) ==
                               pupper(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }

                     if (use_index) {
                        for (int d = 0; d < node_data->getDepth(); d++) {
                           (*node_data)(*ni, d) =
                              (double)(face_bdry[fb].getLocationIndex() + 100);
                        }
                     }
                  }
               }
            }
         }
      }

   }

}

void NodeMultiblockTest::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   const hier::Connector& dst_to_encon,
   const hier::Box& fill_box,
   const hier::BoundaryBox& bbox,
   const boost::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
   const tbox::Dimension& dim = fill_box.getDim();

   const hier::BoxId& dst_mb_id = patch.getBox().getId();
   const hier::BlockId& patch_blk_id = patch.getBox().getBlockId();

   const std::list<hier::BaseGridGeometry::Neighbor>& neighbors =
      grid_geometry->getNeighbors(patch_blk_id);

   for (int i = 0; i < d_variables.getSize(); i++) {

      boost::shared_ptr<pdat::NodeData<double> > node_data(
         patch.getPatchData(d_variables[i], getDataContext()),
         boost::detail::dynamic_cast_tag());

      hier::Box sing_fill_box(node_data->getGhostBox() * fill_box);

      int depth = node_data->getDepth();

      hier::Box pbox(pdat::NodeGeometry::toNodeBox(patch.getBox()));

      hier::Index plower(pbox.lower());
      hier::Index pupper(pbox.upper());

      pdat::NodeIterator niend(sing_fill_box, false);
      for (pdat::NodeIterator ni(sing_fill_box, true); ni != niend; ++ni) {
         bool use_index = true;
         for (int n = 0; n < d_dim.getValue(); n++) {
            if (bbox.getBox().numberCells(n) == 1) {
               if ((*ni)(n) == plower(n) || (*ni)(n) == pupper(n)) {
                  use_index = false;
                  break;
               }
            }
         }
         if (use_index) {
            for (int d = 0; d < depth; d++) {
               (*node_data)(*ni, d) = 0.0;
            }
         }
      }

      int num_encon_used = 0;

      if (grid_geometry->hasEnhancedConnectivity()) {

         hier::Connector::ConstNeighborhoodIterator ni =
            dst_to_encon.findLocal(dst_mb_id);

         if (ni != dst_to_encon.end()) {

            for (hier::Connector::ConstNeighborIterator ei = dst_to_encon.begin(ni);
                 ei != dst_to_encon.end(ni); ++ei) {

               const hier::BlockId& encon_blk_id = ei->getBlockId();
               boost::shared_ptr<hier::Patch> encon_patch(
                  encon_level.getPatch(ei->getId()));

               hier::Transformation::RotationIdentifier rotation =
                  hier::Transformation::NO_ROTATE;
               hier::IntVector offset(dim);

               for (std::list<hier::BaseGridGeometry::Neighbor>::const_iterator
                    nbri(neighbors.begin()); nbri != neighbors.end(); nbri++) {

                  if (nbri->getBlockId() == encon_blk_id) {
                     rotation = nbri->getRotationIdentifier();
                     offset = nbri->getShift();
                     break;
                  }
               }

               offset *= patch.getPatchGeometry()->getRatio();

               hier::Transformation transformation(rotation, offset,
                                                   encon_blk_id,
                                                   patch_blk_id);
               hier::Box encon_patch_box(encon_patch->getBox());
               transformation.transform(encon_patch_box);

               hier::Box encon_fill_box(encon_patch_box * sing_fill_box);
               if (!encon_fill_box.empty()) {

                  const hier::Transformation::RotationIdentifier back_rotate =
                     hier::Transformation::getReverseRotationIdentifier(
                        rotation, dim);

                  hier::IntVector back_shift(dim);

                  hier::Transformation::calculateReverseShift(
                     back_shift, offset, rotation);

                  hier::Transformation back_trans(back_rotate, back_shift,
                                                  patch_blk_id,
                                                  encon_blk_id);

                  boost::shared_ptr<pdat::NodeData<double> > sing_data(
                     encon_patch->getPatchData(d_variables[i], getDataContext()),
                     boost::detail::dynamic_cast_tag());

                  pdat::NodeIterator ciend(sing_fill_box, false);
                  for (pdat::NodeIterator ci(sing_fill_box, true);
                       ci != ciend; ++ci) {
                     bool use_index = true;
                     for (int n = 0; n < d_dim.getValue(); n++) {
                        if (bbox.getBox().numberCells(n) == 1) {
                           if ((*ci)(n) == plower(n) || (*ci)(n) == pupper(n)) {
                              use_index = false;
                              break;
                           }
                        }
                     }
                     if (use_index) {
                        pdat::NodeIndex src_index(*ci);
                        pdat::NodeGeometry::transform(src_index, back_trans);
                        for (int d = 0; d < depth; d++) {
                           (*node_data)(*ci, d) += (*sing_data)(src_index, d);
                        }
                     }
                  }
                  ++num_encon_used;
               }
            }
         }
      }

      if (num_encon_used) {
         pdat::NodeIterator ciend(sing_fill_box, false);
         for (pdat::NodeIterator ci(sing_fill_box, true); ci != ciend; ++ci) {
            bool use_index = true;
            for (int n = 0; n < d_dim.getValue(); n++) {
               if (bbox.getBox().numberCells(n) == 1) {
                  if ((*ci)(n) == plower(n) || (*ci)(n) == pupper(n)) {
                     use_index = false;
                     break;
                  }
               }
            }
            if (use_index) {
               for (int d = 0; d < depth; d++) {
                  (*node_data)(*ci, d) /= num_encon_used;
               }
            }
         }

      } else {

         /*
          * In cases of reduced connectivity, there are no other blocks
          * from which to acquire data.
          */

         pdat::NodeIterator ciend(sing_fill_box, false);
         for (pdat::NodeIterator ci(sing_fill_box, true); ci != ciend; ++ci) {
            bool use_index = true;
            for (int n = 0; n < d_dim.getValue(); n++) {
               if (bbox.getBox().numberCells(n) == 1) {
                  if ((*ci)(n) == plower(n) || (*ci)(n) == pupper(n)) {
                     use_index = false;
                     break;
                  }
               }
            }
            if (use_index) {
               for (int d = 0; d < depth; d++) {
                  (*node_data)(*ci,
                               d) = (double)bbox.getLocationIndex() + 200.0;
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Verify results of communication operations.  This test must be
 * consistent with data initialization and boundary operations above.
 *
 *************************************************************************
 */
bool NodeMultiblockTest::verifyResults(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   const int level_number,
   const hier::BlockId& block_id)
{

   tbox::plog << "\nEntering NodeMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl;

   hier::IntVector tgcw(d_dim, 0);
   for (int i = 0; i < d_variables.getSize(); i++) {
      tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
         getGhostCellWidth());
   }
   hier::Box pbox = patch.getBox();

   boost::shared_ptr<pdat::NodeData<double> > solution(
      new pdat::NodeData<double>(pbox, 1, tgcw));

   hier::Box tbox(pbox);
   tbox.grow(tgcw);

   const std::list<hier::BaseGridGeometry::Neighbor>& neighbors =
      hierarchy->getGridGeometry()->getNeighbors(block_id);
   hier::BoxContainer singularity(
      hierarchy->getGridGeometry()->getSingularityBoxContainer(block_id));

   hier::IntVector ratio =
      hierarchy->getPatchLevel(level_number)->getRatioToLevelZero();

   singularity.refine(ratio);

   bool test_failed = false;

   for (int i = 0; i < d_variables.getSize(); i++) {

      double correct = (double)block_id.getBlockValue();

      boost::shared_ptr<pdat::NodeData<double> > node_data(
         patch.getPatchData(d_variables[i], getDataContext()),
         boost::detail::dynamic_cast_tag());
      int depth = node_data->getDepth();

      hier::Box interior_box(pbox);
      interior_box.grow(hier::IntVector(d_dim, -1));

      pdat::NodeIterator ciend(interior_box, false);
      for (pdat::NodeIterator ci(interior_box, true); ci != ciend; ++ci) {
         for (int d = 0; d < depth; d++) {
            double result = (*node_data)(*ci, d);

            if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
               tbox::perr << "Test FAILED: ...."
                          << " : node index = " << *ci << endl;
               tbox::perr << "    Variable = " << d_variable_src_name[i]
                          << " : depth index = " << d << endl;
               tbox::perr << "    result = " << result
                          << " : correct = " << correct << endl;
               test_failed = true;
            }
         }
      }

      boost::shared_ptr<hier::PatchGeometry> pgeom(patch.getPatchGeometry());

      hier::Box gbox = node_data->getGhostBox();

      hier::Box patch_node_box =
         pdat::NodeGeometry::toNodeBox(pbox);

      hier::BoxContainer sing_node_boxlist;
      for (hier::BoxContainer::iterator si(singularity);
           si != singularity.end(); ++si) {
         sing_node_boxlist.pushFront(pdat::NodeGeometry::toNodeBox(*si));
      }

      hier::BoxContainer tested_neighbors;

      for (std::list<hier::BaseGridGeometry::Neighbor>::const_iterator
           ne(neighbors.begin()); ne != neighbors.end(); ne++) {

         correct = ne->getBlockId().getBlockValue();

         hier::BoxContainer neighbor_ghost(ne->getTransformedDomain());

         hier::BoxContainer neighbor_node_ghost;
         for (hier::BoxContainer::iterator nn(neighbor_ghost);
              nn != neighbor_ghost.end(); ++nn) {
            hier::Box neighbor_ghost_interior(
               pdat::NodeGeometry::toNodeBox(*nn));
            neighbor_ghost_interior.grow(-hier::IntVector::getOne(d_dim));
            neighbor_node_ghost.pushFront(neighbor_ghost_interior);
         }

         neighbor_node_ghost.refine(ratio);

         neighbor_node_ghost.intersectBoxes(
            pdat::NodeGeometry::toNodeBox(gbox));

         neighbor_node_ghost.removeIntersections(sing_node_boxlist);
         neighbor_node_ghost.removeIntersections(tested_neighbors);

         for (hier::BoxContainer::iterator ng(neighbor_node_ghost);
              ng != neighbor_node_ghost.end(); ++ng) {

            hier::Box::iterator ciend(*ng, false);
            for (hier::Box::iterator ci(*ng, true); ci != ciend; ++ci) {
               pdat::NodeIndex ni(*ci, hier::IntVector(d_dim, 0));
               if (!patch_node_box.contains(ni)) {
                  for (int d = 0; d < depth; d++) {
                     double result = (*node_data)(ni, d);

                     if (!tbox::MathUtilities<double>::equalEps(correct,
                            result)) {
                        tbox::perr << "Test FAILED: ...."
                                   << " : node index = " << ni << endl;
                        tbox::perr << "  Variable = " << d_variable_src_name[i]
                                   << " : depth index = " << d << endl;
                        tbox::perr << "    result = " << result
                                   << " : correct = " << correct << endl;
                        test_failed = true;
                     }
                  }
               }
            }
         }
         tested_neighbors.spliceBack(neighbor_node_ghost);
      }

      for (int b = 0; b < d_dim.getValue(); b++) {
         tbox::Array<hier::BoundaryBox> bdry =
            pgeom->getCodimensionBoundaries(b + 1);

         for (int k = 0; k < bdry.size(); k++) {
            hier::Box fill_box = pgeom->getBoundaryFillBox(bdry[k],
                  patch.getBox(),
                  tgcw);
            fill_box = fill_box * gbox;

            if (bdry[k].getIsMultiblockSingularity()) {
               correct = 0.0;

               int num_sing_neighbors = 0;
               for (std::list<hier::BaseGridGeometry::Neighbor>::const_iterator
                    ns(neighbors.begin()); ns != neighbors.end(); ns++) {
                  if (ns->isSingularity()) {
                     hier::BoxContainer neighbor_ghost(
                        ns->getTransformedDomain());
                     neighbor_ghost.refine(ratio);
                     neighbor_ghost.intersectBoxes(fill_box);
                     if (neighbor_ghost.size()) {
                        num_sing_neighbors++;
                        correct += ns->getBlockId().getBlockValue();
                     }
                  }
               }

               if (num_sing_neighbors == 0) {

                  correct = (double)bdry[k].getLocationIndex() + 200.0;

               } else {

                  correct /= (double)num_sing_neighbors;

               }

            } else {
               correct = (double)(bdry[k].getLocationIndex() + 100);
            }

            pdat::NodeIterator ciend(fill_box, false);
            for (pdat::NodeIterator ci(fill_box, true); ci != ciend; ++ci) {

               if (!patch_node_box.contains(*ci)) {

                  bool use_index = true;
                  for (int n = 0; n < d_dim.getValue(); n++) {
                     if (bdry[k].getBox().numberCells(n) == 1) {
                        if ((*ci)(n) == patch_node_box.lower() (n) ||
                            (*ci)(n) == patch_node_box.upper() (n)) {
                           use_index = false;
                           break;
                        }
                     }
                  }

                  if (use_index) {
                     for (int d = 0; d < depth; d++) {
                        double result = (*node_data)(*ci, d);

                        if (!tbox::MathUtilities<double>::equalEps(correct,
                               result)) {
                           tbox::perr << "Test FAILED: ...."
                                      << " : node index = " << *ci << endl;
                           tbox::perr << "  Variable = "
                                      << d_variable_src_name[i]
                                      << " : depth index = " << d << endl;
                           tbox::perr << "    result = " << result
                                      << " : correct = " << correct << endl;
                           test_failed = true;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (!test_failed) {
      tbox::plog << "NodeMultiblockTest Successful!" << endl;
   } else {
      tbox::perr << "Multiblock NodeMultiblockTest FAILED: .\n" << endl;
   }

   solution.reset();   // just to be anal...

   tbox::plog << "\nExiting NodeMultiblockTest::verifyResults..." << endl;
   tbox::plog << "level_number = " << level_number << endl;
   tbox::plog << "Patch box = " << patch.getBox() << endl << endl;

   return !test_failed;
}
