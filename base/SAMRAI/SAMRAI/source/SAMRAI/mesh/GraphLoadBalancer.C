/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/
#include "SAMRAI/mesh/GraphLoadBalancer.h"

#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"

#ifdef HAVE_PTSCOTCH
#include "ptscotch.h"
#endif

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 * GraphLoadBalancer constructor.
 *************************************************************************
 */

GraphLoadBalancer::GraphLoadBalancer(
   const tbox::Dimension& dim,
   const std::string& name,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(name),
   d_target_box_size(dim, 0),
   d_coalesce_boxes(true),
   d_tile_size(dim, 1),
   d_min_size(dim),
   d_cut_factor(dim),
   d_bad_interval(dim)
{
#ifndef HAVE_PTSCOTCH
   TBOX_WARNING(
      "SAMRAI configured without PT-Scotch.  GraphLoadBalancer will not repartition boxes.");
#endif

   TBOX_ASSERT(!name.empty());
   getFromInput(input_db);

}

GraphLoadBalancer::~GraphLoadBalancer()
{
}

/*
 *************************************************************************
 * Load balance and redistribute the level
 *************************************************************************
 */

void
GraphLoadBalancer::loadBalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector* balance_to_anchor,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
   NULL_USE(rank_group);

   // Set effective_cut_factor to least common multiple of cut_factor and d_tile_size.
   hier::IntVector effective_cut_factor = cut_factor;
   if ( d_tile_size != hier::IntVector::getOne(d_dim) ) {
      const size_t nblocks = hierarchy->getGridGeometry()->getNumberBlocks();
      for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
         for ( unsigned int d=0; d<d_dim.getValue(); ++d ) {
            while ( effective_cut_factor(b,d)/d_tile_size[d]*d_tile_size[d] != effective_cut_factor(b,d) ) {
               effective_cut_factor(b,d) += cut_factor[d];
            }
         }
      }
   }

#ifdef HAVE_PTSCOTCH
   d_min_size = min_size;
   d_bad_interval = bad_interval;
   d_cut_factor = effective_cut_factor;
   d_block_domain_boxes.clear();
   size_t nblocks =
      domain_box_level.getGridGeometry()->getNumberBlocks();
   d_block_domain_boxes.resize(nblocks);

   if (nblocks == 1) {
      domain_box_level.getGlobalBoxes(d_block_domain_boxes[0]);
      d_block_domain_boxes[0].refine(balance_box_level.getRefinementRatio());
   } else {
      for (int b = 0; b < nblocks; ++b) {
         d_block_domain_boxes[b] = hier::BoxContainer(
               domain_box_level.getGlobalBoxes(), hier::BlockId(b));

         d_block_domain_boxes[b].refine(balance_box_level.getRefinementRatio());
      }
   }

   /*
    * Periodic image Box should be ignored during load balancing
    * because they have no real work.  The load-balanced results
    * should contain no periodic images.
    *
    * To avoid need for special logic to skip periodic images while
    * load balancing, we just remove periodic images in the
    * balance_box_level and all periodic edges in
    * anchor<==>balance.
    */

   balance_box_level.removePeriodicImageBoxes();
   if (balance_to_anchor) {

      balance_to_anchor->getTranspose().removePeriodicRelationships();
      balance_to_anchor->getTranspose().setHead(balance_box_level, true);

      balance_to_anchor->removePeriodicRelationships();
      balance_to_anchor->setBase(balance_box_level, true);

   }

   hier::Connector& anchor_to_balance = balance_to_anchor->getTranspose();

   SCOTCH_Dgraph* graph = SCOTCH_dgraphAlloc();
   SCOTCH_dgraphInit(graph, balance_box_level.getMPI().getCommunicator());

   const tbox::Dimension& dim = balance_box_level.getDim();

   boost::shared_ptr<hier::Connector> balance_to_balance;

   int dont_do_graph = 0;

   const hier::BoxContainer& boxes = balance_box_level.getBoxes();
   const tbox::SAMRAI_MPI& my_mpi = balance_box_level.getMPI();

   hier::IntVector target_size(dim);

   if (d_target_box_size == hier::IntVector::getZero(dim)) {
      /*
       * Heuristic to choose a target box size if none was given in
       * input.
       */
      for (int d = 0; d < dim.getValue(); ++d) {
         std::set<int> choices;

         choices.insert(min_size(d) * min_size(d));
         choices.insert(max_size(d) / 2);
         choices.insert((min_size(d) + max_size(d)) / 2);

         std::set<int>::iterator med_itr = choices.begin();
         ++med_itr;
         target_size(d) = *med_itr;
      }
   } else {
      target_size = d_target_box_size;
   }

   for (int d = 0; d < dim.getValue(); ++d) {
      if (target_size(d) > max_size(d)) target_size(d) = max_size(d);
      if (target_size(d) < min_size(d)) target_size(d) = min_size(d);
   }

   /*
    * Chop boxes to be at or under target size.
    */
   chopBoxes(balance_box_level, &anchor_to_balance, target_size);

   const hier::MappingConnectorAlgorithm mca;

   /*
    * PT-Scotch requires globally sequenced LocalId values
    */
   renumberBoxes(balance_box_level,
      anchor_to_balance,
      mca);

   const hier::OverlapConnectorAlgorithm oca;

   /*
    * The balance_to_balance Connector describes adjacency within the
    * level. The connector width should be nonzero if this load balancer
    * is called from GriddingAlgorithm. If called from non-standard contexts,
    * the non-scalable findOverlaps method may be called.
    */
   if (balance_to_anchor->getConnectorWidth() > hier::IntVector::getZero(dim)) {
      oca.bridgeWithNesting(
         balance_to_balance,
         *balance_to_anchor,
         anchor_to_balance,
         hier::IntVector::getZero(dim),
         hier::IntVector::getZero(dim),
         hier::IntVector::getOne(dim),
         false);
   } else {
      balance_to_balance.reset(new hier::Connector(dim));
      balance_to_balance->clearNeighborhoods();
      balance_to_balance->setBase(balance_box_level);
      balance_to_balance->setHead(balance_box_level);
      balance_to_balance->setWidth(hier::IntVector::getOne(dim), true);

      oca.findOverlaps(*balance_to_balance);
   }

   std::map<hier::BoxId, bool> has_nabrs;
   for (hier::BoxContainer::const_iterator bi = boxes.begin();
        bi != boxes.end(); ++bi) {
      const hier::BoxId& box_id = bi->getBoxId();

      if (balance_to_balance->hasNeighborSet(box_id)) {
         hier::Connector::ConstNeighborhoodIterator nh =
            balance_to_balance->findLocal(box_id);

         bool has_non_trivial = false;

         for (hier::Connector::ConstNeighborIterator na = balance_to_balance->begin(nh);
              na != balance_to_balance->end(nh); ++na) {
            if (na->getBoxId() != box_id) {
               has_non_trivial = true;
               break;
            }
         }

         if (!has_non_trivial) {
            has_nabrs[box_id] = false;
         } else {
            has_nabrs[box_id] = true;
         }

      } else {
         has_nabrs[box_id] = false;
      }
   }

   std::vector<SCOTCH_Num> edgeloctab;
   std::vector<SCOTCH_Num> edloloctab;
   std::vector<SCOTCH_Num> vertloctab;
   std::vector<SCOTCH_Num> veloloctab;
   std::vector<SCOTCH_Num> vendloctab;
   std::vector<SCOTCH_Num> extra_nabrs;

   vertloctab.push_back(0);

   if (!boxes.empty()) {
      const hier::BoxId& back_id = boxes.back().getBoxId();

      for (hier::Connector::ConstNeighborhoodIterator ei = balance_to_balance->begin();
           ei != balance_to_balance->end(); ++ei) {
         const hier::BoxId& box_id = *ei;
         const hier::Box& box = *balance_box_level.getBox(box_id);
         hier::Box node_box(box);
         node_box.upper() += hier::IntVector::getOne(dim);
         if (has_nabrs[box_id]) {
            /*
             * If the box has neighbors, create graph edges to the neighbors.
             */
            for (hier::Connector::ConstNeighborIterator na = balance_to_balance->begin(ei);
                 na != balance_to_balance->end(ei); ++na) {
               const hier::Box& nbr_box = *na;
               if (nbr_box.getBoxId() != box_id) {
                  const hier::LocalId& local_id = nbr_box.getLocalId();
                  edgeloctab.push_back(local_id.getValue());
                  hier::Box node_nbr(nbr_box);
                  node_nbr.upper() += hier::IntVector::getOne(dim);
                  SCOTCH_Num edge_wgt = 1;
                  if (node_box.getBlockId() == node_nbr.getBlockId())
                     edge_wgt = (node_box * node_nbr).size();
                  edloloctab.push_back(edge_wgt);
               }
            }
            if (box_id == back_id && !extra_nabrs.empty()) {
               for (std::vector<SCOTCH_Num>::const_iterator extra_itr =
                       extra_nabrs.begin(); extra_itr != extra_nabrs.end(); ++extra_itr) {
                  edgeloctab.push_back(*extra_itr);
                  edloloctab.push_back(1);
               }
            }
         } else if (boxes.size() > 1) {
            /*
             * When the box has no neighbors, arbitrarily create an edge
             * to connect it to the graph.
             */
            if (back_id != box_id) {
               edgeloctab.push_back(back_id.getLocalId().getValue());
               edloloctab.push_back(1);
               extra_nabrs.push_back(box_id.getLocalId().getValue());
            } else {
               SCOTCH_Num my_id = box_id.getLocalId().getValue();
               edgeloctab.push_back(my_id);
               edloloctab.push_back(1);
               vertloctab.pop_back();
               vendloctab.pop_back();
               vertloctab.push_back(edgeloctab.size());
               vendloctab.push_back(edgeloctab.size());

               for (std::vector<SCOTCH_Num>::const_iterator extra_itr =
                       extra_nabrs.begin(); extra_itr != extra_nabrs.end(); ++extra_itr) {
                  if (*extra_itr != my_id - 1) {
                     edgeloctab.push_back(*extra_itr);
                     edloloctab.push_back(1);
                  }
               }

            }
         }

         vertloctab.push_back(edgeloctab.size());
         vendloctab.push_back(edgeloctab.size());
         veloloctab.push_back(box.size());
      }
   }
   vertloctab.pop_back();
   if (vertloctab.empty()) {
      vertloctab.push_back(0);
   }
   if (veloloctab.empty()) {
      veloloctab.push_back(0);
   }
   if (vendloctab.empty()) {
      vendloctab.push_back(0);
   }
   int edgelocsize = edgeloctab.size();
   if (edgeloctab.empty()) {
      edgeloctab.push_back(0);
   }
   if (edloloctab.empty()) {
      edloloctab.push_back(0);
   }

   SCOTCH_Num* vertloc_ptr;
   SCOTCH_Num* veloloc_ptr;
   SCOTCH_Num* vendloc_ptr;
   SCOTCH_Num* edgeloc_ptr;
   SCOTCH_Num* edloloc_ptr;
   if (!vertloctab.empty()) {
      vertloc_ptr = &vertloctab[0];
   } else {
      vertloc_ptr = 0;
   }
   if (!veloloctab.empty()) {
      veloloc_ptr = &veloloctab[0];
   } else {
      veloloc_ptr = 0;
   }
   if (!vendloctab.empty()) {
      vendloc_ptr = &vendloctab[0];
   } else {
      vendloc_ptr = 0;
   }
   if (!edgeloctab.empty()) {
      edgeloc_ptr = &edgeloctab[0];
   } else {
      edgeloc_ptr = 0;
   }
   if (!edloloctab.empty()) {
      edloloc_ptr = &edloloctab[0];
   } else {
      edloloc_ptr = 0;
   }

   SCOTCH_Num baseval = 0;

   int nlocvert = vertloctab.size(); //-1;
   if (nlocvert < 0) nlocvert = 0;
   if (boxes.empty()) nlocvert = 0;

   SCOTCH_dgraphBuild(graph,
      baseval,
      nlocvert,
      nlocvert,
      vertloc_ptr,
      vendloc_ptr,
      veloloc_ptr,                 // (node weights
      0,                 // vlblocltab (labels)
      edgelocsize,                 // zero if local proc has no nodes
      edgeloctab.size(),
      edgeloc_ptr,
      0,                 // edgegsttab (ghosts)
      edloloc_ptr);                // edge weights

   SCOTCH_Num partloctab[nlocvert];

   SCOTCH_Strat stradat;
   SCOTCH_stratInit(&stradat);
   SCOTCH_stratDgraphMapBuild(&stradat,
      SCOTCH_STRATBALANCE,
      balance_box_level.getMPI().getSize(),
      0,
      0.00);

   int err2 = SCOTCH_dgraphPart(graph,
         balance_box_level.getMPI().getSize(),
         &stradat,
         partloctab);

   int num_ranks = my_mpi.getSize();
   int my_rank = my_mpi.getRank();
   int boxes_on_rank[num_ranks];
   for (int rank = 0; rank < num_ranks; ++rank) {
      if (rank != my_rank) {
         boxes_on_rank[rank] = 0;
      } else {
         boxes_on_rank[rank] = boxes.size();
      }
   }

   my_mpi.AllReduce(boxes_on_rank,
      num_ranks,
      MPI_SUM);

   int start_box = 0;
   for (int rank = 0; rank < my_rank; ++rank) {
      start_box += boxes_on_rank[rank];
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!boxes.empty()) {
      TBOX_ASSERT(start_box == boxes.begin()->getLocalId().getValue());
   }
#endif

   int global_num_boxes = balance_box_level.getGlobalNumberOfBoxes();

   int old_global_partition[global_num_boxes];
   int new_global_partition[global_num_boxes];
   for (int b = 0; b < global_num_boxes; ++b) {
      old_global_partition[b] = 0;
      new_global_partition[b] = 0;
   }
   for (int b = start_box; b < start_box + nlocvert; ++b) {
      old_global_partition[b] = my_rank;
      new_global_partition[b] = partloctab[b - start_box];
   }

   my_mpi.AllReduce(old_global_partition,
      global_num_boxes,
      MPI_SUM);
   my_mpi.AllReduce(new_global_partition,
      global_num_boxes,
      MPI_SUM);

   hier::BoxLevel graph_level(balance_box_level.getRefinementRatio(),
                              balance_box_level.getGridGeometry(),
                              balance_box_level.getMPI(),
                              hier::BoxLevel::DISTRIBUTED);

   tbox::AsyncCommStage send_stage;
   std::map<int, tbox::AsyncCommPeer<char> *> send_comms;
   std::set<int> send_procs;
   tbox::AsyncCommStage recv_stage;
   std::map<int, tbox::AsyncCommPeer<char> *> recv_comms;

   setupAsyncCommObjects(
      send_stage,
      send_comms,
      send_procs,
      recv_stage,
      recv_comms,
      old_global_partition,
      new_global_partition,
      global_num_boxes,
      my_rank,
      my_mpi);

   for (std::map<int, tbox::AsyncCommPeer<char> *>::iterator ri =
           recv_comms.begin(); ri != recv_comms.end(); ++ri) {
      tbox::AsyncCommPeer<char> *& recv_peer = ri->second;
      recv_peer->beginRecv();
      if (recv_peer->isDone()) {
         recv_peer->pushToCompletionQueue();
      }
   }

   /*
    * Boxes that stay on the local process are added to graph_level.
    * Otherwise BoxInTransit objecst are set up to be communicated.
    */
   std::list<BoxInTransit> transit_boxes;
   for (hier::BoxContainer::const_iterator itr = boxes.begin();
        itr != boxes.end(); ++itr) {
      const int& local_id = itr->getLocalId().getValue();
      TBOX_ASSERT(old_global_partition[local_id] == my_rank);

      if (new_global_partition[local_id] == my_rank) {
         graph_level.addBoxWithoutUpdate(*itr);
      } else {
         transit_boxes.push_back(BoxInTransit(*itr,
               *itr,
               new_global_partition[local_id],
               hier::LocalId(local_id)));
      }
   }

   /*
    * Determine number of boxes to send.
    */
   int num_send_boxes[num_ranks];
   for (int rank = 0; rank < num_ranks; ++rank) {
      num_send_boxes[rank] = 0;
   }
   for (int b = 0; b < global_num_boxes; ++b) {
      if (old_global_partition[b] == my_rank &&
          new_global_partition[b] != my_rank) {
         ++num_send_boxes[new_global_partition[b]];
      }
   }

   tbox::MessageStream* mstreams = new tbox::MessageStream[num_ranks];
   for (int rank = 0; rank < num_ranks; ++rank) {
      if (num_send_boxes[rank]) {
         mstreams[rank] << num_send_boxes[rank];
      }
   }

   /*
    * Place transit boxes in stream and send.
    */
   for (std::list<BoxInTransit>::const_iterator ti = transit_boxes.begin();
        ti != transit_boxes.end(); ++ti) {
      const BoxInTransit& transit_box = *ti;
      const int send_rank = transit_box.d_box.getOwnerRank();
      if (send_rank != transit_box.d_orig_box.getOwnerRank()) {

         transit_box.putToMessageStream(mstreams[send_rank]);

      }
   }

   int send_ct = static_cast<int>(send_procs.size());
   std::set<int>::const_iterator si = send_procs.lower_bound(my_rank + 1);
   for (int num_sent = 0; num_sent < send_ct; ++num_sent) {
      if (si == send_procs.end()) {
         si = send_procs.begin();
      }
      int send_rank = *si;
      tbox::AsyncCommPeer<char> *& send_peer = send_comms[send_rank];
      const tbox::MessageStream& msg = mstreams[send_rank];
      send_peer->beginSend(static_cast<const char *>(msg.getBufferStart()),
         static_cast<int>(msg.getCurrentSize()));
      ++si;
   }

   for (si = send_procs.begin(); si != send_procs.end(); ++si) {

      tbox::AsyncCommPeer<char> *& send_peer = send_comms[*si];
      send_peer->completeCurrentOperation();
   }

   delete[] mstreams;

   hier::MappingConnector balance_to_graph(dim);
   hier::MappingConnector graph_to_balance(dim);

   balance_to_graph.clearNeighborhoods();
   balance_to_graph.setBase(balance_box_level);
   balance_to_graph.setHead(graph_level);
   balance_to_graph.setWidth(hier::IntVector::getZero(dim), true);
   graph_to_balance.clearNeighborhoods();
   graph_to_balance.setBase(graph_level);
   graph_to_balance.setHead(balance_box_level);
   graph_to_balance.setWidth(hier::IntVector::getZero(dim), true);

   /*
    * balance_to_graph connector set up based on the boxes that were sent.
    */
   for (std::list<BoxInTransit>::iterator ti = transit_boxes.begin();
        ti != transit_boxes.end(); ++ti) {

      const BoxInTransit& sent_box = *ti;

      hier::Connector::NeighborhoodIterator base_box_itr =
         balance_to_graph.makeEmptyLocalNeighborhood(
            sent_box.d_orig_box.getBoxId());

      balance_to_graph.insertLocalNeighbor(
         sent_box.getBox(),
         base_box_itr);
   }

   /*
    * Complete receives and unpack streams.
    */
   std::list<BoxInTransit> received_transit_boxes;
   while (recv_stage.hasCompletedMembers() || recv_stage.advanceSome()) {

      tbox::AsyncCommPeer<char>* recv_peer =
         CPP_CAST<tbox::AsyncCommPeer<char> *>(recv_stage.popCompletionQueue());

      TBOX_ASSERT(recv_peer != 0);

      tbox::MessageStream mstream(recv_peer->getRecvSize(),
                                  tbox::MessageStream::Read,
                                  recv_peer->getRecvData(),
                                  false);

      int stream_size = static_cast<int>(mstream.getCurrentSize());
      int num_boxes = 0;
      mstream >> num_boxes;

      BoxInTransit received_box(d_dim);
      for (int i = 0; i < num_boxes; ++i) {
         received_box.getFromMessageStream(mstream);
         received_transit_boxes.push_back(received_box);
      }
   }

   TBOX_ASSERT(!recv_stage.hasPendingRequests());

   /*
    * Add received boxes to graph_level and add graph_to_balance neighbors.
    */
   for (std::list<BoxInTransit>::iterator ri = received_transit_boxes.begin();
        ri != received_transit_boxes.end(); ++ri) {

      const BoxInTransit& received_box = *ri;
      graph_level.addBoxWithoutUpdate(received_box.getBox());

      hier::Connector::NeighborhoodIterator base_box_itr =
         graph_to_balance.makeEmptyLocalNeighborhood(
            received_box.getBox().getBoxId());

      graph_to_balance.insertLocalNeighbor(
         received_box.d_orig_box,
         base_box_itr);
   }

   graph_level.finalize();

   balance_to_graph.setTranspose(&graph_to_balance, false);

   /*
    * balance_box_level modified to become the graph-based partition.
    */
   mca.modify(anchor_to_balance,
      balance_to_graph,
      &balance_box_level);

   if (d_coalesce_boxes) {
      /*
       * Coalesce local boxes and then chop to enforce max size.
       */

      coalesceBoxLevel(balance_box_level,
         anchor_to_balance,
         mca);

      hier::IntVector maxvector(d_dim, tbox::MathUtilities<int>::getMax());
      if (max_size != maxvector) {
         chopBoxes(balance_box_level, &anchor_to_balance, max_size);
      }
   }

   for (std::map<int, tbox::AsyncCommPeer<char> *>::iterator c_itr = send_comms.begin();
        c_itr != send_comms.end();
        ++c_itr) {

      if (c_itr->second)
         delete c_itr->second;

   }
   for (std::map<int, tbox::AsyncCommPeer<char> *>::iterator c_itr = recv_comms.begin();
        c_itr != recv_comms.end();
        ++c_itr) {

      if (c_itr->second)
         delete c_itr->second;

   }

#else
   NULL_USE(balance_box_level);
   NULL_USE(balance_to_anchor);
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(min_size);
   NULL_USE(max_size);
   NULL_USE(domain_box_level);
   NULL_USE(bad_interval);
   NULL_USE(cut_factor);
   TBOX_WARNING(
      "SAMRAI configured without PT-Scotch library:  GraphLoadBalancer calls will do nothing");

#endif

}

/*
 *  ***********************************************************************
 *  0 to N-1 renumbering of all boxes on the level
 *  ***********************************************************************
 */
void
GraphLoadBalancer::renumberBoxes(
   hier::BoxLevel& new_box_level,
   hier::Connector& anchor_to_balance,
   const hier::MappingConnectorAlgorithm& mca) const
{
   boost::shared_ptr<hier::MappingConnector> sorting_map;
   boost::shared_ptr<hier::BoxLevel> seq_box_level;
   hier::BoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.makeSortingMap(
      seq_box_level,
      sorting_map,
      new_box_level,
      false,          //sort_by_corners = false
      true);          //sequentialize_global_indices = true

   mca.modify(anchor_to_balance,
      *sorting_map,
      &new_box_level);

}

/*
 *  ***********************************************************************
 *  Coalesce local boxes
 *  ***********************************************************************
 */
void
GraphLoadBalancer::coalesceBoxLevel(
   hier::BoxLevel& level,
   hier::Connector& anchor_to_level,
   const hier::MappingConnectorAlgorithm& mca) const
{
   const hier::BoxContainer& level_boxes = level.getBoxes();

   hier::BoxLevel coalesced(level.getRefinementRatio(),
                            level.getGridGeometry(),
                            level.getMPI());
   hier::MappingConnector level_to_coalesced(level,
                                             coalesced,
                                             hier::IntVector::getZero(d_dim));

   if (!level_boxes.empty()) {

      const int my_rank = level_boxes.begin()->getOwnerRank();

      size_t nblocks = level.getGridGeometry()->getNumberBlocks();
      hier::LocalId local_id(0);

      for (hier::BlockId::block_t b = 0; b < nblocks; ++b) {
         hier::BlockId block_id(b);

         hier::BoxContainer block_boxes(level_boxes, block_id);

         if (!block_boxes.empty()) {
            block_boxes.unorder();
            block_boxes.coalesce();

            for (hier::BoxContainer::iterator bi = block_boxes.begin();
                 bi != block_boxes.end(); ++bi) {

               const hier::Box new_box(*bi,
                                       local_id++,
                                       my_rank);

               coalesced.addBoxWithoutUpdate(new_box);

            }
         }
      }
   }

   coalesced.finalize();

   const hier::BoxContainer& coalesced_boxes = coalesced.getBoxes();

   int was_coalesced = 0;
   if (coalesced_boxes.size() != level_boxes.size()) {
      was_coalesced = 1;
   }
   const tbox::SAMRAI_MPI& mpi = level.getMPI();
   mpi.AllReduce(&was_coalesced,
      1,
      MPI_SUM);

   if (was_coalesced != 0) {

      for (hier::BoxContainer::const_iterator li = level_boxes.begin();
           li != level_boxes.end(); ++li) {

         const hier::Box& old_box = *li;
         const hier::BlockId& block_id = old_box.getBlockId();

         hier::Connector::NeighborhoodIterator base_box_itr =
            level_to_coalesced.makeEmptyLocalNeighborhood(
               old_box.getBoxId());

         for (hier::BoxContainer::const_iterator bi = coalesced_boxes.begin();
              bi != coalesced_boxes.end(); ++bi) {
            if (bi->getBlockId() == block_id) {

               if (old_box.intersects(*bi)) {

                  level_to_coalesced.insertLocalNeighbor(
                     *bi,
                     base_box_itr);

               }
            }
         }
      }

      mca.modify(anchor_to_level,
         level_to_coalesced,
         &level);
   }
}

/*
 *  ***********************************************************************
 *  Chop boxes which are greater than max_size
 *  ***********************************************************************
 */

void
GraphLoadBalancer::chopBoxes(
   hier::BoxLevel& box_level,
   hier::Connector* anchor_to_level,
   const hier::IntVector& max_size) const
{
   const tbox::Dimension& dim = max_size.getDim();
   TBOX_ASSERT(!anchor_to_level || anchor_to_level->hasTranspose());
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY1(dim, box_level);

   const hier::IntVector& zero_vector(hier::IntVector::getZero(dim));

   hier::BoxLevel constrained(box_level.getRefinementRatio(),
                              box_level.getGridGeometry(),
                              box_level.getMPI());
   hier::MappingConnector unconstrained_to_constrained(box_level,
                                                       constrained,
                                                       zero_vector);

   const hier::BoxContainer& unconstrained_boxes = box_level.getBoxes();

   hier::LocalId next_available_index = box_level.getLastLocalId() + 1;

   for (hier::BoxContainer::const_iterator ni = unconstrained_boxes.begin();
        ni != unconstrained_boxes.end(); ++ni) {

      const hier::Box& box = *ni;

      const hier::IntVector box_size = box.numberCells();

      /*
       * If box already conform to max size constraint, keep it.
       * Else chop it up and keep the parts.
       */

      if (box_size <= max_size) {

         constrained.addBoxWithoutUpdate(box);

      } else {

         hier::BoxContainer chopped(box);
         hier::BoxUtilities::chopBoxes(
            chopped,
            max_size,
            d_min_size,
            d_cut_factor,
            d_bad_interval,
            d_block_domain_boxes[box.getBlockId().getBlockValue()]);
         TBOX_ASSERT(!chopped.empty());

         if (chopped.size() != 1) {

            hier::Connector::NeighborhoodIterator base_box_itr =
               unconstrained_to_constrained.makeEmptyLocalNeighborhood(
                  box.getBoxId());

            for (hier::BoxContainer::iterator li = chopped.begin();
                 li != chopped.end(); ++li) {

               const hier::Box new_box(*li,
                                       next_available_index++,
                                       box.getOwnerRank());
               TBOX_ASSERT(new_box.getBlockId() == ni->getBlockId());

               constrained.addBoxWithoutUpdate(new_box);

               unconstrained_to_constrained.insertLocalNeighbor(
                  new_box,
                  base_box_itr);

            }

         } else {
            TBOX_ASSERT(box.isSpatiallyEqual(chopped.front()));
            constrained.addBoxWithoutUpdate(box);
         }

      }

   }

   constrained.finalize();

   int was_chopped = 0;
   if (constrained.getBoxes().size() != box_level.getBoxes().size()) {
      was_chopped = 1;
   }
   const tbox::SAMRAI_MPI& mpi = box_level.getMPI();
   mpi.AllReduce(&was_chopped,
      1,
      MPI_SUM);

   if (was_chopped != 0) {
      hier::MappingConnectorAlgorithm mca;
      mca.setTimerPrefix(d_object_name);
      mca.modify(*anchor_to_level,
         unconstrained_to_constrained,
         &box_level);
   }
}

/*
 *  ***********************************************************************
 *  BoxInTransit constructors
 *  ***********************************************************************
 */

GraphLoadBalancer::BoxInTransit::BoxInTransit(
   const tbox::Dimension& dim):
   d_box(dim),
   d_orig_box(dim)
{
}

GraphLoadBalancer::BoxInTransit::BoxInTransit(
   const hier::Box& other,
   const hier::Box& box,
   int rank,
   hier::LocalId local_id):
   d_box(box, local_id, rank),
   d_orig_box(other)
{
}

/*
 *  ***********************************************************************
 *  Set up communications
 *  ***********************************************************************
 */

void
GraphLoadBalancer::setupAsyncCommObjects(
   tbox::AsyncCommStage& send_stage,
   std::map<int, tbox::AsyncCommPeer<char> *>& send_comms,
   std::set<int>& send_procs,
   tbox::AsyncCommStage& recv_stage,
   std::map<int, tbox::AsyncCommPeer<char> *>& recv_comms,
   const int* old_partition,
   const int* new_partition,
   const int num_boxes,
   const int my_rank,
   const tbox::SAMRAI_MPI& mpi) const
{
   std::set<int> recv_procs;
   for (int b = 0; b < num_boxes; ++b) {
      if (old_partition[b] == my_rank &&
          new_partition[b] != my_rank) {
         send_procs.insert(new_partition[b]);
         if (send_comms[new_partition[b]] == 0) {
            send_comms[new_partition[b]] = new tbox::AsyncCommPeer<char>();
         }
      }
      if (old_partition[b] != my_rank &&
          new_partition[b] == my_rank) {
         recv_procs.insert(old_partition[b]);
         if (recv_comms[old_partition[b]] == 0) {
            recv_comms[old_partition[b]] = new tbox::AsyncCommPeer<char>();
         }
      }
   }

   const int num_sends = static_cast<int>(send_procs.size());
   const int num_recvs = static_cast<int>(recv_procs.size());

   if (num_sends > 0) {

      for (std::set<int>::const_iterator s_itr = send_procs.begin();
           s_itr != send_procs.end(); ++s_itr) {
         int send_num = *s_itr;
         send_comms[send_num]->initialize(&send_stage);
         send_comms[send_num]->setPeerRank(*s_itr);
         send_comms[send_num]->setMPI(mpi);
         send_comms[send_num]->setMPITag(GraphLoadBalancer_LOADTAG0,
            GraphLoadBalancer_LOADTAG1);
         send_comms[send_num]->limitFirstDataLength(
            sizeof(BoxInTransit) * GraphLoadBalancer_FIRSTDATALEN);

      }
   }

   if (num_recvs > 0) {

      for (std::set<int>::const_iterator s_itr = recv_procs.begin();
           s_itr != recv_procs.end(); ++s_itr) {

         int recv_num = *s_itr;
         recv_comms[recv_num]->initialize(&recv_stage);
         recv_comms[recv_num]->setPeerRank(*s_itr);
         recv_comms[recv_num]->setMPI(mpi);
         recv_comms[recv_num]->setMPITag(GraphLoadBalancer_LOADTAG0,
            GraphLoadBalancer_LOADTAG1);
         recv_comms[recv_num]->limitFirstDataLength(
            sizeof(BoxInTransit) * GraphLoadBalancer_FIRSTDATALEN);

      }

   }

}

/*
 *  ***********************************************************************
 *  BoxInTransit stream methods
 *  ***********************************************************************
 */

void
GraphLoadBalancer::BoxInTransit::putToMessageStream(
   tbox::MessageStream& mstream) const
{
   d_box.putToMessageStream(mstream);
   d_orig_box.putToMessageStream(mstream);
}

void
GraphLoadBalancer::BoxInTransit::getFromMessageStream(
   tbox::MessageStream& mstream)
{
   d_box.getFromMessageStream(mstream);
   d_orig_box.getFromMessageStream(mstream);
}

/*
 *  ***********************************************************************
 *  Read from input
 *  ***********************************************************************
 */

void
GraphLoadBalancer::getFromInput(
   const boost::shared_ptr<tbox::Database>& input_db)
{

   if (input_db) {

      if (input_db->isInteger("target_box_size")) {
         int target_box_size[d_dim.getValue()];
         input_db->getIntegerArray("target_box_size",
            target_box_size,
            d_dim.getValue());
         d_target_box_size = hier::IntVector(d_dim, target_box_size);
         for (int i = 0; i < d_dim.getValue(); ++i) {
            if (d_target_box_size[i] <= 0) {
               INPUT_RANGE_ERROR("target_box_size");
            }
         }
      }

      d_coalesce_boxes = input_db->getBoolWithDefault("coalesce_boxes", true);

      if (input_db->isInteger("tile_size")) {
         input_db->getIntegerArray("tile_size", &d_tile_size[0], d_tile_size.getDim().getValue());
         for (int i = 0; i < d_dim.getValue(); ++i) {
            if (!(d_tile_size[i] >= 1)) {
               TBOX_ERROR("CascadePartitioner tile_size must be >= 1 in all directions.\n"
                  << "Input tile_size is " << d_tile_size);
            }
         }
      }

   }
}

}
}
