/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#ifndef included_tbox_BalancedDepthFirstTree_C
#define included_tbox_BalancedDepthFirstTree_C

#include "SAMRAI/tbox/BalancedDepthFirstTree.h"

#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

BalancedDepthFirstTree::BalancedDepthFirstTree()
{
}

BalancedDepthFirstTree::BalancedDepthFirstTree(
   unsigned int first_rank,
   unsigned int last_rank,
   unsigned int my_rank,
   bool do_left_leaf_switch)
{
   initialize(first_rank, last_rank, my_rank, do_left_leaf_switch);
}

BalancedDepthFirstTree::~BalancedDepthFirstTree()
{
}

void
BalancedDepthFirstTree::initialize(
   unsigned int first_rank,
   unsigned int last_rank,
   unsigned int rank,
   bool do_left_leaf_switch)
{
   TBOX_ASSERT(first_rank <= rank);
   TBOX_ASSERT(rank <= last_rank);
#if defined(BalancedDepthFirstTree_ExtraDebug)
   plog
   << "BalancedDepthFirstTree::initialize with first_rank,last_rank,x="
   << first_rank << " " << last_rank << " " << rank << std::endl;
#endif

   unsigned int rbeg = first_rank;
   unsigned int rend = last_rank;
   unsigned int up = getInvalidRank();          // Temporary guess for parent.
   unsigned int upp = getInvalidRank(); // Temporary guess for grandparent.
   unsigned int cl, cr;         // Temporary guesses for children.
   bool is_switchable = false;    // Part of a left-leaf switchable trio.

   size_t nr;           // Number of nodes on right branch
   size_t nl;           // Number of nodes on left branch

   while (1) {

      /*
       * Walk from root to leaf to find the position of rank, its
       * parent and its children.
       */

      unsigned int node = rbeg; // Node being examined
      size_t nrem = rend - rbeg;  // Number or nodes remaining, excluding node.

      nr = nrem / 2;      // Number on right branch
      nl = nrem - nr;     // Number on left branch

      /*
       * Both children are leaves => parent and children make
       * a switchable trio.
       */
      if (nrem == 2) {
         is_switchable = true;
      }

      cl = getInvalidRank();
      cr = getInvalidRank();
      if (nl > 0) cl = node + 1;        // left child
      if (nr > 0) cr = cl + static_cast<int>(nl);         // right child

#if defined(BalancedDepthFirstTree_ExtraDebug)
      plog << "There are" << " " << nrem << " "
           << "remaining nodes.  nl,nr=" << " " << nl << " " << nr
           << std::endl;
      plog << "cl=" << " " << cl << " " << "cr=" << " " << cr
           << std::endl;
#endif

      if (node == rank) break;
      else {
         TBOX_ASSERT(nl > 0);
         TBOX_ASSERT(cl != getInvalidRank());
         upp = up;
         up = node;
         if (nr < 1 || rank < cr) {
#if defined(BalancedDepthFirstTree_ExtraDebug)
            plog << "Going left to" << " " << cl << std::endl;
#endif
            rbeg = cl;
            rend = cl + static_cast<int>(nl) - 1;
         } else {
            TBOX_ASSERT(nr > 0);
            TBOX_ASSERT(cr != getInvalidRank());
#if defined(BalancedDepthFirstTree_ExtraDebug)
            plog << "Going right to" << " " << cr << std::endl;
#endif
            rbeg = cr;
            rend = cr + static_cast<int>(nr) - 1;
         }
      }
   }

   const int gparent = upp;
   d_rank = rank;
   d_parent = up;
   d_children[0] = cl;
   d_children[1] = cr;
   d_num_children = 0;
#if defined(BalancedDepthFirstTree_ExtraDebug)
   plog << "  " << d_parent << std::endl;
   plog << "   |" << std::endl;
   plog << "  " << d_rank << std::endl;
   plog << "  /  \\ " << std::endl;
   plog << d_children[0] << "   " << d_children[1] << std::endl;
#endif

   if (do_left_leaf_switch) {
      if (is_switchable) {
         /*
          * Trios of a parent and 2 leaf children are subject to
          * switching, in which the parent and left child switch
          * places:
          *
          * Before:    parent              After:   left
          *            /  \                        /  \
          *        left    right             parent    right
          */
         if (nl == 1) {
            // Parent in a left-leaf switchable.
            d_parent = cl;
            d_children[0] = getInvalidRank();
            d_children[1] = getInvalidRank();
         } else if (rank == d_parent + 1) {
            // Left child in a left-leaf switchable.
            d_children[0] = d_parent;
            d_parent = gparent;
            d_children[1] = rank + 1;
         } else {
            // Right child in a left-leaf switchable.
            d_parent = d_parent + 1;
         }
      } else {
         /*
          * Rank is not in a switchable trio, but its children
          * may be.  Example:
          *
          * Before:       rank                   After:       rank
          *             /    \                              /    \
          *            /      \                            /      \
          *      rank+1        rank+4                rank+2        rank+5
          *     /      \      /      \              /      \      /      \
          *    /        \    /        \            /        \    /        \
          * rank+2   rank+3  rank+5   rank+6    rank+1   rank+3  rank+4   rank+6
          */
         if (nl == 3) {
            ++d_children[0];
         }
         if (nr == 3) {
            ++d_children[1];
         }
      }
   }

   for (int i = 0; i < 2; ++i) {
      if (d_children[i] != getInvalidRank()) {
         ++d_num_children;
      }
   }

#if defined(BalancedDepthFirstTree_ExtraDebug)
   plog << "  " << d_parent << std::endl;
   plog << "   |" << std::endl;
   plog << "  " << d_rank << std::endl;
   plog << "  /  \\ " << std::endl;
   plog << d_children[0] << "   " << d_children[1] << std::endl;
#endif
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
