/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Utility for building efficient communication tree.
 *
 ************************************************************************/
#ifndef included_tbox_BalancedDepthFirstTree
#define included_tbox_BalancedDepthFirstTree

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Utility to compute neighbors in a balanced depth-first tree.
 *
 * This class is a tool to create a processor tree appropriate for use
 * in collective communication operations.  An example of a tree created
 * is
 * @verbatim
 *                     0
 *                    / \
 *                   /   \
 *                  /     \
 *                 1       8
 *                / \     / \
 *               2   5   9   12
 *              /|  /|  / \  |\
 *             3 4 6 7 10 11 13...
 * @endverbatim
 *
 * The tree is as balanced as possible.  Nodes that are close together
 * in the tree tends to be close together in natural ordering.  Without
 * knowing about the underlying message passing network structure, we
 * assume that close natural ordering usually means close together on
 * the network.  Thus nodes close together in the tree are also close
 * together on the network.  Thus, communication between nearest
 * neighbors in the tree tend to be faster.
 *
 * The tree formed by this class has the property that each and every
 * subtree is composed nodes with contiguous natural ordering.  This
 * again benefits communication.
 */
class BalancedDepthFirstTree
{

public:
   typedef int LocalId;

   /*!
    * @brief Constructor.
    */
   BalancedDepthFirstTree();

   /*!
    * @brief Initializing constructor.
    *
    * @see initialize().
    */
   BalancedDepthFirstTree(
      unsigned int first_rank,
      unsigned int last_rank,
      unsigned int rank,
      bool do_left_leaf_switch);

   /*!
    * @brief Destructor.
    *
    * Deallocate internal data.
    */
   ~BalancedDepthFirstTree();

   /*!
    * @brief Construct the tree.
    *
    * Initialize the parent, left and right children, given a
    * contiguous range of processor ranks and the desired rank.
    * Initialization has complexity ln(length of range).
    *
    * The tree is set up for collective communications with the root
    * rank being the first_rank.
    *
    * @param first_rank The first in a contiguous range of ranks in the
    * communication group.
    *
    * @param last_rank The last in a contiguous range of ranks in the
    * communication group.
    *
    * @param rank The rank whose parent and children are sought.
    *
    * @param do_left_leaf_switch
    */
   void
   initialize(
      unsigned int first_rank,
      unsigned int last_rank,
      unsigned int rank,
      bool do_left_leaf_switch);

   /*!
    * @brief Access the rank used to initialize.
    */
   unsigned int
   getRank() const
   {
      return d_rank;
   }

   /*!
    * @brief Access the parent rank.
    */
   unsigned int
   getParentRank() const
   {
      return d_parent;
   }

   /*!
    * @brief Access a child rank.
    *
    * Currently, child_number must be 0 or 1 because we only support
    * a binary tree.
    */
   unsigned int
   getChildRank(
      unsigned int child_number) const
   {
      return d_children[child_number];
   }

   unsigned short int
   getNumberOfChildren() const
   {
      return d_num_children;
   }

   /*!
    * @brief What this class considers an invalid rank.
    *
    * When a parent or child does not exist, this value is returned for
    * the rank.
    */
   unsigned int
   getInvalidRank() const
   {
      return 1 << (8 * sizeof(unsigned int) - 2);
   }

private:
   /*!
    * @brief Initialized rank.
    *
    * @see initialize();
    */
   unsigned int d_rank;

   unsigned int d_parent;

   unsigned int d_children[2];

   unsigned short int d_num_children;

};

}
}

#endif  // included_tbox_BalancedDepthFirstTree
