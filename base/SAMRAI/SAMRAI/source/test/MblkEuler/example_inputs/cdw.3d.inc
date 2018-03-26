/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

Multiblock {

   num_blocks = 4

   // Geometries of the hierachy blocks
   BlockGeometry0 {
       domain_boxes = [ (0,0,0) , (2,2,2) ]
   }

   BlockGeometry1 {
       domain_boxes = [ (0,0,0) , (2,2,7) ]
   }

   BlockGeometry2 {
       domain_boxes = [ (0,0,0) , (2,7,2) ]
   }

   BlockGeometry3 {
       domain_boxes = [ (0,0,0) , (7,2,2) ]
   }

   // Singularity edge definition
   Singularity0 {
       blocks = 0, 1, 2
       sing_box_0 = [ (0,3,3), (2,3,3) ]
       sing_box_1 = [ (0,3,-1), (2,3,-1) ]
       sing_box_2 = [ (0,-1,3), (2,-1,3) ]
   }

   Singularity1 {
       blocks = 0, 1, 3
       sing_box_0 = [ (3,0,3), (3,2,3) ]
       sing_box_1 = [ (3,0,-1), (3,2,-1) ]
       sing_box_3 = [ (-1,0,3), (-1,2,3) ]
   }

   Singularity2 {
       blocks = 0, 2, 3
       sing_box_0 = [ (3,3,0), (3,3,2) ]
       sing_box_2 = [ (3,-1,0), (3,-1,2) ]
       sing_box_3 = [ (-1,3,0), (-1,3,2) ]
   }

   Singularity3 {
       blocks = 1, 2, 3
       sing_box_1 = [ (3,3,0), (3,3,7) ]
       sing_box_2 = [ (3,0,3), (3,7,3) ]
       sing_box_3 = [ (0,3,3), (7,3,3) ]
   }

   // define face connectivity
   BlockNeighbors0 {
       block_a = 0
       block_b = 1
       rotation_b_to_a = "I_UP", "J_UP", "K_UP"
       point_in_a_space = 0, 0, 0
       point_in_b_space = 0, 0, -3
   }

   BlockNeighbors1 {
       block_a = 0
       block_b = 2
       rotation_b_to_a = "I_UP", "J_UP", "K_UP"
       point_in_a_space = 0, 0, 0
       point_in_b_space = 0, -3, 0
   }

   BlockNeighbors2 {
       block_a = 0
       block_b = 3
       rotation_b_to_a = "I_UP", "J_UP", "K_UP"
       point_in_a_space = 0, 0, 0
       point_in_b_space = -3, 0, 0
   }

   BlockNeighbors3 {
       block_a = 1
       block_b = 2
       rotation_b_to_a = "I_UP", "K_DOWN", "J_UP"
       point_in_a_space = 0, 0, 0
       point_in_b_space = 0, 0, 5
   }

   BlockNeighbors4 {
       block_a = 1
       block_b = 3
       rotation_b_to_a = "K_DOWN", "J_UP", "I_UP"
       point_in_a_space = 0, 0, 0
       point_in_b_space = 0, 0, 5
   }

   BlockNeighbors5 {
       block_a = 2
       block_b = 3
       rotation_b_to_a = "J_DOWN", "I_UP", "K_UP"
       point_in_a_space = 0, 0, 0
       point_in_b_space = 0, 5, 0
   }

}
