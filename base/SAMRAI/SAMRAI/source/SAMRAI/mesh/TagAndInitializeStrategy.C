/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for params, tagging, init for gridding.
 *
 ************************************************************************/

#ifndef included_mesh_TagAndInitializeStrategy_C
#define included_mesh_TagAndInitializeStrategy_C

#include "SAMRAI/mesh/TagAndInitializeStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

#include <stdio.h>

namespace SAMRAI {
namespace mesh {

TagAndInitializeStrategy::TagAndInitializeStrategy(
   const tbox::Dimension& dim,
   const std::string& object_name):
   d_dim(dim),
   d_object_name(object_name)
{
}

TagAndInitializeStrategy::~TagAndInitializeStrategy()
{
}

/*
 *************************************************************************
 *
 * Sets refine boxes for case where refine region is specified by the
 * user.  The bool return value specifies whether or not the refine
 * boxes have been reset from the last time the method was called
 * (true = they are reset, false = they have NOT changed).
 *
 * Note that if any method which invokes tagging is performed there
 * is always potential that the boxes have changed so this method will
 * always return true in this case.
 *
 *************************************************************************
 */
bool
TagAndInitializeStrategy::getUserSuppliedRefineBoxes(
   hier::BoxContainer& refine_boxes,
   const int level_num,
   const double time)
{
   TBOX_ASSERT(level_num >= 0);
   TBOX_ASSERT(time >= 0.);

   /*
    * The cycle counter and boolean array specifying whether times
    * are used are initialially set based on inputs. There could be
    * circumstances where a set of refine boxes is requested for a
    * level number greater than the number of entries the user
    * supplied in input.  If this occurs, resize the cycle counter
    * and boolean times arrays and set appropriate defaults to
    * avoid logic errors elsewhere.
    */
   if (level_num >= d_refine_boxes_cycle_counter.getSize()) {

      d_refine_boxes_times.resizeArray(level_num + 1);
      d_refine_boxes_cycles.resizeArray(level_num + 1);
      d_refine_boxes_use_times.resizeArray(level_num + 1);
      d_refine_boxes_cycle_counter.resizeArray(level_num + 1);
      d_refine_boxes_old_seq_num.resizeArray(level_num + 1);

      d_refine_boxes_times[level_num].resizeArray(1);
      d_refine_boxes_cycles[level_num].resizeArray(1);
      d_refine_boxes_times[level_num][0] = 0.;
      d_refine_boxes_cycles[level_num][0] = 0;
      d_refine_boxes_use_times[level_num] = false;
      d_refine_boxes_cycle_counter[level_num] = 0;
      d_refine_boxes_old_seq_num[level_num] = -1;

   }

   /*
    * Increment step counter.
    */
   d_refine_boxes_cycle_counter[level_num]++;

   /*
    * Determine which sequence entry in the refine_boxes box array
    * to use.
    */
   int seq_num = 0;
   if (d_refine_boxes_use_times[level_num]) {

      /*
       * If we are using times, the user has supplied an array
       * times = {time1, time2, ...} that corresponds with the
       * refine box array boxes = {boxarr1, boxarr2, ...}. Pick
       * the appropriate sequence number based on the specified
       * time.
       */
      for (int i = 0; i < d_refine_boxes_times[level_num].getSize();
           i++) {
         if (time > d_refine_boxes_times[level_num][i]) seq_num = i;
      }
   } else {

      /*
       * If we are using steps, the user has supplied an array
       * cycles = {cycle1, cycle2, ...} that corresponds with the
       * refine box array boxes = {boxarr1, boxarr2, ...}.  Pick
       * the appropriate seq number based on the counter.
       */
      for (int i = 0; i < d_refine_boxes_cycles[level_num].getSize();
           i++) {
         if (d_refine_boxes_cycle_counter[level_num] >
             d_refine_boxes_cycles[level_num][i])
            seq_num = i;
      }
   }

   /*
    * Print some warnings if the user-supplied entries will not
    * generate any refined boxes for the level.
    */
   hier::BoxContainer empty_boxes;
   if ((d_refine_boxes.getSize() <= level_num) ||
       (d_refine_boxes[level_num][seq_num].size() == 0)) {

      TBOX_WARNING(
         d_object_name << ": getRefineBoxes\n"
                       << "No refine boxes specified for level "
                       << level_num);
      refine_boxes = empty_boxes;

   } else if (d_refine_boxes[level_num].getSize() <= seq_num) {

      if (d_refine_boxes_use_times[level_num]) {
         TBOX_WARNING(
            d_object_name << ": getRefineBoxes\n"
                          << "No refine boxes specified for time sequence "
                          << seq_num << " on level " << level_num
                          << ".\n No refinement will be performed.");
      } else {
         TBOX_WARNING(
            d_object_name << ": getRefineBoxes\n"
                          << "No refine boxes specified for step sequence "
                          << seq_num << " on level " << level_num
                          << ".\n No refinement will be performed.");
      }
      refine_boxes = empty_boxes;

   } else {

      refine_boxes = d_refine_boxes[level_num][seq_num];

   }

   /*
    * If the user has requested their own particular set of refine
    * boxes (i.e. by calling resetRefineBoxes()), overwrite any previously
    * determined refine boxes with their requested set.
    */
   bool use_reset = false;
   if (d_refine_boxes_reset.getSize() > level_num) {
      if (d_refine_boxes_reset[level_num]) {
         use_reset = true;
         refine_boxes = d_reset_refine_boxes[level_num];
      }
   }

   /*
    * If we have not moved to a new sequence number, or otherwise reset
    * boxes from the last time this method was called, then we return
    * "false", indicating boxes have NOT been reset.  If we have reset
    * the boxes, return "true".
    */
   bool modified_refine_boxes;
   if (!use_reset && d_refine_boxes_old_seq_num[level_num] == seq_num) {
      modified_refine_boxes = false;
   } else {
      d_refine_boxes_old_seq_num[level_num] = seq_num;
      modified_refine_boxes = true;
   }

   /*
    * If one of the tagging methods (e.g. gradient detector or
    * Richardson extrapolation) is used, the boxes may be modified
    * even if they have not been modified by user input.
    */
   if (!refineUserBoxInputOnly()) modified_refine_boxes = true;

   return modified_refine_boxes;

}

/*
 *************************************************************************
 *
 * Resets refine boxes for specified level.
 *
 *************************************************************************
 */

void
TagAndInitializeStrategy::resetRefineBoxes(
   const hier::BoxContainer& refine_boxes,
   const int level_num)
{
   TBOX_ASSERT(level_num >= 0);

   int i = d_reset_refine_boxes.getSize();
   if (i <= level_num) {
      d_reset_refine_boxes.resizeArray(level_num + 1);
      d_refine_boxes_reset.resizeArray(level_num + 1);
      for ( ; i < d_reset_refine_boxes.getSize(); ++i) {
         d_refine_boxes_reset[i] = false;
      }
   }

   d_refine_boxes_reset[level_num] = true;
   d_reset_refine_boxes[level_num] = refine_boxes;

}

/*
 *************************************************************************
 *
 * Read specified refinement boxes, if any.
 *
 *************************************************************************
 */

void
TagAndInitializeStrategy::getFromInput(
   const boost::shared_ptr<tbox::Database>& db)
{
   TBOX_ASSERT(db);

   /*
    * Read refine boxes.
    */
   if (!db->keyExists("RefineBoxes")) {
      TBOX_ERROR("If REFINE_BOXES is used as a tagging_method, you must\n"
         << "provide a `RefineBoxes' database entry in the input\n"
         << "to specify the level boxes to be refined.  See header\n"
         << "for TagAndInitializeStrategy class for \n"
         << "discussion of the entry format." << std::endl);
   }
   boost::shared_ptr<tbox::Database> refine_box_db(
      db->getDatabase("RefineBoxes"));
   tbox::Array<std::string> box_keys = refine_box_db->getAllKeys();
   int nkeys = box_keys.getSize();

   d_refine_boxes.resizeArray(nkeys);
   d_refine_boxes_cycles.resizeArray(nkeys);
   d_refine_boxes_times.resizeArray(nkeys);
   d_refine_boxes_use_times.resizeArray(nkeys);
   d_refine_boxes_cycle_counter.resizeArray(nkeys);
   d_refine_boxes_old_seq_num.resizeArray(nkeys);

   /*
    * Do a check here to see if we are using the "old" version of input
    * (i.e. v1.3.1 or before) which only allows you to specify a single
    * set of refine boxes, or the "new" input format which allows you
    * to specify a specified sequence of refine boxes.
    */
   bool use_new_input = false;
   for (int ln = 0; ln < nkeys; ln++) {
      std::string level_boxes_name = "level_" + tbox::Utilities::intToString(ln);
      if (refine_box_db->isDatabase(level_boxes_name)) {
         use_new_input = true;
      }
   }

   for (int ln = 0; ln < nkeys; ln++) {
      /*
       * Set counter for each level to zero.
       */
      d_refine_boxes_cycle_counter[ln] = 0;

      /*
       * Set old sequence number to -1, assuring it will always
       * interpret a new sequence at initialization, when the
       * counter is zero.
       */
      d_refine_boxes_old_seq_num[ln] = -1;
   }

   if (use_new_input) {
      /*
       * We are using the updated input format that allows multiple
       * sequence entries.
       */
      for (int ln = 0; ln < nkeys; ln++) {
         std::string level_boxes_name = "level_" + tbox::Utilities::intToString(
               ln);
         if (!refine_box_db->keyExists(level_boxes_name)) {
            TBOX_ERROR(
               d_object_name << "\n"
                             << ": Expected sub-database level entries in the\n"
                             << " 'RefineBoxes' database to specify boxes for\n"
                             << "  different time or cycle sequences: \n"
                             << "  e.g.  Level0 { boxes_0 = <box array> \n"
                             << "                 boxes_1 = <box array> \n"
                             << "                 ...} \n"
                             << "        Level1 { boxes_0 = <box array> \n"
                             << "                 ...} \n"
                             << "See header for this class for further discussion\n"
                             << "of the expected input format."
                             << std::endl);
         }
         boost::shared_ptr<tbox::Database> level_refine_box_db(
            refine_box_db->getDatabase(level_boxes_name));

         /*
          * Read cycles.
          */
         bool use_cycles = false;
         if (level_refine_box_db->keyExists("cycles")) {
            d_refine_boxes_cycles[ln] =
               level_refine_box_db->getIntegerArray("cycles");
            use_cycles = true;
         } else {
            d_refine_boxes_cycles[ln].resizeArray(1);
            d_refine_boxes_cycles[ln][0] = 0;
         }

         /*
          * Read times.
          */
         if (level_refine_box_db->keyExists("times")) {
            d_refine_boxes_use_times[ln] = true;
            d_refine_boxes_times[ln] =
               level_refine_box_db->getDoubleArray("times");
            if (use_cycles) {
               TBOX_WARNING(
                  d_object_name << ": You have entries for "
                                << "both 'cycles' and 'times' for level " << ln
                                << ".\n Because 'times' takes precedence, the "
                                << "'cycles' entries will be ignored." << std::endl);
            }
         } else {
            d_refine_boxes_times[ln].resizeArray(1);
            d_refine_boxes_times[ln][0] = 0.;
            d_refine_boxes_use_times[ln] = false;
         }

         /*
          * Use the size of the "cycles" or "times" arrays to govern
          * how many boxes entries there may be.
          */
         int max_seq = 1;
         if (d_refine_boxes_use_times[ln]) {
            max_seq = d_refine_boxes_times[ln].getSize();
         } else {
            max_seq = d_refine_boxes_cycles[ln].getSize();
         }
         d_refine_boxes[ln].resizeArray(max_seq);

         /*
          * Read boxes.
          */
         for (int i = 0; i < max_seq; i++) {
            std::string boxes_name = "boxes_" + tbox::Utilities::intToString(i);
            if (level_refine_box_db->keyExists(boxes_name)) {
               d_refine_boxes[ln][i] =
                  level_refine_box_db->getDatabaseBoxArray(boxes_name);
            }
         }

      } // loop over levels
   } else {

      for (int ln = 0; ln < nkeys; ln++) {
         std::string level_boxes_name = "level_" + tbox::Utilities::intToString(
               ln);
         d_refine_boxes[ln].resizeArray(1);
         d_refine_boxes_cycles[ln].resizeArray(1);
         d_refine_boxes_times[ln].resizeArray(1);

         if (refine_box_db->keyExists(level_boxes_name)) {
            d_refine_boxes[ln][0] =
               refine_box_db->getDatabaseBoxArray(level_boxes_name);
         }
         d_refine_boxes_cycles[ln][0] = 0;
         d_refine_boxes_times[ln][0] = 0.;
         d_refine_boxes_use_times[ln] = false;
      }
   } // not using new input format
}

}
}
#endif
