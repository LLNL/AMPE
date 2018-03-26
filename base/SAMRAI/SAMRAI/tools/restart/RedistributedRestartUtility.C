/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   $Description
 *
 ************************************************************************/

#include "RedistributedRestartUtility.h"

#ifdef HAVE_HDF5

#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <cassert>
#include <list>

#define NAME_BUF_SIZE (32)

/*
 **************************************************************************
 * writeRedistributedRestartFiles
 **************************************************************************
 */

void RedistributedRestartUtility::writeRedistributedRestartFiles(
   const string& output_dirname,
   const string& input_dirname,
   const int total_input_files,
   const int total_output_files,
   const tbox::Array<tbox::Array<int> >& file_mapping,
   const int restore_num)
{
   int num_files_written = 0;
   int num_iterations = tbox::MathUtilities<int>::Min(total_input_files,
         total_output_files);

   for (int icount = 0; icount < num_iterations; icount++) {

      //We are writing to one file or reading only one file
      int num_files_to_read = (total_input_files < total_output_files) ?
         1 : file_mapping[icount].size();
      int num_files_to_write = (total_input_files < total_output_files) ?
         file_mapping[icount].size() : 1;

      string restore_buf;
      string nodes_buf;

      restore_buf = "/restore." + tbox::Utilities::intToString(restore_num, 6);
      nodes_buf = "/nodes." + tbox::Utilities::nodeToString(total_output_files);

      string restart_dirname = output_dirname + restore_buf + nodes_buf;

      //Make the subdirectories if this is the first iteration.
      if (icount == 0) {
         tbox::Utilities::recursiveMkdir(restart_dirname);
      }

      //Mount the output files on an array of output databases
      tbox::Array<boost::shared_ptr<tbox::Database> >
      output_dbs(num_files_to_write);

      string proc_buf;
      for (int j = 0; j < num_files_to_write; j++) {

         proc_buf = "/proc." + tbox::Utilities::processorToString(
               num_files_written + j);

         string output_filename = restart_dirname + proc_buf;

         output_dbs[j].reset(new tbox::HDFDatabase(output_filename));

         int open_success = output_dbs[j]->create(output_filename);

         if (open_success < 0) {
            TBOX_ERROR(
               "Failed to open output file " << output_filename
                                             << "  HDF return code:  "
                                             << open_success);
         }

      }

      //Mount the input files on an array of input databases.
      tbox::Array<boost::shared_ptr<tbox::Database> >
      input_dbs(num_files_to_read);

      nodes_buf = "/nodes." + tbox::Utilities::nodeToString(total_input_files);

      tbox::Array<string> input_keys(0);
      tbox::Array<string> test_keys(0);
      int num_keys = 0;

      for (int i = 0; i < num_files_to_read; i++) {

         int cur_in_file_id;
         if (total_input_files < total_output_files) {
            //cur_in_file_id = num_files_written;
            cur_in_file_id = icount;
         } else {
            cur_in_file_id = file_mapping[icount][i];
         }

         proc_buf = "/proc." + tbox::Utilities::processorToString(
               cur_in_file_id);

         string restart_filename = input_dirname + restore_buf + nodes_buf
            + proc_buf;

         input_dbs[i].reset(new tbox::HDFDatabase(restart_filename));

         int open_success = input_dbs[i]->open(restart_filename);

         if (open_success < 0) {
            TBOX_ERROR(
               "Failed to open input file " << restart_filename
                                            << "  HDF return code:  "
                                            << open_success);
         }

         //Get the array of input keys.
         if (i == 0) {
            input_keys = input_dbs[i]->getAllKeys();
            num_keys = input_keys.size();
         } else {
            test_keys = input_dbs[i]->getAllKeys();
            if (test_keys.size() != num_keys) {
               TBOX_ERROR("Input files contain differing number of keys");
            }
         }
      }

      //For every input key, call the recursive function that reads from the
      //input databases and writes to output databases.
      for (int j = 0; j < num_keys; j++) {
         readAndWriteRestartData(output_dbs,
            input_dbs,
            input_keys[j],
            &file_mapping,
            num_files_written,
            icount,
            total_input_files,
            total_output_files);
      }

      //Unmount the databases.  This closes the files.
      int k;
      for (k = 0; k < num_files_to_read; k++) {
         input_dbs[k]->close();
      }
      for (k = 0; k < num_files_to_write; k++) {
         output_dbs[k]->close();
      }

      num_files_written += num_files_to_write;
   }
}

/*
 **************************************************************************
 * readAndWriteRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWriteRestartData(
   tbox::Array<boost::shared_ptr<tbox::Database> >& output_dbs,
   const tbox::Array<boost::shared_ptr<tbox::Database> >& input_dbs,
   const string& key,
   const tbox::Array<tbox::Array<int> >* file_mapping,  // = NULL
   const int num_files_written, // = -1,
   const int input_proc_num, // = -1
   const int total_input_files, // = -1,
   const int total_output_files) // = -1););
{
#ifdef DEBUG_CHECK_ASSERTIONS
   //One of the database arrays must be of size 1, and the other must be of
   //size >= 1.
   assert(output_dbs.size() >= 1);
   assert(input_dbs.size() >= 1);
   assert(input_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   //This function works under the assumption that all of the input databases
   //contain the same keys, so we only need to check the type associated with
   //the key with one input database.

   //If the key is associated with any type other than Database, then the data
   //can be read from input_dbs[0] and written to every element of output_dbs.
   //The Database case is handled separately.

   if (input_dbs[0]->isDatabase(key)) {
      boost::shared_ptr<tbox::Database> db = input_dbs[0]->getDatabase(key);

      if (db->keyExists("d_is_patch_level") &&
          db->getBool("d_is_patch_level")) {

         //Here we are handling the input database(s) for a PatchLevel.

         //Create array of level input databases.
         tbox::Array<boost::shared_ptr<tbox::Database> > level_in_dbs;
         level_in_dbs.resizeArray(input_dbs.size());

         for (int i = 0; i < input_dbs.size(); i++) {
            level_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         //input_proc_nums is an array that contains all of the processor
         //numbers that created the input databases that are currently
         //being processed.
         tbox::Array<int> input_proc_nums;
         if (total_input_files < total_output_files) {
            input_proc_nums.resizeArray(1);
            input_proc_nums[0] = input_proc_num;
         } else {
            input_proc_nums = (*file_mapping)[num_files_written];
         }

         //Call routine to write output according to the new processor mapping
         readAndWritePatchLevelRestartData(output_dbs, level_in_dbs,
            key,
            num_files_written,
            input_proc_nums,
            total_output_files);

      } else if (db->keyExists("d_is_mapped_box_level") &&
                 db->getBool("d_is_mapped_box_level")) {

         //Create array of level input databases.
         tbox::Array<boost::shared_ptr<tbox::Database> > level_in_dbs;
         level_in_dbs.resizeArray(input_dbs.size());

         for (int i = 0; i < input_dbs.size(); i++) {
            level_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         //input_proc_nums is an array that contains all of the processor
         //numbers that created the input databases that are currently
         //being processed.
         tbox::Array<int> input_proc_nums;
         if (total_input_files < total_output_files) {
            input_proc_nums.resizeArray(1);
            input_proc_nums[0] = input_proc_num;
         } else {
            input_proc_nums = (*file_mapping)[num_files_written];
         }

         readAndWriteBoxLevelRestartData(output_dbs, level_in_dbs,
            key,
            num_files_written,
            input_proc_nums,
            total_output_files);

      } else if (db->keyExists("d_is_edge_set") &&
                 db->getBool("d_is_edge_set")) {
         // don't write out edge sets.  They will be reconstructed from
         // PatchHierarchy during restarted run initialization.

         // no operation needed here.

      } else {
         //If this block is entered, then the key represents a database that
         //is not a patch level database.  We created child database arrays
         //for input_dbs and output_dbs, and then call readAndWriteRestartData
         //recursively.

         tbox::Array<boost::shared_ptr<tbox::Database> > child_in_dbs;
         child_in_dbs.resizeArray(input_dbs.size());

         for (int i = 0; i < input_dbs.size(); i++) {
            child_in_dbs[i] = input_dbs[i]->getDatabase(key);
         }

         tbox::Array<boost::shared_ptr<tbox::Database> > child_out_dbs;
         child_out_dbs.resizeArray(output_dbs.size());

         for (int i = 0; i < output_dbs.size(); i++) {
            child_out_dbs[i] = output_dbs[i]->putDatabase(key);
         }

         tbox::Array<string> child_keys = db->getAllKeys();

         for (int j = 0; j < child_keys.size(); j++) {
            readAndWriteRestartData(child_out_dbs,
               child_in_dbs,
               child_keys[j],
               file_mapping,
               num_files_written,
               input_proc_num,
               total_input_files,
               total_output_files);
         }
      }
   } else if (input_dbs[0]->isInteger(key)) {

      tbox::Array<int> int_array(0);

      int_array = input_dbs[0]->getIntegerArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putIntegerArray(key, int_array);
      }

   } else if (input_dbs[0]->isDouble(key)) {

      tbox::Array<double> double_array(0);

      double_array = input_dbs[0]->getDoubleArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putDoubleArray(key, double_array);
      }

   } else if (input_dbs[0]->isBool(key)) {

      tbox::Array<bool> bool_array(0);

      bool_array = input_dbs[0]->getBoolArray(key);

      if (key == "d_write_edges_for_restart") {
         for (int j = 0; j < bool_array.size(); j++) {
            bool_array[j] = false;
         }
      }

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putBoolArray(key, bool_array);
      }

   } else if (input_dbs[0]->isDatabaseBox(key)) {

      tbox::Array<tbox::DatabaseBox> box_array(0);

      box_array = input_dbs[0]->getDatabaseBoxArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putDatabaseBoxArray(key, box_array);
      }

   } else if (input_dbs[0]->isString(key)) {

      tbox::Array<string> string_array(0);

      string_array = input_dbs[0]->getStringArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putStringArray(key, string_array);
      }

   } else if (input_dbs[0]->isComplex(key)) {

      tbox::Array<dcomplex> complex_array(0);

      complex_array = input_dbs[0]->getComplexArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putComplexArray(key, complex_array);
      }

   } else if (input_dbs[0]->isChar(key)) {

      tbox::Array<char> char_array(0);

      char_array = input_dbs[0]->getCharArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putCharArray(key, char_array);
      }

   } else if (input_dbs[0]->isFloat(key)) {

      tbox::Array<float> float_array(0);

      float_array = input_dbs[0]->getFloatArray(key);

      for (int i = 0; i < output_dbs.size(); i++) {
         output_dbs[i]->putFloatArray(key, float_array);
      }

   } else {

      TBOX_ERROR(
         "The key " << key
                    << " is invalid or not associated with a supported datatype.");

   }
}

/*
 **************************************************************************
 * readAndWritePatchLevelRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWritePatchLevelRestartData(
   tbox::Array<boost::shared_ptr<tbox::Database> >& output_dbs,
   const tbox::Array<boost::shared_ptr<tbox::Database> >& level_in_dbs,
   const string& key,
   const int num_files_written,
   const tbox::Array<int>& input_proc_nums,
   const int total_output_files)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(output_dbs.size() >= 1);
   assert(level_in_dbs.size() >= 1);
   assert(level_in_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   //Create an array of level output databases
   tbox::Array<boost::shared_ptr<tbox::Database> > level_out_dbs;
   level_out_dbs.resizeArray(output_dbs.size());

   for (int i = 0; i < output_dbs.size(); i++) {
      level_out_dbs[i] = output_dbs[i]->putDatabase(key);
   }

   //Read in data that is global to every processor
   bool is_patch_level = level_in_dbs[0]->getBool("d_is_patch_level");
   int version = level_in_dbs[0]->getInteger("HIER_PATCH_LEVEL_VERSION");
   tbox::Array<tbox::DatabaseBox> box_array(0);
   if (level_in_dbs[0]->keyExists("d_boxes")) {
      level_in_dbs[0]->getDatabaseBoxArray("d_boxes");
   }

   tbox::Array<int> ratio_to_zero =
      level_in_dbs[0]->getIntegerArray("d_ratio_to_level_zero");
   int number_blocks = level_in_dbs[0]->getInteger("d_number_blocks");
   tbox::Array<tbox::Array<tbox::DatabaseBox> > physical_domain(number_blocks);
   for (int nb = 0; nb < number_blocks; nb++) {
      std::string domain_name = "d_physical_domain_"
         + tbox::Utilities::blockToString(nb);
      physical_domain[nb] =
         level_in_dbs[0]->getDatabaseBoxArray(domain_name);
   }
   int level_number = level_in_dbs[0]->getInteger("d_level_number");
   int next_coarser_level =
      level_in_dbs[0]->getInteger("d_next_coarser_level_number");
   bool in_hierarchy = level_in_dbs[0]->getBool("d_in_hierarchy");
   tbox::Array<int> ratio_to_coarser =
      level_in_dbs[0]->getIntegerArray("d_ratio_to_coarser_level");

   const int out_size = level_out_dbs.size();

   //Write out global data.
   for (int i = 0; i < out_size; i++) {
      level_out_dbs[i]->putBool("d_is_patch_level", is_patch_level);
      level_out_dbs[i]->putInteger("HIER_PATCH_LEVEL_VERSION", version);
      if (box_array.size() > 0) {
         level_out_dbs[i]->putDatabaseBoxArray("d_boxes", box_array);
      }
      level_out_dbs[i]->putIntegerArray("d_ratio_to_level_zero",
         ratio_to_zero);
      level_out_dbs[i]->putInteger("d_number_blocks", number_blocks);
      for (int nb = 0; nb < number_blocks; nb++) {
         std::string domain_name = "d_physical_domain_"
            + tbox::Utilities::blockToString(nb);
         level_out_dbs[i]->putDatabaseBoxArray(domain_name,
            physical_domain[nb]);
      }
      level_out_dbs[i]->putInteger("d_level_number", level_number);
      level_out_dbs[i]->putInteger("d_next_coarser_level_number",
         next_coarser_level);
      level_out_dbs[i]->putBool("d_in_hierarchy", in_hierarchy);
      level_out_dbs[i]->putIntegerArray("d_ratio_to_coarser_level",
         ratio_to_coarser);

   }

   tbox::Array<boost::shared_ptr<tbox::Database> >
   mapped_box_level_dbs_in(input_proc_nums.size());

   std::list<int> local_indices_used;
   int max_index_used = 0;

   //Each iteration of this loop processes the patches from one input
   //database.
   for (int i = 0; i < input_proc_nums.size(); i++) {

      boost::shared_ptr<tbox::Database> mbl_database =
         level_in_dbs[i]->getDatabase("mapped_box_level");

      mapped_box_level_dbs_in[i] = mbl_database;

      boost::shared_ptr<tbox::Database> mapped_boxes_db =
         mbl_database->getDatabase("mapped_boxes");

      tbox::Array<int> local_indices(0);
      if (mapped_boxes_db->keyExists("local_indices")) {
         local_indices = mapped_boxes_db->getIntegerArray("local_indices");
      }
      tbox::Array<int> block_ids(0);
      if (mapped_boxes_db->keyExists("block_ids")) {
         block_ids = mapped_boxes_db->getIntegerArray("block_ids");
      }

      //This list will contain all of the patch numbers that came from a
      //single processor.
      int mbs_size = local_indices.size();
      std::list<int> input_local_patch_nums;
      std::list<int> input_local_block_ids;
      std::list<int> output_local_patch_nums;
      std::list<int> output_local_block_ids;
      int max_local_indices = 0;
      int min_local_indices = tbox::MathUtilities<int>::getMax();

      for (int j = 0; j < mbs_size; j++) {
         bool new_patch_num = true;
         for (std::list<int>::iterator p(input_local_patch_nums.begin());
              p != input_local_patch_nums.end(); p++) {
            if (*p == local_indices[j]) {
               new_patch_num = false;
               break;
            }
         }
         if (new_patch_num) {
            input_local_patch_nums.push_front(local_indices[j]);
            input_local_block_ids.push_front(block_ids[j]);
         }
      }

      if (out_size == 1) {
         bool recompute_local_patch_nums = false;
         if (local_indices_used.size() == 0) {
            for (std::list<int>::iterator ni(input_local_patch_nums.begin());
                ni != input_local_patch_nums.end(); ni++) {
               local_indices_used.push_front(*ni);
               max_index_used =
                  tbox::MathUtilities<int>::Max(max_index_used, *ni);
            }
         } else {
            for (std::list<int>::iterator ni(input_local_patch_nums.begin());
                 ni != input_local_patch_nums.end(); ni++) {
               bool repeat_found = false;
               for (std::list<int>::iterator li(local_indices_used.begin());
                    li != local_indices_used.end(); li++) {
                  if (*ni == *li) {
                     repeat_found = true;
                     break;
                  }
               }
               if (repeat_found) {
                  recompute_local_patch_nums = true;
                  int new_value = max_index_used + 1;
                  for (int a = 0; a < mbs_size; a++) {
                     if (local_indices[a] == *ni) {
                        local_indices[a] = new_value;
                     }
                  }
                  local_indices_used.push_front(new_value);
                  max_index_used = new_value;
               } else {
                  local_indices_used.push_front(*ni);
                  max_index_used =
                     tbox::MathUtilities<int>::Max(max_index_used, *ni);
               }
            }
         }

         if (recompute_local_patch_nums) {
            for (int j = 0; j < mbs_size; j++) {
               bool new_patch_num = true;
               for (std::list<int>::iterator p(output_local_patch_nums.begin());
                    p != output_local_patch_nums.end();
                    p++) {
                  if (*p == local_indices[j]) {
                     new_patch_num = false;
                     break;
                  }
               }
               if (new_patch_num) {
                  output_local_patch_nums.push_front(local_indices[j]);
                  output_local_block_ids.push_front(block_ids[j]);
               }
            }
         } else {
            output_local_patch_nums = input_local_patch_nums;
            output_local_block_ids = input_local_block_ids;
         }
      } else {
         output_local_patch_nums = input_local_patch_nums;
         output_local_block_ids = input_local_block_ids;
      }

      for (int j = 0; j < mbs_size; j++) {
         max_local_indices = tbox::MathUtilities<int>::Max(max_local_indices,
               local_indices[j]);
         min_local_indices = tbox::MathUtilities<int>::Min(min_local_indices,
               local_indices[j]);
      }

      int indices_range;
      if (mbs_size) {
         indices_range = max_local_indices - min_local_indices;
      } else {
         indices_range = 0;
         min_local_indices = 0;
      }

      tbox::Array<int> output_dist_cutoff(out_size);
      for (int j = 0; j < out_size; j++) {
         output_dist_cutoff[j] = min_local_indices + j
            * (indices_range / out_size);
      }

      //For every patch number, get the patch database from input,
      //create a database for output, and call routine to read and
      //write patch database data.
      std::list<int>::iterator olp(output_local_patch_nums.begin());
      std::list<int>::iterator ilb(input_local_block_ids.begin());
      std::list<int>::iterator olb(output_local_block_ids.begin());
      for (std::list<int>::iterator ilp(input_local_patch_nums.begin());
           ilp != input_local_patch_nums.end(); ) {
         int output_id = 0;
         for (int a = 1; a < out_size; a++) {
            if (*olp > output_dist_cutoff[a]) {
               output_id = a;
            }
         }

         string in_patch_name;
         string out_patch_name;
         in_patch_name = "level_" + tbox::Utilities::levelToString(level_number)
            + "-patch_" + tbox::Utilities::patchToString(*ilp)
            + "-block_" + tbox::Utilities::blockToString(*ilb);
         out_patch_name = "level_" + tbox::Utilities::levelToString(
               level_number)
            + "-patch_" + tbox::Utilities::patchToString(*olp)
            + "-block_" + tbox::Utilities::blockToString(*olb);

         boost::shared_ptr<tbox::Database> patch_in_db =
            level_in_dbs[i]->getDatabase(in_patch_name);

         boost::shared_ptr<tbox::Database> patch_out_db =
            level_out_dbs[output_id]->putDatabase(out_patch_name);

         int output_proc = num_files_written + output_id;
         readAndWritePatchRestartData(patch_out_db, patch_in_db, output_proc);

         ilp++;
         olp++;
         ilb++;
         olb++;
      }
   }

   readAndWriteBoxLevelRestartData(
      level_out_dbs, mapped_box_level_dbs_in,
      "mapped_box_level", num_files_written,
      input_proc_nums, total_output_files);

}

/*
 **************************************************************************
 * readAndWritePatchRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWritePatchRestartData(
   boost::shared_ptr<tbox::Database>& patch_out_db,
   const boost::shared_ptr<tbox::Database>& patch_in_db,
   const int output_proc)
{
   //Get the keys in the patch input database.
   tbox::Array<string> keys = patch_in_db->getAllKeys();

   //Place the database on arrays of length 1.
   tbox::Array<boost::shared_ptr<tbox::Database> > in_db_array(1);
   tbox::Array<boost::shared_ptr<tbox::Database> > out_db_array(1);

   in_db_array[0] = patch_in_db;
   out_db_array[0] = patch_out_db;

   //Call recursive function to read and write the data associated with each
   //key.
   for (int i = 0; i < keys.size(); i++) {
      if (keys[i] == "d_patch_owner") {
         patch_out_db->putInteger(keys[i], output_proc);
      } else {
         readAndWriteRestartData(out_db_array, in_db_array, keys[i]);
      }
   }
}

/*
 **************************************************************************
 * readAndWriteBoxLevelRestartData
 **************************************************************************
 */

void RedistributedRestartUtility::readAndWriteBoxLevelRestartData(
   tbox::Array<boost::shared_ptr<tbox::Database> >& output_dbs,
   const tbox::Array<boost::shared_ptr<tbox::Database> >& level_in_dbs,
   const string& key,
   const int num_files_written,
   const tbox::Array<int>& input_proc_nums,
   const int total_output_files)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(output_dbs.size() >= 1);
   assert(level_in_dbs.size() >= 1);
   assert(level_in_dbs.size() == 1 || output_dbs.size() == 1);
#endif

   const int out_size = output_dbs.size();

   //Create an array of level output databases
   tbox::Array<boost::shared_ptr<tbox::Database> > level_out_dbs;
   level_out_dbs.resizeArray(out_size);

   for (int i = 0; i < out_size; i++) {
      level_out_dbs[i] = output_dbs[i]->putDatabase(key);
   }

   bool is_mapped_box_level = level_in_dbs[0]->getBool("d_is_mapped_box_level");
   int version = level_in_dbs[0]->getInteger("HIER_MAPPED_BOX_LEVEL_VERSION");
   int dim = level_in_dbs[0]->getInteger("dim");
   tbox::Array<int> ratio = level_in_dbs[0]->getIntegerArray("d_ratio");

   for (int i = 0; i < out_size; i++) {
      level_out_dbs[i]->putBool("d_is_mapped_box_level", is_mapped_box_level);
      level_out_dbs[i]->putInteger("HIER_MAPPED_BOX_LEVEL_VERSION", version);
      level_out_dbs[i]->putInteger("dim", dim);
      level_out_dbs[i]->putIntegerArray("d_ratio", ratio);

      level_out_dbs[i]->putInteger("d_nproc", total_output_files);
      level_out_dbs[i]->putInteger("d_rank", num_files_written + i);
   }

   std::vector<int>* out_local_indices = new std::vector<int>[out_size];
   std::vector<int>* out_ranks = new std::vector<int>[out_size];
   std::vector<int>* out_periodic_ids =
      new std::vector<int>[out_size];

   tbox::Array<tbox::DatabaseBox>* out_box_array =
      new tbox::Array<tbox::DatabaseBox>[out_size];

   int out_vec_size = 0;
   tbox::Array<int> out_mbs_size(out_size, 0);

   std::list<int> local_indices_used;
   int max_index_used = 0;

   //Each iteration of this loop processes the patches from one input
   //database.
   for (int i = 0; i < input_proc_nums.size(); i++) {

      boost::shared_ptr<tbox::Database> mapped_boxes_in_db =
         level_in_dbs[i]->getDatabase("mapped_boxes");

      int mbs_size = mapped_boxes_in_db->getInteger("mapped_box_set_size");

      tbox::Array<int> local_indices;
      tbox::Array<int> ranks;
      tbox::Array<int> periodic_ids;
      tbox::Array<tbox::DatabaseBox> boxes;

      if (mapped_boxes_in_db->keyExists("local_indices")) {
         local_indices = mapped_boxes_in_db->getIntegerArray("local_indices");
      }
      if (mapped_boxes_in_db->keyExists("ranks")) {
         ranks = mapped_boxes_in_db->getIntegerArray("ranks");
      }
      if (mapped_boxes_in_db->keyExists("periodic_ids")) {
         periodic_ids =
            mapped_boxes_in_db->getIntegerArray("periodic_ids");
      }
      if (mapped_boxes_in_db->keyExists("boxes")) {
         boxes = mapped_boxes_in_db->getDatabaseBoxArray("boxes");
      }

      if (out_size == 1) {
         std::list<int> new_indices;
         for (int k = 0; k < mbs_size; k++) {
            bool is_new_index = true;
            for (std::list<int>::iterator ni(new_indices.begin());
                 ni != new_indices.end(); ni++) {
               if (local_indices[k] == *ni) {
                  is_new_index = false;
                  break;
               }
            }
            if (is_new_index) {
               new_indices.push_front(local_indices[k]);
            }
         }
         if (local_indices_used.size() == 0) {
            for (std::list<int>::iterator ni(new_indices.begin());
                 ni != new_indices.end(); ni++) {
               local_indices_used.push_front(*ni);
               max_index_used =
                  tbox::MathUtilities<int>::Max(max_index_used, *ni);
            }
         } else {
            for (std::list<int>::iterator ni(new_indices.begin());
                 ni != new_indices.end(); ni++) {
               bool repeat_found = false;
               for (std::list<int>::iterator li(local_indices_used.begin());
                    li != local_indices_used.end();
                    li++) {
                  if (*ni == *li) {
                     repeat_found = true;
                     break;
                  }
               }
               if (repeat_found) {
                  int new_value = max_index_used + 1;
                  for (int a = 0; a < mbs_size; a++) {
                     if (local_indices[a] == *ni) {
                        local_indices[a] = new_value;
                     }
                  }
                  local_indices_used.push_front(new_value);
                  max_index_used = new_value;
               } else {
                  local_indices_used.push_front(*ni);
                  max_index_used =
                     tbox::MathUtilities<int>::Max(max_index_used, *ni);
               }
            }
         }
      }

      int max_local_indices = 0;
      int min_local_indices = tbox::MathUtilities<int>::getMax();
      for (int j = 0; j < mbs_size; j++) {
         max_local_indices = tbox::MathUtilities<int>::Max(max_local_indices,
               local_indices[j]);
         min_local_indices = tbox::MathUtilities<int>::Min(min_local_indices,
               local_indices[j]);
      }

      int indices_range = max_local_indices - min_local_indices;

      tbox::Array<int> output_dist_cutoff(out_size);
      for (int j = 0; j < out_size; j++) {
         output_dist_cutoff[j] = min_local_indices + j
            * (indices_range / out_size);
      }

      if (mbs_size > 0) {
         out_vec_size += mbs_size;

         tbox::Array<int> output_ids(mbs_size);

         for (int k = 0; k < mbs_size; k++) {
            output_ids[k] = 0;
            for (int a = 1; a < out_size; a++) {
               if (local_indices[k] > output_dist_cutoff[a]) {
                  output_ids[k] = a;
               }
            }
         }

         for (int j = 0; j < out_size; j++) {
            out_local_indices[j].reserve(out_vec_size);
            out_ranks[j].reserve(out_vec_size);
            out_periodic_ids[j].reserve(out_vec_size);

            out_box_array[j].resizeArray(out_vec_size);
            for (int k = 0; k < mbs_size; k++) {
               if (output_ids[k] == j) {
                  int output_rank = num_files_written + output_ids[k];
                  out_local_indices[j].push_back(local_indices[k]);
                  out_ranks[j].push_back(output_rank);
                  out_periodic_ids[j].push_back(
                     periodic_ids[k]);
                  out_box_array[j][out_mbs_size[j]++] = boxes[k];
               }
            }
         }
      }
   }

   for (int j = 0; j < out_size; j++) {
      boost::shared_ptr<tbox::Database> mapped_boxes_out_db =
         level_out_dbs[j]->putDatabase("mapped_boxes");

      mapped_boxes_out_db->
      putInteger("HIER_MAPPED_BOX_SET_VERSION", version);
      mapped_boxes_out_db->
      putInteger("mapped_box_set_size", out_mbs_size[j]);

      if (out_mbs_size[j]) {
         mapped_boxes_out_db->
         putIntegerArray("local_indices",
            &out_local_indices[j][0],
            out_mbs_size[j]);
         mapped_boxes_out_db->
         putIntegerArray("periodic_ids",
            &out_periodic_ids[j][0],
            out_mbs_size[j]);
         mapped_boxes_out_db->
         putIntegerArray("ranks", &out_ranks[j][0], out_mbs_size[j]);
         mapped_boxes_out_db->
         putDatabaseBoxArray("boxes", &out_box_array[j][0], out_mbs_size[j]);

      }
   }

   delete[] out_local_indices;
   delete[] out_ranks;
   delete[] out_periodic_ids;
   delete[] out_box_array;
/*
 *    for (int j = 0; j < out_size; j++) {
 *       std::vector<int> out_local_indices;
 *       std::vector<int> out_ranks;
 *       std::vector<int> out_periodic_ids;
 *       out_local_indices.reserve(mbs_size);
 *       out_ranks.reserve(mbs_size);
 *       out_periodic_ids.reserve(mbs_size);
 *
 *       tbox::Array<tbox::DatabaseBox> out_box_array(mbs_size);
 *
 *       int out_mbs_size = 0;
 *       for (int k = 0; k < mbs_size; k++) {
 *          int output_id = 0;
 *          for (int a = 1; a < out_size; a++) {
 *             if (local_indices[k] > output_dist_cutoff[a]) {
 *                output_id = a;
 *             }
 *          }
 *          if (output_id == j) {
 *             int output_rank = num_files_written + output_id;
 *             out_local_indices.push_back(local_indices[k]);
 *             out_ranks.push_back(output_rank);
 *             out_periodic_ids.push_back(periodic_ids[k]);
 *             out_box_array[out_mbs_size++] = boxes[k];
 *          }
 *       }
 *
 *       boost::shared_ptr<tbox::Database> mapped_boxes_out_db =
 *          level_out_dbs[j]->putDatabase("mapped_boxes");
 *
 *       mapped_boxes_out_db->
 *          putInteger("HIER_MAPPED_BOX_SET_VERSION", version);
 *       mapped_boxes_out_db->
 *          putInteger("mapped_box_set_size", out_mbs_size);
 *
 *       if (out_mbs_size) {
 *          mapped_boxes_out_db->
 *             putIntegerArray("local_indices",
 *                             &out_local_indices[0],
 *                             out_mbs_size);
 *          mapped_boxes_out_db->
 *             putIntegerArray("periodic_ids",
 *                             &out_periodic_ids[0],
 *                             out_mbs_size);
 *          mapped_boxes_out_db->
 *             putIntegerArray("ranks", &out_ranks[0], out_mbs_size);
 *          mapped_boxes_out_db->
 *             putDatabaseBoxArray("boxes", &out_box_array[0], out_mbs_size);
 *
 *       }
 *    }
 */
}

#endif
