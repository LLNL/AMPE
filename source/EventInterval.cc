// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
//
#include "EventInterval.h"

#include "SAMRAI/tbox/Utilities.h"


EventInterval::EventInterval(std::shared_ptr<tbox::Database> input_db,
                             const std::string name, const double default_value,
                             const std::string default_type,
                             const bool include_first, const bool include_last)
    : d_name(name)
{
   std::string interval_type = default_type;
   d_dbl_interval_value = default_value;
   d_int_interval_value = int(default_value + 0.5);
   d_include_first_step = include_first;
   d_include_last_step = include_last;

   std::shared_ptr<tbox::Database> interval_db;

   if (input_db) {
      if (input_db->isDatabase(name)) {
         interval_db = input_db->getDatabase(name);

         if (interval_db->keyExists("interval_type")) {
            interval_type = interval_db->getString("interval_type");
         }

         if (interval_db->keyExists("include_first_step")) {
            d_include_first_step = interval_db->getBool("include_first_step");
         }

         if (interval_db->keyExists("include_last_step")) {
            d_include_last_step = interval_db->getBool("include_last_step");
         }
      }
   }

   if (interval_type != "step" && interval_type != "time" &&
       interval_type != "dt" && interval_type != "cycle") {

      TBOX_ERROR("Error in EventInterval \"" << name << "\":\n"
                                             << "  invalid interval_type");
   }

   if (interval_type == "step" || interval_type == "cycle") {

      d_interval_type = CYCLE;
      d_dbl_interval_value = 0.0;
      if (interval_db) {
         if (interval_db->keyExists("interval")) {
            d_int_interval_value = interval_db->getInteger("interval");
         }
      }

   } else {

      d_interval_type = TIME;
      d_int_interval_value = 0;
      if (interval_db) {
         if (interval_db->keyExists("interval")) {
            d_dbl_interval_value = interval_db->getDouble("interval");
         }
      }
   }

   d_previous_time = -1.0;
   d_time_at_last_event = -1.0;
}

bool EventInterval::isActive()
{
   bool active = false;
   if (d_interval_type == CYCLE && d_int_interval_value > 0) {
      active = true;
   } else if (d_interval_type == TIME && d_dbl_interval_value > 0.) {
      active = true;
   }
   return active;
}

void EventInterval::recordEvent(const double current_time)
{
   d_time_at_last_event = current_time;
}

bool EventInterval::eventOccurredAtTime(const double current_time)
{
   bool occurred = false;

   if (isActive() && d_time_at_last_event == current_time) {

      occurred = true;
   }

   return occurred;
}

bool EventInterval::includeInitial(const double current_time)
{
   bool do_initial = false;

   if (isActive() && d_include_first_step &&
       !eventOccurredAtTime(current_time)) {

      do_initial = true;
      recordEvent(current_time);
   }

   d_previous_time = current_time;

   return do_initial;
}

bool EventInterval::includeFinal(const double current_time)
{
   bool do_final = false;

   if (isActive() && d_include_last_step &&
       !eventOccurredAtTime(current_time)) {

      do_final = true;
      recordEvent(current_time);
   }

   d_previous_time = current_time;

   return do_final;
}

bool EventInterval::hasIntervalPassed(const int current_step,
                                      const double current_time)
{
   bool has_it_passed = false;

   if (isActive() && !eventOccurredAtTime(current_time)) {

      if (d_include_first_step && d_previous_time < 0.0) {

         has_it_passed = true;

      } else if (d_interval_type == CYCLE) {

         if ((current_step % d_int_interval_value) == 0) {
            has_it_passed = true;
         }

      } else if (d_previous_time >= 0.0) {

         int n_prev = (int)(d_previous_time / d_dbl_interval_value);
         int n = (int)(current_time / d_dbl_interval_value);
         if (n > n_prev) {
            has_it_passed = true;
         }
      }
   }

   d_previous_time = current_time;

   if (has_it_passed) {
      recordEvent(current_time);
   }

   return has_it_passed;
}
