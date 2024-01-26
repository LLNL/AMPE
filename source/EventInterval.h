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
#ifndef included_EventInterval
#define included_EventInterval

#include "SAMRAI/tbox/Database.h"

#include <string>

using namespace SAMRAI;

class EventInterval
{
 public:
   enum INTERVAL_TYPE { CYCLE = 0, TIME = 1 };

   EventInterval(std::shared_ptr<tbox::Database> input_db,
                 const std::string name, const double default_value = 0.0,
                 const std::string default_type = "step",
                 const bool include_first = false,
                 const bool include_last = true);

   virtual ~EventInterval() { ; }

   bool isActive();

   bool eventOccurredAtTime(const double current_time);

   bool includeInitial(const double current_time);
   bool includeFinal(const double current_time);

   bool hasIntervalPassed(const int current_step, const double current_time);

 private:
   void recordEvent(const double current_time);

   int d_int_interval_value;
   double d_dbl_interval_value;
   INTERVAL_TYPE d_interval_type;
   double d_previous_time;
   double d_time_at_last_event;
   bool d_include_first_step;
   bool d_include_last_step;
   const std::string d_name;
};

#endif
