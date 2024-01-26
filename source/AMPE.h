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
#ifndef included_AMPE
#define included_AMPE

#include "SAMRAI/tbox/TimerManager.h"

#include <cstring>

class PFModel;

class AMPE
{
 public:
   AMPE(MPI_Comm comm);

   ~AMPE();

   std::string gitCommitID() const
   {
#define xstr2(x) #x
#define xstr(x) xstr2(x)

      return xstr(GITVERSION);
   }

   void initialize(const std::string input_filename,
                   const std::string restart_read_dirname = "",
                   const int restore_num = -1);

   void run();

 private:
   PFModel* d_pfm;

   SAMRAI::tbox::TimerManager* d_time_man;
};

#endif
