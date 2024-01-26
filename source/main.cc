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

#include "AMPE.h"

#include <string>


int main(int argc, char *argv[])
{
   // Initialize MPI
   MPI_Init(&argc, &argv);

   {
      // create main object
      AMPE ampe(MPI_COMM_WORLD);

      /*
       * Process command line arguments.
       * For non-restarted case, command line is:
       *
       *    executable <input file name>
       *
       * For restarted run, command line is:
       *
       *    executable <input file name> <restart directory> \
       *               <restart number>
       */
      if ((argc != 2) && (argc != 4)) {
         std::cerr << "USAGE:  " << argv[0] << " <input filename> "
                   << "<restart dir> <restore number> [options]\n"
                   << "  options:\n"
                   << "  none at this time" << std::endl;
         return (-1);
      }

      std::string input_filename = argv[1];
      int restore_num = -1;
      std::string restart_read_dirname = "";
      if (argc == 4) {
         restart_read_dirname = argv[2];
         restore_num = atoi(argv[3]);
      }

      ampe.initialize(input_filename, restart_read_dirname, restore_num);

      ampe.run();
   }

   MPI_Finalize();

   return (0);
}
