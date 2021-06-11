// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
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
