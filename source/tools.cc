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
#include "tools.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/InputManager.h"

using namespace SAMRAI;

//-----------------------------------------------------------------------

void listLocalToGlobal(std::map<int, double>& lcl_map,
                       std::map<int, double>& gbl_map)
{
   int size_in = static_cast<int>(lcl_map.size());
   int* lcl_n = new int[size_in];
   double* lcl_v = new double[size_in];

   int ii = 0;
   for (std::map<int, double>::const_iterator it = lcl_map.begin();
        it != lcl_map.end(); it++) {
      lcl_n[ii] = it->first;
      lcl_v[ii] = it->second;
      ii++;
   }

   int size_out = sumReduction((int)size_in);
   int* gbl_n = new int[size_out];
   double* gbl_v = new double[size_out];

   allGatherv(lcl_n, size_in, gbl_n, size_out);
   allGatherv(lcl_v, size_in, gbl_v, size_out);

   for (ii = 0; ii < size_out; ii++) {
      gbl_map[gbl_n[ii]] += gbl_v[ii];
   }

   delete[] lcl_n;
   delete[] lcl_v;
   delete[] gbl_n;
   delete[] gbl_v;
}

void sumReduction(double* x, const int count)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   mpi.AllReduce(x, count, MPI_SUM);
}

double sumReduction(const double x)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   double val = x;
   mpi.AllReduce(&val, 1, MPI_SUM);
   return val;
}

int sumReduction(const int x)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   int val = x;
   mpi.AllReduce(&val, 1, MPI_SUM);
   return val;
}

// common setup funtion for all-to-all functions                         *
void allGatherSetup(int size_in, const int size_out, int*& rcounts, int*& disps)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   int np = mpi.getSize();

   rcounts = new int[np];
   disps = new int[np];

   /* figure out where where each processor's input will be placed */
   MPI_Allgather(&size_in, 1, MPI_INT, rcounts, 1, MPI_INT,
                 mpi.getCommunicator());

   disps[0] = 0;
   for (int p = 1; p < np; ++p) {
      disps[p] = disps[p - 1] + rcounts[p - 1];
   }

   /* verify that the x_out array is the appropriate size! */
   int c = 0;
   for (int x = 0; x < np; ++x) {
      c += rcounts[x];
   }
   if (c != size_out) {
      TBOX_ERROR("allGatherSetup error..."
                 << "\n   size_out =" << size_out << "appears to be incorrect; "
                 << "should be: " << c << std::endl);
   }
}

void allGatherv(const int* x_in, int size_in, int* x_out, int size_out)
{
   int* rcounts = (int*)NULL;
   int* disps = (int*)NULL;
   allGatherSetup(size_in, size_out, rcounts, disps);

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   MPI_Allgatherv((void*)x_in, size_in, MPI_INT, x_out, rcounts, disps, MPI_INT,
                  mpi.getCommunicator());

   if (rcounts) {
      delete[] rcounts;
   }
   if (disps) {
      delete[] disps;
   }
}

void allGatherv(const double* x_in, int size_in, double* x_out, int size_out)
{
   int* rcounts = (int*)NULL;
   int* disps = (int*)NULL;
   allGatherSetup(size_in, size_out, rcounts, disps);

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   MPI_Allgatherv((void*)x_in, size_in, MPI_DOUBLE, x_out, rcounts, disps,
                  MPI_DOUBLE, mpi.getCommunicator());

   if (rcounts) {
      delete[] rcounts;
   }
   if (disps) {
      delete[] disps;
   }
}

void printDeprecated(const std::string& s_old, const std::string& s_new)
{
   tbox::pout << "Input " << s_old << " is deprecated.  Use " << s_new << "."
              << std::endl;
}
