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
#include <iostream>

#include "NormalNoise.h"
#include "UniformNoise.h"


int main(int argc, char* argv[])
{
   std::cout << "Test noise generation." << std::endl;

   // 1st test
   {
      std::cout << "Normal distributed noise..." << std::endl;
      NormalNoise& noise(*(NormalNoise::instance()));

      int ntotal = 1000;
      int count = 0;
      double avg = 0.;
      for (int i = 0; i < ntotal; i++) {
         double val = noise.gen();
         // std::cout<<val<<endl;
         if (fabs(val) <= 1.) count++;
         avg += val;
      }
      avg /= (double)ntotal;

      std::cout << count << " out of " << ntotal
                << " values are between -1 and 1." << std::endl;
      std::cout << "Average value is " << avg << std::endl;

      // 68% of values should be between -1 and 1
      if ((double)count > 0.71 * (double)ntotal ||
          (double)count < 0.65 * (double)ntotal) {
         std::cerr << "TEST NormalNoise failed!!!" << std::endl;
         return 1;
      }

      if (fabs(avg) > 0.04) {
         std::cerr << "TEST average NormalNoise failed!!!" << std::endl;
         return 1;
      }
   }

   // 2nd test
   {
      std::cout << "Uniform distributed noise..." << std::endl;
      UniformNoise& noise(*(UniformNoise::instance(42u)));

      int ntotal = 1000;
      int count = 0;
      double avg = 0.;
      for (int i = 0; i < ntotal; i++) {
         double val = noise.gen();
         if (fabs(val) <= 0.5) count++;
         avg += val;
      }
      avg /= ntotal;

      std::cout << count << " out of " << ntotal
                << " values are between -0.5 and 0.5" << std::endl;
      std::cout << "Average value is " << avg << std::endl;

      if (count < ntotal) {
         std::cerr << "TEST UniformNoise failed!!!" << std::endl;
         return 1;
      }

      if (fabs(avg) > 0.02) {
         std::cerr << "TEST average UniformNoise failed!!!" << std::endl;
         return 1;
      }
   }

   std::cout << "TEST successful!" << std::endl;

   return 0;
}
