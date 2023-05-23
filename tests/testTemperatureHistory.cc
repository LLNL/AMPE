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
#include "TemperatureHistory.h"

#include <iostream>
#include <cmath>

int main(int argc, char* argv[])
{
   int ret = 0;
   double tol = 1.e-8;

   const std::string filename("temperature.csv");
   TemperatureHistory temp_history(std::cout);
   temp_history.readCSV(filename);

   double temperature = 0.;
   double gradT = 0.;
   int th = 0;

   // test 1
   th = temp_history.getTandGradT(0.0014, -1.4e-2, temperature, gradT);
   std::cout << "Test 1: value for gradT: " << gradT << std::endl;
   if (th > 0) {
      std::cerr << "t,x outside domain..." << std::endl;
      ret = 1;
   } else if (std::fabs(gradT - 150.) > tol) {
      std::cerr << "Incorrect value for gradT" << std::endl;
      ret = 1;
   }

   // test 2
   th = temp_history.getTandGradT(0.0034, -3.4e-2, temperature, gradT);
   std::cout << "Test 2: value for gradT: " << gradT << std::endl;
   if (th > 0) {
      std::cerr << "t,x outside domain..." << std::endl;
      ret = 1;
   } else if (std::fabs(gradT - 150.) > tol) {
      std::cerr << "Incorrect value for gradT" << std::endl;
      ret = 1;
   }

   // test 3
   th = temp_history.getTandGradT(0.0034, -2.1e-2, temperature, gradT);
   std::cout << "Test 3: value for gradT: " << gradT << std::endl;
   if (th > 0) {
      std::cerr << "t,x outside domain..." << std::endl;
      ret = 1;
   } else if (std::fabs(gradT - 150.) > tol) {
      std::cerr << "Incorrect value for gradT" << std::endl;
      ret = 1;
   }

   // test 4: x outside domain
   th = temp_history.getTandGradT(0.0034, -9.4e-2, temperature, gradT);
   if (th == 0) {
      std::cerr << "Should return 1: x outside domain..." << std::endl;
      ret = 1;
   }

   // test 5: x outside domain
   th = temp_history.getTandGradT(0.0034, 1.4e-2, temperature, gradT);
   if (th == 0) {
      std::cerr << "Should return 1: x outside domain..." << std::endl;
      ret = 1;
   }

   if (ret == 0)
      std::cout << "\nPASSED" << std::endl;
   else
      std::cout << "\nFAILED" << std::endl;

   return ret;
}
