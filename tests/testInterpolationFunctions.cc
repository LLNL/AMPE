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
#include "FuncFort.h"

#include <iostream>
#include <math.h>

using namespace std;


int main(int argc, char *argv[])
{
   cout << "Test interpolation functions..." << endl;

   const double tol = 1.e-8;

   {
      std::string interp_func_type = "pbg";

      double val = INTERP_FUNC(0.5, interp_func_type.c_str());
      if (fabs(val - 0.5) > tol) {
         cerr << "Test failed!" << endl;
         return 1;
      }

      val = INTERP_FUNC(-0.5, interp_func_type.c_str());
      if (fabs(val - 0.) > tol) {
         cerr << "Test failed!" << endl;
         return 1;
      }
   }

   /*
    * Ratio pbg/lin
    */
   {
      std::string interp_func_type1 = "pbg";
      std::string interp_func_type2 = "lin";
      double phi = 0.05;
      double val1 = INTERP_RATIO_FUNC(phi, interp_func_type1.c_str(),
                                      interp_func_type2.c_str());
      double val2 = phi * phi * (10. - 15. * phi + 6 * phi * phi);
      if (fabs(val1 - val2) > tol) {
         cerr << "Test failed!" << endl;
         return 1;
      }
   }

   {
      std::string interp_func_type1 = "pbg";
      std::string interp_func_type2 = "lin";
      double phi = 0.05;
      double val1 = COMPL_INTERP_RATIO_FUNC(phi, interp_func_type1.c_str(),
                                            interp_func_type2.c_str());
      double a = 1. - INTERP_FUNC(phi, interp_func_type1.c_str());
      double b = 1. - INTERP_FUNC(phi, interp_func_type2.c_str());

      if (fabs(val1 - (a / b)) > tol) {
         cerr << "Test failed!" << endl;
         return 1;
      }
   }

   cout << "TEST successful!" << endl;

   return (0);
}
