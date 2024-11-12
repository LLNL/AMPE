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
#include "TiltingFolchPlapp2005.h"
#include "FuncFort.h"

#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;


int main(int argc, char *argv[])
{
   std::cout << "Test 3 phases interpolation functions..." << endl;
   std::cout << setprecision(12);

   const double tol = 2.e-8;

   cout << "Test FolchPlapp2005..." << endl;
   {
      double val = TiltingFolchPlapp2005::g0(1., 0., 0.);
      if (fabs(val - 1.) > tol) {
         cerr << "Test failed!" << endl;
         return 1;
      }

      double phi = 0.2;
      double val0 = TiltingFolchPlapp2005::g0(1. - phi, phi, 0.);
      double val1 = 1. - TiltingFolchPlapp2005::g0(phi, 1. - phi, 0.);
      if (fabs(val1 - val0) > tol) {
         cerr << "Test failed: g0(1-phi,phi,0)!=1-g0(phi,1-phi,0)!" << endl;
         return 1;
      }
   }

   /*
    * Derivatives
    */
   {
      double p0 = 0.45;
      double p1 = 0.33;
      double p2 = 0.22;
      // compare with finite difference
      double eps = 1.e-8;
      double val0 = TiltingFolchPlapp2005::g0(p0, p1, p2);
      double val1 = TiltingFolchPlapp2005::g0(p0 + 2. * eps / 3., p1 - eps / 3.,
                                              p2 - eps / 3.);
      double deriv0 = (val1 - val0) / eps;
      double deriv1 = TiltingFolchPlapp2005::dg0dp0(p0, p1, p2);
      cout << "FolchPlapp2005 FD deriv=" << deriv0
           << ", Analytical deriv=" << deriv1 << std::endl;
      if (fabs(deriv1 - deriv0) > tol) {
         return 1;
      }

      val0 = TiltingFolchPlapp2005::g0(p0, p1, p2);
      val1 = TiltingFolchPlapp2005::g0(p0 - eps / 3., p1 + 2. * eps / 3.,
                                       p2 - eps / 3.);
      deriv0 = (val1 - val0) / eps;
      deriv1 = TiltingFolchPlapp2005::dg0dp1(p0, p1, p2);
      cout << "FolchPlapp2005 FD deriv=" << deriv0
           << ", Analytical deriv=" << deriv1 << std::endl;
      if (fabs(deriv1 - deriv0) > tol) {
         return 1;
      }
   }

   // verify derivatives are 0 at (1,0,0)
   {
      double p0 = 1.;
      double p1 = 0.;
      double p2 = 0.;

      double val01 = TiltingFolchPlapp2005::dg0dp1(p0, p1, p2);
      cout << "val01 = " << val01 << std::endl;
      if (fabs(val01) > tol) {
         return 1;
      }

      double val00 = TiltingFolchPlapp2005::dg0dp0(p0, p1, p2);
      cout << "val00 = " << val00 << std::endl;
      if (fabs(val00) > tol) {
         return 1;
      }
   }

   {
      double p0 = 0.9;
      double p1 = 0.;
      double p2 = 0.1;

      double val01 = TiltingFolchPlapp2005::dg0dp1(p0, p1, p2);
      cout << "val01 = " << val01 << std::endl;
      if (fabs(val01) > tol) {
         return 1;
      }
   }

   {
      double p0 = 0.;
      double p1 = 1.;
      double p2 = 0.;

      double val01 = TiltingFolchPlapp2005::dg0dp1(p0, p1, p2);
      cout << "val01 = " << val01 << std::endl;
      if (fabs(val01) > tol) {
         return 1;
      }
   }


   std::cout << "Test Triple Well..." << endl;
   {
      double p0 = 0.1;
      double p1 = 0.7;
      double p2 = 0.2;
      double eps = 1.e-8;

      double val0 = TRIPLE_WELL_FUNC(p0, p1, p2);
      double val1 =
          TRIPLE_WELL_FUNC(p0 + 2. * eps / 3., p1 - eps / 3., p2 - eps / 3.);
      double deriv1 = DERIV_TRIPLE_WELL_FUNC(p0, p1, p2);
      double deriv0 = (val1 - val0) / eps;

      std::cout << "Triple Well FD deriv=" << deriv0
                << ", Analytical deriv=" << deriv1 << std::endl;
      if (fabs(deriv1 - deriv0) > 3.e-8) {
         return 1;
      }
   }

   cout << "TEST successful!" << endl;

   return (0);
}
