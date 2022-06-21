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
//
#include "TiltingFolchPlapp2005.h"

#include <iostream>
#include <math.h>

using namespace std;


int main(int argc, char *argv[])
{
   cout << "Test 3 phases interpolation functions..." << endl;

   const double tol = 1.e-8;

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
      double p0 = 0.22;
      double p1 = 0.33;
      double p2 = 0.45;
      // compare with finite difference
      double eps = 1.e-8;
      double val0 = TiltingFolchPlapp2005::g0(p0, p1, p2);
      double val1 = TiltingFolchPlapp2005::g0(p0 + 2. * eps / 3., p1 - eps / 3.,
                                              p2 - eps / 3.);
      double deriv0 = (val1 - val0) / eps;
      double deriv1 = TiltingFolchPlapp2005::dg0dp0(p0, p1, p2);
      cerr << "FD deriv=" << deriv0 << ", Analytical deriv=" << deriv1
           << std::endl;
      if (fabs(deriv1 - deriv0) > tol) {
         return 1;
      }

      val0 = TiltingFolchPlapp2005::g0(p0, p1, p2);
      val1 = TiltingFolchPlapp2005::g0(p0 - eps / 3., p1 + 2. * eps / 3.,
                                       p2 - eps / 3.);
      deriv0 = (val1 - val0) / eps;
      deriv1 = TiltingFolchPlapp2005::dg0dp1(p0, p1, p2);
      cerr << "FD deriv=" << deriv0 << ", Analytical deriv=" << deriv1
           << std::endl;
      if (fabs(deriv1 - deriv0) > tol) {
         return 1;
      }
   }

   cout << "TEST successful!" << endl;

   return (0);
}
