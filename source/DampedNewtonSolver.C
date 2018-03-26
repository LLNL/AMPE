// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
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
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
#include "DampedNewtonSolver.h"

#include <iostream>
#include <cmath>
#include <cassert>

#define DEBUG_CONVERGENCE
#ifdef DEBUG_CONVERGENCE
#include <vector>
#endif

#include <iomanip>

using namespace std;

//=======================================================================

DampedNewtonSolver::DampedNewtonSolver() :
      NewtonSolver(),
      d_alpha( 1. )
{};

   
//=======================================================================

void DampedNewtonSolver::UpdateSolution(
   double* const c,
   const double* const fvec,
   double** const fjac )
{
   int nn=size();

   static double* mwork[3];
   static double mtmp[9];
   for ( int ii = 0; ii < nn; ii++ ) {
      mwork[ii] = &mtmp[ii*nn];
   }

   const double D = Determinant( fjac );
   const double D_inv = 1.0 / D;

   //cout << "D = " << D << endl;

   static double del_c[3];

   // use Cramer's rule to solve linear system
   for ( int jj = 0; jj < nn; jj++ ) {

      CopyMatrix( mwork, fjac );
      for ( int ii = 0; ii < nn; ii++ ) {
         mwork[ii][jj] = fvec[ii];
      }

      del_c[jj] = D_inv * Determinant( mwork );

      //cout << "del_c[" << jj << "] = " << del_c[jj] << endl;

   }

#if 0
   double w = 1.0e20;

   for ( int ii = 0; ii < nn; ii++ ) {

      double ctilde = c[ii] - del_c[ii];

      double w_c = 1.0;

      if ( ctilde < 0. ) {
         w_c = ( del_c[ii] + ctilde ) / del_c[ii];
         //cout << "ctilde[" << ii << "]  = " << ctilde << endl;
      }
      else if ( ctilde > 1 ) {
         w_c = ( del_c[ii] + ctilde - 1.0 ) / del_c[ii];
         //cout << "ctilde[" << ii << "]  = " << ctilde << endl;
      }

      if ( w_c < w ) {
         w = w_c;
      }

   }
#else
   //double w = 1.0;
   double w = d_alpha;
#endif

   // make sure c remains between tol and 1-tol
   double tol=1.e-16;
   for ( int ii = 0; ii < nn; ii++ ) {

      if( c[ii] <tol )
      {
         c[ii]=tol;
         if( del_c[ii]<0.) del_c[ii]=0.;
      }
      if( c[ii] >1.-tol )
      {
         c[ii]=1.-tol;
         if( del_c[ii]>0.) del_c[ii]=0.;
      }

   }
   for ( int ii = 0; ii < nn; ii++ ) {

      c[ii] = c[ii] - w * del_c[ii];

   }
   
   bool flag;
   do
   {
      flag=false;
      for ( int ii = 0; ii < nn; ii++ ) {
         if(c[ii]<0. || c[ii]>1. )
         {
            w*=0.5;
            //cout<<"c="<<c[ii]<<", rescale w..."<<w<<endl;
            
            for ( int jj = 0; jj < nn; jj++ )
            {
               c[jj] += 2.*w * del_c[jj];
               c[jj] -=    w * del_c[jj];
            }
            flag=true;
            break;
         }
      }
   
   }while( flag );

   //cout << endl;
}
