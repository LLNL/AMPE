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
#include "FuncFort.h"

#include <iostream>
#include <math.h>

using namespace std;


int main( int argc, char *argv[] )
{
   cout<<"Test interpolation functions..."<<endl;

   const double tol = 1.e-8;

   {
      std::string interp_func_type = "pbg";

      double val = FORT_INTERP_FUNC( 0.5,  interp_func_type.c_str() );
      if( fabs( val-0.5 )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }

      val = FORT_INTERP_FUNC( -0.5,  interp_func_type.c_str() );
      if( fabs( val-0. )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }
   }

   /*
    * Ratio pbg/lin
    */
   {
      std::string interp_func_type1 = "pbg";
      std::string interp_func_type2 = "lin";
      double phi=0.05;
      double val1 = FORT_INTERP_RATIO_FUNC( phi, interp_func_type1.c_str(),
                                                 interp_func_type2.c_str() );
      double val2 = phi*phi*(10.-15.*phi+6*phi*phi);
      if( fabs( val1-val2 )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }

   }

   {
      std::string interp_func_type1 = "pbg";
      std::string interp_func_type2 = "lin";
      double phi=0.05;
      double val1 = FORT_COMPL_INTERP_RATIO_FUNC( phi,
                                                  interp_func_type1.c_str(),
                                                  interp_func_type2.c_str() );
      double a = 1.-FORT_INTERP_FUNC( phi, interp_func_type1.c_str() );
      double b = 1.-FORT_INTERP_FUNC( phi, interp_func_type2.c_str() );

      if( fabs( val1-(a/b) )>tol ){
         cerr<<"Test failed!"<<endl;
         return 1;
      }

   }

   cout<<"TEST successful!"<<endl;

   return(0);
}
