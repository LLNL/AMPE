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
#include "math_utilities.h"
#include <math.h>

using namespace std;

//-----------------------------------------------------------------------

double Determinant3( double** const m )
{
   double d =
         m[0][0] * m[1][1] * m[2][2] -
         m[0][0] * m[1][2] * m[2][1] -
         m[0][1] * m[1][0] * m[2][2] +
         m[0][1] * m[1][2] * m[2][0] +
         m[0][2] * m[1][0] * m[2][1] -
         m[0][2] * m[1][1] * m[2][0];

   return d;
}

//=======================================================================

double Determinant4( double** const m )
{
   double d =
         m[0][0] * (m[1][1] * m[2][2] * m[3][3] -
                    m[1][1] * m[2][3] * m[3][2] -
                    m[1][2] * m[2][1] * m[3][3] +
                    m[1][2] * m[2][3] * m[3][1] +
                    m[1][3] * m[2][1] * m[3][2] -
                    m[1][3] * m[2][2] * m[3][1])

      - m[0][1] * (m[1][0] * m[2][2] * m[3][3] -
                   m[1][0] * m[2][3] * m[3][2] -
                   m[1][2] * m[2][0] * m[3][3] +
                   m[1][2] * m[2][3] * m[3][0] +
                   m[1][3] * m[2][0] * m[3][2] -
                   m[1][3] * m[2][2] * m[3][0] )

      + m[0][2] * (m[1][0] * m[2][1] * m[3][3] -
                   m[1][0] * m[2][3] * m[3][1] -
                   m[1][1] * m[2][0] * m[3][3] +
                   m[1][1] * m[2][3] * m[3][0] +
                   m[1][3] * m[2][0] * m[3][1] -
                   m[1][3] * m[2][1] * m[3][0] )

      - m[0][3] * (m[1][0] * m[2][1] * m[3][2] -
                   m[1][0] * m[2][2] * m[3][1] -
                   m[1][1] * m[2][0] * m[3][2] +
                   m[1][1] * m[2][2] * m[3][0] +
                   m[1][2] * m[2][0] * m[3][1] -
                   m[1][2] * m[2][1] * m[3][0] );

   return d;
}

//=======================================================================

double DeterminantN(double** mat, const short n)
{
   double* submat[5];  
   if (n == 2) 
   {
      return( (mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
   }
   else
   {  
      double d=0.;
      for(short c = 0; c < n; c++)
      {  
         //loop over rows
         short subi = 0;
         for(short i = 0; i < n; i++)
         {
            if (i==c)continue;

            submat[subi]=&mat[i][1];
            subi++;
         }
         d += (pow(-1 ,c) * mat[c][0] * DeterminantN(submat,n-1));
      }
      return d;
   }
}
