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
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* Long period (> 2 × 1018) random number generator of L'Ecuyer with Bays-Durham shuffle
   and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
   the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
   idum between successive deviates in a sequence. RNMX should approximate the largest floating
   value that is less than 1. */

double ran2(long *idum)
{
   int j;
   long k;
   static long idum2=123456789;
   static long iy=0;
   static long iv[NTAB];
   double temp;
   if (*idum <= 0) {  // Initialize.
      if (-(*idum) < 1) *idum=1; // Be sure to prevent idum = 0.
      else *idum = -(*idum);
      idum2=(*idum);
      for (j=NTAB+7;j>=0;j--) {  // Load the shuffle table (after 8 warm-ups).
         k=(*idum)/IQ1;
         *idum=IA1*(*idum-k*IQ1)-k*IR1;
         if (*idum < 0) *idum += IM1;
         if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ1;  // Start here when not initializing.
   *idum=IA1*(*idum-k*IQ1)-k*IR1; // Compute idum=(IA1*idum) % IM1 without
                                  //overflows by Schrage\u2019s if (*idum < 0) *idum += IM1; method.
   k=idum2/IQ2;
   idum2=IA2*(idum2-k*IQ2)-k*IR2;  // Compute idum2=(IA2*idum) % IM2 likewise.
   if (idum2 < 0) idum2 += IM2;
   j=iy/NDIV;  // Will be in the range 0..NTAB-1.
   iy=iv[j]-idum2;  // Here idum is shuffled, idum and idum2 are
   iv[j] = *idum;  // combined to generate output.
   if (iy < 1) iy += IMM1;
   if ((temp=AM*iy) > RNMX) return RNMX;  // Because users don't expect endpoint values.
   else return temp;
}

