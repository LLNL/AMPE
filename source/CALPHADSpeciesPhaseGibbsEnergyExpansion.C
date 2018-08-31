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
#include "CALPHADSpeciesPhaseGibbsEnergyExpansion.h"

#include <math.h>

CALPHADSpeciesPhaseGibbsEnergyExpansion::
   CALPHADSpeciesPhaseGibbsEnergyExpansion(
      const double a,
      const double b,
      const double c,
      const double d2,
      const double d3,
      const double d4,
      const double d7,
      const double dm1,
      const double dm9):
         d_a(a),
         d_b(b),
         d_c(c),
         d_d2(d2),
         d_d3(d3),
         d_d4(d4),
         d_d7(d7),
         d_dm1(dm1),
         d_dm9(dm9)
{
}

double CALPHADSpeciesPhaseGibbsEnergyExpansion::value(
   const double temperature)const
{
   const double t2=temperature*temperature;
   const double t4=t2*t2;

   return d_a
         +d_b*temperature
         +d_c*temperature*log(temperature)
         +d_d2*t2
         +d_d3*t2*temperature
         +d_d4*t4
         +d_d7*t4*t2*temperature
         +d_dm1/temperature
         +d_dm9/(t4*t4*temperature);
}

