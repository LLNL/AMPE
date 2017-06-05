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
#ifndef included_FreeEnergyFunctions
#define included_FreeEnergyFunctions 

#include "Phases.h"

class FreeEnergyFunctions
{
public:
   FreeEnergyFunctions(){};

   ~FreeEnergyFunctions()
   {
   };
   
   virtual void energyVsPhiAndC(const double temperature, 
                                const double* const ceq,
                                const bool found_ceq,
                                const bool third_phase,
                                const int npts_phi=51,
                                const int npts_c=50)=0;
   virtual void printEnergyVsComposition(const double temperature, 
                                         const int npts=100 )=0;

   virtual bool computeCeqT(
      const double temperature,
      const PHASE_INDEX pi0, const PHASE_INDEX pi1,
      double* ceq,
      const bool verbose = false ){ return false; };

   virtual double computeFreeEnergy(
      const double temperature,
      const double* const conc,
      const PHASE_INDEX pi,
      const bool gp=false  )=0;
   virtual double computeDerivFreeEnergy(
      const double temperature,
      const double* const conc,
      const PHASE_INDEX pi )=0;
      
   virtual void preRunDiagnostics(std::ostream& os)=0;
};

#endif

