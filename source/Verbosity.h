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
#ifndef included_Verbosity
#define included_Verbosity

class Verbosity
{
 public:
   Verbosity(){};
   virtual ~Verbosity(){};
   virtual bool isSilent(void) = 0;
   virtual bool notSilent(void) = 0;
   virtual void setBasicLevel(int l) = 0;
   virtual void setAmrLevel(int l) = 0;
   virtual int basicLevel(void) = 0;
   virtual int amrLevel(void) = 0;
};

class QuatVerbosity : public Verbosity
{
 public:
   QuatVerbosity(void) : i_basic_level(1), i_amr_level(0){};
   ~QuatVerbosity(){};
   bool isSilent(void) { return basicLevel() == 0; }
   bool notSilent(void) { return !isSilent(); }
   void setBasicLevel(int l) { i_basic_level = l; }
   void setAmrLevel(int l)
   {
      i_amr_level = l;
      if (isSilent()) i_amr_level = 0;
   }
   int basicLevel(void) { return i_basic_level; }
   int amrLevel(void) { return i_amr_level; }

 private:
   int i_basic_level;
   int i_amr_level;
};

#endif
