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
#ifndef FuncFort_H
#define FuncFort_H

// Link between C/C++ and Fortran files
//       name in             name in
//      C/C++ code            Fortran code
//      ----------            ------------
#define FORT_WELL_FUNC well_func_
#define FORT_DERIV_WELL_FUNC deriv_well_func_
#define FORT_SECOND_DERIV_WELL_FUNC second_deriv_well_func_
#define FORT_INTERP_FUNC interp_func_
#define FORT_DERIV_INTERP_FUNC deriv_interp_func_
#define FORT_SECOND_DERIV_INTERP_FUNC second_deriv_interp_func_
#define FORT_AVERAGE_FUNC average_func_
#define FORT_DERIV_AVERAGE_FUNC deriv_average_func_

// Function argument list interfaces
extern "C" {
   double FORT_WELL_FUNC( const double&, const char* );
   double FORT_DERIV_WELL_FUNC( const double&, const char* );
   double FORT_SECOND_DERIV_WELL_FUNC( const double&, const char* );
   double FORT_INTERP_FUNC( const double&, const char* );
   double FORT_DERIV_INTERP_FUNC( const double&, const char* );
   double FORT_SECOND_DERIV_INTERP_FUNC( const double&, const char* );
   double FORT_AVERAGE_FUNC( const double&, const double&, const char* );
   double FORT_DERIV_AVERAGE_FUNC( const double&, const double&, const char* );
}

#endif

