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
#include "Noise.h"

#include "SAMRAI/pdat/CellData.h"

Noise* Noise::s_pinstance=0;
int Noise::s_data_id=-1;
std::unique_ptr<boost::variate_generator<rng_type,distribution_type> >
   Noise::s_gen=0;

Noise::Noise()
{
   rng_type          rng;
   distribution_type nd(0.0, 1.0);

   s_gen = std::unique_ptr<gen_type>(new gen_type(rng,nd));
}

void Noise::setup(const int data_id)
{
   s_data_id=data_id;
}

void Noise::update(
   boost::shared_ptr<hier::Patch > patch)
{
   boost::shared_ptr< pdat::CellData<double> > cd(
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>(
         patch->getPatchData( s_data_id ) ) );
   assert( cd );

   double* ptr_data = cd->getPointer();

   const hier::Box& gbox = cd->getGhostBox();

   int imin[3] = {gbox.lower(0),gbox.lower(1),0};
   int imax[3] = {gbox.upper(0),gbox.upper(1),0};
#if (NDIM == 3)
   imin[2] = gbox.lower(2);
   imax[2] = gbox.upper(2);
#endif

   int idx=0;

   for ( int kk = imin[2]; kk <= imax[2]; kk++ ) {
      for ( int jj = imin[1]; jj <= imax[1]; jj++ ) {
         for ( int ii = imin[0]; ii <= imax[0]; ii++ ) {

            ptr_data[idx]=(*s_gen)();

            idx++;

         } // ii

      } // jj

   } // kk

}

