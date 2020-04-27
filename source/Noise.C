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


void Noise::setField(std::shared_ptr<hier::Patch> patch, const int data_id,
                     const int phi_id)
{
   assert(data_id >= 0);

   std::shared_ptr<pdat::CellData<double> > cd(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(data_id)));
   assert(cd);

   std::shared_ptr<pdat::CellData<double> > phi(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(phi_id)));
   assert(phi);
   assert(phi->getGhostCellWidth()[0] >= cd->getGhostCellWidth()[0]);

   double* ptr_data = cd->getPointer();
   double* ptr_phi = phi->getPointer();

   const hier::Box& gbox = cd->getGhostBox();
   int imin[3] = {gbox.lower(0), gbox.lower(1), 0};
   int imax[3] = {gbox.upper(0), gbox.upper(1), 0};
#if (NDIM == 3)
   imin[2] = gbox.lower(2);
   imax[2] = gbox.upper(2);
#endif

   const hier::Box& phi_gbox = phi->getGhostBox();
   int min_phi[3] = {phi_gbox.lower(0), phi_gbox.lower(1), 0};
   int inc_j_phi = phi_gbox.numberCells(0);
#if (NDIM == 3)
   min_phi[2] = phi_gbox.lower(2);
   int inc_k_phi = inc_j_phi * phi_gbox.numberCells(1);
#else
   int inc_k_phi = 0;
#endif

   int idx = 0;
   int idx_phi = (imin[0] - min_phi[0]) + (imin[1] - min_phi[1]) * inc_j_phi +
                 (imin[2] - min_phi[2]) * inc_k_phi;

   for (int kk = imin[2]; kk <= imax[2]; kk++) {
      for (int jj = imin[1]; jj <= imax[1]; jj++) {
         for (int ii = imin[0]; ii <= imax[0]; ii++) {

            ptr_data[idx] +=
                gen() * 4. * ptr_phi[idx_phi] * (1. - ptr_phi[idx_phi]);
            idx++;
            idx_phi++;
         }  // ii
         idx_phi +=
             2 * (phi->getGhostCellWidth()[0] - cd->getGhostCellWidth()[0]);
      }  // jj
      idx_phi += 2 * inc_j_phi *
                 (phi->getGhostCellWidth()[1] - cd->getGhostCellWidth()[1]);
   }  // kk
}
