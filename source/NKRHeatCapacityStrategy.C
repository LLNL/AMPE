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
#include "NKRHeatCapacityStrategy.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"

#include <set>
using namespace std;

void NKRHeatCapacityStrategy::setCurrentValue(
   boost::shared_ptr<hier::PatchHierarchy > patch_hierarchy)
{
   assert( d_cp_id>=0 );
   assert( d_temperature_id>=0 );
   assert( d_composition_id>=0 );
   assert( !d_cp.empty() );
   
   //build a set with all powers used in T expansion, for all species
   set<short> powers;
   for(vector< map<short,double> >::iterator it=d_cp.begin(); it!=d_cp.end(); ++it)
   {
      for(map<short,double>::iterator itm= it->begin(); itm!=it->end(); ++itm)
      {
         powers.insert(itm->first);
      }
   }
   
   //store these powers in an array we can give as an argument to fortran function
   int* cp_powers=new int[powers.size()];
   short itp=0;
   for(set<short>::iterator its=powers.begin();its!=powers.end();its++)
   {
      cp_powers[itp]=*its;
      itp++;
   }
   
   short npowers=powers.size();
   double* cp_coeffs = new double[d_cp.size()*npowers];
   memset(cp_coeffs,0,d_cp.size()*npowers*sizeof(double));

   //loop over species
   //tbox::pout<<"Setting up Heat capacity using "<<d_cp.size()<<"species"<<endl;
   short isp=0;
   for(vector< map<short,double> >::iterator it=d_cp.begin(); it!=d_cp.end(); ++it)
   {
      //loop over expansion powers
      short ip=0;
      for(set<short>::iterator its=powers.begin();its!=powers.end();its++)
      {
         for(map<short,double>::iterator itm= it->begin(); itm!=it->end(); ++itm)
         {
            if( *its==itm->first ){
                cp_coeffs[isp*npowers+ip]=itm->second;
            }
         }
         ip++;
      }
      isp++;
   }
   
   int maxln = patch_hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++ ) {

      boost::shared_ptr<hier::PatchLevel > level = patch_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p!=level->end(); p++) {
         boost::shared_ptr<hier::Patch > patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast =  pbox.upper();

         boost::shared_ptr< pdat::CellData<double> > conc (
            patch->getPatchData( d_composition_id ), boost::detail::dynamic_cast_tag());
         boost::shared_ptr< pdat::CellData<double> > temp (
            patch->getPatchData( d_temperature_id ), boost::detail::dynamic_cast_tag());
         boost::shared_ptr< pdat::CellData<double> > cp (
            patch->getPatchData( d_cp_id ), boost::detail::dynamic_cast_tag());
            
         FORT_HEAT_CAPACITY_NKR(
            ifirst(0), ilast(0),
            ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            conc->getPointer(), conc->getGhostCellWidth()[0],
            temp->getPointer(), temp->getGhostCellWidth()[0],
            cp->getPointer(),  cp->getGhostCellWidth()[0],
            cp_powers, npowers,
            cp_coeffs, d_cp.size()*npowers);
      }
   }
   
   delete[] cp_coeffs;
   delete[] cp_powers;
}
