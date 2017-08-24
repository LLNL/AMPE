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
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include "QuatParams.h"
#include "ConcFort.h"
#include "QuatFort.h"
#include "CompositionRHSStrategy.h"

#include <cassert>

using namespace std;

CompositionRHSStrategy::CompositionRHSStrategy(
   const string& avg_func_type)
   : d_avg_func_type( avg_func_type )
{
};

void CompositionRHSStrategy::addFluxFromGradTonPatch(
   hier::Patch& patch,
   const int temperature_id,
   const int flux_id)
{
   TBOX_ERROR( "CompositionRHSStrategy::addFluxFromGradTonPatch() not implemented..."<< endl );
}

void CompositionRHSStrategy::addFluxFromAntitrappingonPatch(
   hier::Patch& patch,
   const int phase_scratch_id,
   const int dphidt_id,
   const double alpha,
   const int flux_id)
{
   TBOX_ERROR( "CompositionRHSStrategy::addFluxFromAntitrappingonPatch() not implemented..."<< endl );
}

//-----------------------------------------------------------------------

void CompositionRHSStrategy::setZeroFluxAtBoundaryOnPatch(
   hier::Patch& patch,
   const int flux_id)
{
   boost::shared_ptr< geom::CartesianPatchGeometry > pg (
      BOOST_CAST< geom::CartesianPatchGeometry , hier::PatchGeometry>(patch.getPatchGeometry()) );

   // Get face boundary boxes.
   const std::vector< hier::BoundaryBox > bdry =
      pg->getCodimensionBoundaries(1);

   boost::shared_ptr< pdat::SideData<double> > flux (
      BOOST_CAST< pdat::SideData<double>, hier::PatchData>(patch.getPatchData( flux_id) ) );
   assert( flux );

   const hier::Box& gbox=flux->getGhostBox();
   const hier::Box& pbox=patch.getBox();


   const hier::Index plo(pbox.lower());
   const hier::Index pup(pbox.upper());

   const hier::Index glo(gbox.lower());
   for(int i=0;i<NDIM;i++)assert( glo(i)==plo(i) );
   
   for ( int i = 0; i < bdry.size(); i++ ) {
   
      const int locind=bdry[i].getLocationIndex();
      const int dir =locind>>1;
      const int side=locind&1;
      if( pg->getTouchesRegularBoundary(dir,side) )
      {
         //cout<<"EBSCompositionRHSStrategy::setZeroFluxAtBoundaryOnPatch"
         //          <<" for side "<<side<<" in direction "<<dir<<endl;

         hier::Index gup(gbox.upper());
         gup(dir)++;

         hier::Index lo(plo);
         hier::Index up(pup);
         
         // set direction perpendicular to boundary
         if( side==0 ){ // lower side
            up(dir)=plo(dir);
            lo(dir)=plo(dir);
         }else{         // upper side
            lo(dir)=pup(dir)+1;
            up(dir)=pup(dir)+1;
         }

#if (NDIM==3)
         FORT_SETTOZERO(lo(0), lo(1), lo(2), up(0), up(1), up(2),
                   glo(0), glo(1), glo(2), gup(0), gup(1), gup(2),
                   flux->getPointer(dir));
#endif
#if (NDIM==2)
         FORT_SETTOZERO(lo(0), lo(1), up(0), up(1),
                   glo(0), glo(1), gup(0), gup(1), 
                   flux->getPointer(dir));
#endif
      }
   } // loop over bdry boxes

}

