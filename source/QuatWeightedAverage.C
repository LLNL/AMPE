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
#include "QuatWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cassert>
using namespace std;

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
#if (NDIM == 2)
   void quatcoarsen_(
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int*, const double*, const double*,
      const double*, double*, const int&,
      const int&,
      const int*, const int*,
      const int&, const int&, const int&, const int&
      );
#endif
#if (NDIM == 3)
   void quatcoarsen_(
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int*, const double*, const double*,
      const double*, double*, const int&,
      const int&,
      const int*, const int*, const int*,
      const int&, const int&, const int&,
      const int&, const int&, const int&
      );
#endif
}

using namespace SAMRAI;

QuatWeightedAverage::QuatWeightedAverage(
   const bool symmetry_aware,
   const int quat_symm_rotation_id )
   : hier::CoarsenOperator("BASE_QUAT_COARSEN")
{
   d_name_id = "QUAT_COARSEN";
   d_symmetry_aware = symmetry_aware;
   d_quat_symm_rotation_id = quat_symm_rotation_id;
}

QuatWeightedAverage::~QuatWeightedAverage()
{
}

bool QuatWeightedAverage::findCoarsenOperator(
   const boost::shared_ptr< hier::Variable >& var,
   const string &op_name) const
{
   const boost::shared_ptr< pdat::CellVariable<double> > cast_var(
      BOOST_CAST<pdat::CellVariable<double>, hier::Variable>( var ) );
   if ( cast_var && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

const string&
QuatWeightedAverage::getOperatorName() const
{
   return(d_name_id);
}

int QuatWeightedAverage::getOperatorPriority() const
{
   return(0);
}

hier::IntVector
QuatWeightedAverage::getStencilWidth(const tbox::Dimension& dim) const {
   return(hier::IntVector(dim,0));
}

void QuatWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   boost::shared_ptr< pdat::CellData<double> > fdata ( 
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>( fine.getPatchData(src_component) ) );
   boost::shared_ptr< pdat::CellData<double> > cdata (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>( coarse.getPatchData(dst_component) ) );
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(fdata);
   assert(cdata);
   assert(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const boost::shared_ptr<geom::CartesianPatchGeometry > fgeom (
      BOOST_CAST<geom::CartesianPatchGeometry , hier::PatchGeometry>(fine.getPatchGeometry()) );
   const boost::shared_ptr<geom::CartesianPatchGeometry > cgeom (
      BOOST_CAST<geom::CartesianPatchGeometry , hier::PatchGeometry>(coarse.getPatchGeometry()) );

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   
   const int sym_flag = (int) d_symmetry_aware;

   if ( d_symmetry_aware ) {
      boost::shared_ptr< pdat::SideData<int> > rotation_index (
         BOOST_CAST< pdat::SideData<int>, hier::PatchData>(coarse.getPatchData( d_quat_symm_rotation_id) ) );
      assert( rotation_index );
      const hier::Box & rot_gbox = rotation_index->getGhostBox();
      const hier::Index r_lower = rot_gbox.lower();
      const hier::Index r_upper = rot_gbox.upper();

#if (NDIM == 2)
      quatcoarsen_(
         ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
         filo(0),filo(1),fihi(0),fihi(1),
         cilo(0),cilo(1),cihi(0),cihi(1),
         &ratio[0],
         fgeom->getDx(),
         cgeom->getDx(),
         fdata->getPointer(),
         cdata->getPointer(),
         cdata->getDepth(),
         sym_flag,
         rotation_index->getPointer(0),
         rotation_index->getPointer(1),
         r_lower[0], r_upper[0],
         r_lower[1], r_upper[1]
         );
#else
      quatcoarsen_(
         ifirstc(0),ifirstc(1),ifirstc(2),
         ilastc(0),ilastc(1),ilastc(2),
         filo(0),filo(1),filo(2),
         fihi(0),fihi(1),fihi(2),
         cilo(0),cilo(1),cilo(2),
         cihi(0),cihi(1),cihi(2),
         &ratio[0],
         fgeom->getDx(),
         cgeom->getDx(),
         fdata->getPointer(),
         cdata->getPointer(),
         cdata->getDepth(),
         sym_flag,
         rotation_index->getPointer(0),
         rotation_index->getPointer(1),
         rotation_index->getPointer(2),
         r_lower[0], r_upper[0],
         r_lower[1], r_upper[1],
         r_lower[2], r_upper[2]
         );
#endif
   }
   else {
#if (NDIM == 2)
      quatcoarsen_(
         ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
         filo(0),filo(1),fihi(0),fihi(1),
         cilo(0),cilo(1),cihi(0),cihi(1),
         &ratio[0],
         fgeom->getDx(),
         cgeom->getDx(),
         fdata->getPointer(),
         cdata->getPointer(),
         cdata->getDepth(),
         sym_flag,
         0, 0, 0, 0, 0, 0 );
#else
      quatcoarsen_(
         ifirstc(0),ifirstc(1),ifirstc(2),
         ilastc(0),ilastc(1),ilastc(2),
         filo(0),filo(1),filo(2),
         fihi(0),fihi(1),fihi(2),
         cilo(0),cilo(1),cilo(2),
         cihi(0),cihi(1),cihi(2),
         &ratio[0],
         fgeom->getDx(),
         cgeom->getDx(),
         fdata->getPointer(),
         cdata->getPointer(),
         cdata->getDepth(),
         sym_flag,
         0, 0, 0, 0, 0, 0, 0, 0, 0 );
#endif
   }
}
