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
#include "QuatLinearRefine.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/CellOverlap.h"


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
   void quatlinrefine2d_(
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int*, const double*, const double*,
      const double*, double*,
      const int&,
      const int*, const int*,
      const int&, const int&, const int&, const int&
      );
#endif
#if (NDIM == 3)
   void quatlinrefine3d_(
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int&, const int&, const int&,
      const int*, const double*, const double*,
      const double*, double*,
      const int&,
      const int*, const int*, const int*,
      const int&, const int&, const int&,
      const int&, const int&, const int&
      );
#endif
}

using namespace SAMRAI;

QuatLinearRefine::
QuatLinearRefine( const int quat_symm_rotation_id )
: hier::RefineOperator("BASE_QUAT_LINEAR_REFINE")
{
   d_name_id = "QUAT_LINEAR_REFINE";
   d_quat_symm_rotation_id = quat_symm_rotation_id;
}

QuatLinearRefine::
~QuatLinearRefine()
{
}

bool QuatLinearRefine::findRefineOperator(
   const boost::shared_ptr<hier::Variable >& var,
   const string &op_name) const
{
   const boost::shared_ptr< pdat::CellVariable<double> > cast_var(
      BOOST_CAST<pdat::CellVariable<double>, hier::Variable>(var) );
   if ( cast_var && (op_name == d_name_id) ) {
      return(true);
   } else {
      return(false);
   }
}

const string&
QuatLinearRefine::getOperatorName() const
{
   return(d_name_id);
}

int QuatLinearRefine::getOperatorPriority() const
{
   return(0);
}

hier::IntVector 
QuatLinearRefine::getStencilWidth(const tbox::Dimension& dim) const {
   return(hier::IntVector(dim,1));
}

void QuatLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
      dynamic_cast<const pdat::CellOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b=boxes.begin(); b != boxes.end(); ++b) {
      hier::Box fine_box(*b);
      fine_box.growUpper(hier::IntVector(ratio.getDim(), -1));
      refine(fine,
         coarse,
         dst_component,
         src_component,
         fine_box,
         ratio);
   }
}

void QuatLinearRefine::refine(
   hier::Patch& fine, 
   const hier::Patch& coarse, 
   const int dst_component, 
   const int src_component, 
   const hier::Box& fine_box, 
   const hier::IntVector& ratio) const
{
   //const tbox::Dimension& dim(fine.getDim());

   boost::shared_ptr< pdat::CellData<double> > cdata (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>( coarse.getPatchData(src_component) ) );
   boost::shared_ptr< pdat::CellData<double> > fdata (
      BOOST_CAST< pdat::CellData<double>, hier::PatchData>( fine.getPatchData(dst_component) ) );

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_ASSERT(cdata->getDepth() == 4);
   TBOX_ASSERT_OBJDIM_EQUALITY3(fine, coarse, ratio);

   const hier::Box& cgbox=cdata->getGhostBox();

   const hier::Index& cilo = cgbox.lower();
   const hier::Index& cihi = cgbox.upper();
   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();

   const boost::shared_ptr<geom::CartesianPatchGeometry > cgeom(
                                                    BOOST_CAST<geom::CartesianPatchGeometry , hier::PatchGeometry>(coarse.getPatchGeometry()) );
   const boost::shared_ptr<geom::CartesianPatchGeometry > fgeom(
                                                    BOOST_CAST<geom::CartesianPatchGeometry , hier::PatchGeometry>(fine.getPatchGeometry()) );

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index& ifirstf = fine_box.lower();
   const hier::Index& ilastf = fine_box.upper();

   boost::shared_ptr< pdat::SideData<int> > rotation_index (
      BOOST_CAST< pdat::SideData<int>, hier::PatchData>(coarse.getPatchData( d_quat_symm_rotation_id) ) );
   assert( rotation_index );
   const hier::Box & rot_gbox = rotation_index->getGhostBox();
   const hier::Index& r_lower = rot_gbox.lower();
   const hier::Index& r_upper = rot_gbox.upper();

#if (NDIM==3)
   quatlinrefine3d_(
      ifirstf(0),ifirstf(1),ifirstf(2),
      ilastf(0),ilastf(1),ilastf(2),
      cilo(0),cilo(1),cilo(2),
      cihi(0),cihi(1),cihi(2),
      filo(0),filo(1),filo(2),
      fihi(0),fihi(1),fihi(2),
      &ratio[0],
      cgeom->getDx(),
      fgeom->getDx(),
      cdata->getPointer(),
      fdata->getPointer(),
      fdata->getDepth(),
      rotation_index->getPointer(0),
      rotation_index->getPointer(1),
      rotation_index->getPointer(2),
      r_lower[0], r_upper[0],
      r_lower[1], r_upper[1],
      r_lower[2], r_upper[2]
      );
#endif
#if (NDIM==2)
   quatlinrefine2d_(
      ifirstf(0),ifirstf(1),
      ilastf(0),ilastf(1),
      cilo(0),cilo(1),
      cihi(0),cihi(1),
      filo(0),filo(1),
      fihi(0),fihi(1),
      &ratio[0],
      cgeom->getDx(),
      fgeom->getDx(),
      cdata->getPointer(),
      fdata->getPointer(),
      fdata->getDepth(),
      rotation_index->getPointer(0),
      rotation_index->getPointer(1),
      r_lower[0], r_upper[0],
      r_lower[1], r_upper[1]
      );
#endif
}

