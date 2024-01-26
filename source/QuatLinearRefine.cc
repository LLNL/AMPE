// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
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

#include "fc_internal_mangle.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
#if (NDIM == 2)
void QUATLINREFINE2D(const int&, const int&, const int&, const int&, const int&,
                     const int&, const int&, const int&, const int&, const int&,
                     const int&, const int&, const int*, const double*,
                     const double*, const double*, double*, const int&,
                     const int*, const int*, const int&, const int&, const int&,
                     const int&);
#endif
#if (NDIM == 3)
void QUATLINREFINE3D(const int&, const int&, const int&, const int&, const int&,
                     const int&, const int&, const int&, const int&, const int&,
                     const int&, const int&, const int&, const int&, const int&,
                     const int&, const int&, const int&, const int*,
                     const double*, const double*, const double*, double*,
                     const int&, const int*, const int*, const int*, const int&,
                     const int&, const int&, const int&, const int&,
                     const int&);
#endif
}

using namespace SAMRAI;

QuatLinearRefine::QuatLinearRefine(const int quat_symm_rotation_id)
    : hier::RefineOperator("BASE_QUAT_LINEAR_REFINE")
{
   d_name_id = "QUAT_LINEAR_REFINE";
   d_quat_symm_rotation_id = quat_symm_rotation_id;
}

QuatLinearRefine::~QuatLinearRefine() {}

bool QuatLinearRefine::findRefineOperator(
    const std::shared_ptr<hier::Variable>& var,
    const std::string& op_name) const
{
   const std::shared_ptr<pdat::CellVariable<double> > cast_var(
       SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(var));
   if (cast_var && (op_name == d_name_id)) {
      return (true);
   } else {
      return (false);
   }
}

const std::string& QuatLinearRefine::getOperatorName() const
{
   return (d_name_id);
}

int QuatLinearRefine::getOperatorPriority() const { return (0); }

hier::IntVector QuatLinearRefine::getStencilWidth(
    const tbox::Dimension& dim) const
{
   return (hier::IntVector(dim, 1));
}

void QuatLinearRefine::refine(hier::Patch& fine, const hier::Patch& coarse,
                              const int dst_component, const int src_component,
                              const hier::BoxOverlap& fine_overlap,
                              const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
       dynamic_cast<const pdat::CellOverlap*>(&fine_overlap);

   TBOX_ASSERT(t_overlap != nullptr);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b = boxes.begin(); b != boxes.end();
        ++b) {
      hier::Box fine_box(*b);
      fine_box.growUpper(hier::IntVector(ratio.getDim(), -1));
      refine(fine, coarse, dst_component, src_component, fine_box, ratio);
   }
}

void QuatLinearRefine::refine(hier::Patch& fine, const hier::Patch& coarse,
                              const int dst_component, const int src_component,
                              const hier::Box& fine_box,
                              const hier::IntVector& ratio) const
{
   // const tbox::Dimension& dim(fine.getDim());

   std::shared_ptr<pdat::CellData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           coarse.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           fine.getPatchData(dst_component)));

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_ASSERT(cdata->getDepth() == 4);
   TBOX_ASSERT_OBJDIM_EQUALITY3(fine, coarse, ratio);

   const hier::Box& cgbox = cdata->getGhostBox();

   const hier::Index& cilo = cgbox.lower();
   const hier::Index& cihi = cgbox.upper();
   const hier::Index& filo = fdata->getGhostBox().lower();
   const hier::Index& fihi = fdata->getGhostBox().upper();

   const std::shared_ptr<geom::CartesianPatchGeometry> cgeom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(coarse.getPatchGeometry()));
   const std::shared_ptr<geom::CartesianPatchGeometry> fgeom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(fine.getPatchGeometry()));

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index& ifirstf = fine_box.lower();
   const hier::Index& ilastf = fine_box.upper();

   std::shared_ptr<pdat::SideData<int> > rotation_index(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
           coarse.getPatchData(d_quat_symm_rotation_id)));
   assert(rotation_index);
   const hier::Box& rot_gbox = rotation_index->getGhostBox();
   const hier::Index& r_lower = rot_gbox.lower();
   const hier::Index& r_upper = rot_gbox.upper();

#if (NDIM == 3)
   QUATLINREFINE3D(ifirstf(0), ifirstf(1), ifirstf(2), ilastf(0), ilastf(1),
                   ilastf(2), cilo(0), cilo(1), cilo(2), cihi(0), cihi(1),
                   cihi(2), filo(0), filo(1), filo(2), fihi(0), fihi(1),
                   fihi(2), &ratio[0], cgeom->getDx(), fgeom->getDx(),
                   cdata->getPointer(), fdata->getPointer(), fdata->getDepth(),
                   rotation_index->getPointer(0), rotation_index->getPointer(1),
                   rotation_index->getPointer(2), r_lower[0], r_upper[0],
                   r_lower[1], r_upper[1], r_lower[2], r_upper[2]);
#endif
#if (NDIM == 2)
   QUATLINREFINE2D(ifirstf(0), ifirstf(1), ilastf(0), ilastf(1), cilo(0),
                   cilo(1), cihi(0), cihi(1), filo(0), filo(1), fihi(0),
                   fihi(1), &ratio[0], cgeom->getDx(), fgeom->getDx(),
                   cdata->getPointer(), fdata->getPointer(), fdata->getDepth(),
                   rotation_index->getPointer(0), rotation_index->getPointer(1),
                   r_lower[0], r_upper[0], r_lower[1], r_upper[1]);
#endif
}
