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
#include "QuatWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/Utilities.h"

#include "fc_internal_mangle.h"

#include <cassert>


/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN  routines.                          *
*                                                                       *
*************************************************************************
*/
extern "C" {
#if (NDIM == 2)
void QUATCOARSEN(const int&, const int&, const int&, const int&, const int&,
                 const int&, const int&, const int&, const int&, const int&,
                 const int&, const int&, const int*, const double*,
                 const double*, const double*, double*, const int&, const int&,
                 const int*, const int*, const int&, const int&, const int&,
                 const int&);
#endif
#if (NDIM == 3)
void QUATCOARSEN(const int&, const int&, const int&, const int&, const int&,
                 const int&, const int&, const int&, const int&, const int&,
                 const int&, const int&, const int&, const int&, const int&,
                 const int&, const int&, const int&, const int*, const double*,
                 const double*, const double*, double*, const int&, const int&,
                 const int*, const int*, const int*, const int&, const int&,
                 const int&, const int&, const int&, const int&);
#endif
}

using namespace SAMRAI;

QuatWeightedAverage::QuatWeightedAverage(const bool symmetry_aware,
                                         const int quat_symm_rotation_id)
    : hier::CoarsenOperator("BASE_QUAT_COARSEN")
{
   d_name_id = "QUAT_COARSEN";
   d_symmetry_aware = symmetry_aware;
   d_quat_symm_rotation_id = quat_symm_rotation_id;
}

QuatWeightedAverage::~QuatWeightedAverage() {}

bool QuatWeightedAverage::findCoarsenOperator(
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

const std::string& QuatWeightedAverage::getOperatorName() const
{
   return (d_name_id);
}

int QuatWeightedAverage::getOperatorPriority() const { return (0); }

hier::IntVector QuatWeightedAverage::getStencilWidth(
    const tbox::Dimension& dim) const
{
   return (hier::IntVector(dim, 0));
}

void QuatWeightedAverage::coarsen(hier::Patch& coarse, const hier::Patch& fine,
                                  const int dst_component,
                                  const int src_component,
                                  const hier::Box& coarse_box,
                                  const hier::IntVector& ratio) const
{
   std::shared_ptr<pdat::CellData<double> > fdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           fine.getPatchData(src_component)));
   std::shared_ptr<pdat::CellData<double> > cdata(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           coarse.getPatchData(dst_component)));
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(fdata);
   assert(cdata);
   assert(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const std::shared_ptr<geom::CartesianPatchGeometry> fgeom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(fine.getPatchGeometry()));
   const std::shared_ptr<geom::CartesianPatchGeometry> cgeom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(coarse.getPatchGeometry()));

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   const int sym_flag = (int)d_symmetry_aware;

   if (d_symmetry_aware) {
      std::shared_ptr<pdat::SideData<int> > rotation_index(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
              coarse.getPatchData(d_quat_symm_rotation_id)));
      assert(rotation_index);
      const hier::Box& rot_gbox = rotation_index->getGhostBox();
      const hier::Index r_lower = rot_gbox.lower();
      const hier::Index r_upper = rot_gbox.upper();

#if (NDIM == 2)
      QUATCOARSEN(ifirstc(0), ifirstc(1), ilastc(0), ilastc(1), filo(0),
                  filo(1), fihi(0), fihi(1), cilo(0), cilo(1), cihi(0), cihi(1),
                  &ratio[0], fgeom->getDx(), cgeom->getDx(),
                  fdata->getPointer(), cdata->getPointer(), cdata->getDepth(),
                  sym_flag, rotation_index->getPointer(0),
                  rotation_index->getPointer(1), r_lower[0], r_upper[0],
                  r_lower[1], r_upper[1]);
#else
      QUATCOARSEN(ifirstc(0), ifirstc(1), ifirstc(2), ilastc(0), ilastc(1),
                  ilastc(2), filo(0), filo(1), filo(2), fihi(0), fihi(1),
                  fihi(2), cilo(0), cilo(1), cilo(2), cihi(0), cihi(1), cihi(2),
                  &ratio[0], fgeom->getDx(), cgeom->getDx(),
                  fdata->getPointer(), cdata->getPointer(), cdata->getDepth(),
                  sym_flag, rotation_index->getPointer(0),
                  rotation_index->getPointer(1), rotation_index->getPointer(2),
                  r_lower[0], r_upper[0], r_lower[1], r_upper[1], r_lower[2],
                  r_upper[2]);
#endif
   } else {
#if (NDIM == 2)
      QUATCOARSEN(ifirstc(0), ifirstc(1), ilastc(0), ilastc(1), filo(0),
                  filo(1), fihi(0), fihi(1), cilo(0), cilo(1), cihi(0), cihi(1),
                  &ratio[0], fgeom->getDx(), cgeom->getDx(),
                  fdata->getPointer(), cdata->getPointer(), cdata->getDepth(),
                  sym_flag, 0, 0, 0, 0, 0, 0);
#else
      QUATCOARSEN(ifirstc(0), ifirstc(1), ifirstc(2), ilastc(0), ilastc(1),
                  ilastc(2), filo(0), filo(1), filo(2), fihi(0), fihi(1),
                  fihi(2), cilo(0), cilo(1), cilo(2), cihi(0), cihi(1), cihi(2),
                  &ratio[0], fgeom->getDx(), cgeom->getDx(),
                  fdata->getPointer(), cdata->getPointer(), cdata->getDepth(),
                  sym_flag, 0, 0, 0, 0, 0, 0, 0, 0, 0);
#endif
   }
}
