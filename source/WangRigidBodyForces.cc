// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "WangRigidBodyForces.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

WangRigidBodyForces::WangRigidBodyForces(
    const QuatModelParameters& model_parameters, const int phase_id,
    const int conc_id, const int weight_id)
    : RigidBodyForces(model_parameters, phase_id, weight_id), d_conc_id(conc_id)
{
   tbox::plog << "WangRigidBodyForces..." << std::endl;

   d_norderp = model_parameters.norderpA();
}

void WangRigidBodyForces::evaluatePairForces(std::shared_ptr<hier::Patch> patch,
                                             std::vector<double>& forces)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_id)));
   std::shared_ptr<pdat::CellData<double> > conc(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_conc_id)));
   std::shared_ptr<pdat::CellData<double> > weight(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_weight_id)));

   assert(phase);
   assert(weight);
   assert(conc);
   assert(forces.size() == NDIM * d_norderp * d_norderp);
   assert(weight->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   // printf("orient_interp_func_type2 =%s\n",
   // d_model_parameters.orient_interp_func_type2().c_str()[0]);
   WANG_FORCES(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
               ifirst(2), ilast(2),
#endif
               dx, phase->getPointer(), phase->getGhostCellWidth()[0],
               d_norderp, conc->getPointer(), conc->getGhostCellWidth()[0],
               conc->getDepth(), d_model_parameters.rbEquilGB(),
               d_model_parameters.rbThreshold(), weight->getPointer(),
               forces.data());
}
