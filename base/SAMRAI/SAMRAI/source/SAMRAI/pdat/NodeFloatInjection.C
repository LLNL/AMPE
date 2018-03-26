/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Constant averaging operator for node-centered float data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_NodeFloatInjection_C
#define included_pdat_NodeFloatInjection_C

#include "SAMRAI/pdat/NodeFloatInjection.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"

/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */
extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in concoarsen1d.f:
void F77_FUNC(conavgnodeflot1d, CONAVGNODEFLOT1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const float *, float *);
// in concoarsen2d.f:
void F77_FUNC(conavgnodeflot2d, CONAVGNODEFLOT2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const float *, float *);
// in concoarsen3d.f:
void F77_FUNC(conavgnodeflot3d, CONAVGNODEFLOT3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const float *, float *);
}

namespace SAMRAI {
namespace pdat {

NodeFloatInjection::NodeFloatInjection(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "CONSTANT_COARSEN")
{

}

NodeFloatInjection::~NodeFloatInjection()
{
}

int
NodeFloatInjection::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
NodeFloatInjection::getStencilWidth() const
{
   return hier::IntVector::getZero(getDim());
}

void
NodeFloatInjection::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   boost::shared_ptr<NodeData<float> > fdata(
      fine.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<NodeData<float> > cdata(
      coarse.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_DIM_ASSERT_CHECK_ARGS5(*this, coarse, fine, coarse_box, ratio);

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      if (getDim() == tbox::Dimension(1)) {
         F77_FUNC(conavgnodeflot1d, CONAVGNODEFLOT1D) (ifirstc(0), ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (getDim() == tbox::Dimension(2)) {
         F77_FUNC(conavgnodeflot2d, CONAVGNODEFLOT2D) (ifirstc(0), ifirstc(1),
            ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (getDim() == tbox::Dimension(3)) {
         F77_FUNC(conavgnodeflot3d, CONAVGNODEFLOT3D) (ifirstc(0), ifirstc(1),
            ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else {
         TBOX_ERROR(
            "NodeFloatConstantRefine::coarsen dimension > 3 not supported"
            << std::endl);
      }
   }
}

}
}
#endif
