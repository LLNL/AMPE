// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
#ifndef included_MovingFrameRHS
#define included_MovingFrameRHS

#include "SAMRAI/hier/PatchHierarchy.h"

using namespace SAMRAI;

class MovingFrameRHS
{
 public:
   MovingFrameRHS(const int data_scratch_id, const bool upwind);

   void addRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
               const int ydot_id, const double frame_velocity);

 private:
   const int d_data_id;

#if (NDIM == 3)
   void (*d_add_vdphidx)(const int ifirst0, const int ilast0, const int ifirst1,
                         const int ilast1, const int ifirst2, const int ilast2,
                         const double* dx, const double* phi, const int ngphi,
                         const double frame_velocity, double* data_rhs,
                         const int ngrhs, const int physb);
#else
   void (*d_add_vdphidx)(const int ifirst0, const int ilast0, const int ifirst1,
                         const int ilast1, const double* dx, const double* phi,
                         const int ngphi, const double frame_velocity,
                         double* data_rhs, const int ngrhs, const int physb);
#endif
};

#endif
