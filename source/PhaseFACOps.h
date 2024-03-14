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
#ifndef included_PhaseFACOps
#define included_PhaseFACOps

#include "EllipticFACOps.h"
#include "InterpolationType.h"

#include <string>

using namespace SAMRAI;

class PhaseFACOps : public EllipticFACOps
{

 public:
   PhaseFACOps(const std::string& object_name,
               const std::shared_ptr<tbox::Database> database =
                   std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(const int phase_id, const int phase_mobility_id,
                                const double epsilon_phase, const double gamma,
                                const double phase_well_scale,
                                const std::string phase_well_func_type);

   void setC(const int phi_id, const double gamma, const double phi_well_scale,
             const std::string phi_well_func_type);

 private:
   // Timers
   std::shared_ptr<tbox::Timer> t_setcoeffs_timer;

   void setCOnPatchPrivate(std::shared_ptr<pdat::CellData<double> > cd_phi,
                           std::shared_ptr<pdat::CellData<double> > cd_m,
                           std::shared_ptr<pdat::CellData<double> > cd_c,
                           const double gamma, const double phi_well_scale,
                           const char* phi_well_func_type,
                           const hier::Box& pbox);
};

#endif  // included_PhaseFACOps
