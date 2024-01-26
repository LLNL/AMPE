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
#ifndef included_EtaFACOps
#define included_EtaFACOps

#include "EllipticFACOps.h"
#include "InterpolationType.h"

#include <string>

using namespace SAMRAI;
#ifdef HAVE_THERMO4PFM
using namespace Thermo4PFM;
#else
using namespace ampe_thermo;
#endif

class EtaFACOps : public EllipticFACOps
{

 public:
   EtaFACOps(const std::string& object_name = std::string(),
             const std::shared_ptr<tbox::Database> database =
                 std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(
       const int phase_id, const int eta_id, const int eta_mobility_id,
       const double epsilon_eta, const double gamma,
       const EnergyInterpolationType phase_interp_func_type,
       const double eta_well_scale, const std::string eta_well_func_type);

   void setC(const int phi_id, const int eta_id, const double gamma,
             const EnergyInterpolationType phi_interp_func_type,
             const double eta_well_scale, const std::string eta_well_func_type);

 private:
   void setCOnPatchPrivate(std::shared_ptr<pdat::CellData<double> > cd_phi,
                           std::shared_ptr<pdat::CellData<double> > cd_eta,
                           std::shared_ptr<pdat::CellData<double> > cd_m,
                           std::shared_ptr<pdat::CellData<double> > cd_c,
                           const double gamma,
                           const EnergyInterpolationType phi_interp_func_type,
                           const double eta_well_scale,
                           const char* eta_well_func_type,
                           const hier::Box& pbox);
};

#endif  // included_EtaFACOps
