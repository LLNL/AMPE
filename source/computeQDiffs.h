// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

void computeQGradSide(
    std::shared_ptr<SAMRAI::pdat::SideData<double> > diff_data,
    std::shared_ptr<SAMRAI::pdat::SideData<double> > grad_data,
    const double dx[NDIM], const bool isSymmetryAware,
    const bool isotropicStencil,
    std::shared_ptr<SAMRAI::pdat::SideData<int> > rotation_index);

void computeQGrad(std::shared_ptr<SAMRAI::pdat::SideData<double> > diff_data,
                  std::shared_ptr<SAMRAI::pdat::CellData<double> > grad_data,
                  const double dx[NDIM], const bool isSymmetryAware,
                  std::shared_ptr<SAMRAI::pdat::SideData<int> > rotation_index);

void computeQDiffs(
    std::shared_ptr<SAMRAI::pdat::CellData<double> > quat_data,
    std::shared_ptr<SAMRAI::pdat::SideData<double> > diff_data,
    const bool isSymmetryAware,
    std::shared_ptr<SAMRAI::pdat::SideData<int> > rotation_index);
