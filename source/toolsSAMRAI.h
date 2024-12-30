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
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"

using namespace SAMRAI;


void copyDepthSideData(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                       const int dst_id, const int dst_depth, const int src_id,
                       const int src_depth);

void copyDepthCellData(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                       const int dst_id, const int dst_depth, const int src_id,
                       const int src_depth);

void sideToCell(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const int cdata_id, const int cdepth, const int sdata_id,
                const int sdepth);

double integralDepthCellData(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, const int data_id,
    const int depth, const int weight_id);

int checkForNans(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                 const int data_id);
int checkSideDataForNans(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                         const int data_id);

void checkCellDataFieldValues(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, const int data_id);
