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
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <memory>

void readParabolicData(std::shared_ptr<SAMRAI::tbox::Database> input_db,
                       const std::string datasetname, double coeff[][2]);
