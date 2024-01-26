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
#include "QuatModelParameters.h"

#include "SAMRAI/tbox/InputManager.h"

#include <memory>

void preRunDiagnosticsMobilityInPhases(const double temperature,
                                       QuatModelParameters& model_parameters,
                                       std::shared_ptr<tbox::Database> db);
