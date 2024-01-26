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
#ifndef FuncAvgFort_H
#define FuncAvgFort_H

#include "fc_functions_mangle.h"

// Function argument list interfaces
extern "C" {
double AVERAGE_FUNC(const double&, const double&, const char*);
double DERIV_AVERAGE_FUNC(const double&, const double&, const char*);
}

#endif
