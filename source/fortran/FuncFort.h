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
#ifndef FuncFort_H
#define FuncFort_H

#include "fc_functions_mangle.h"

// Function argument list interfaces
extern "C" {
double WELL_FUNC(const double&, const char*);
double DERIV_WELL_FUNC(const double&, const char*);
double SECOND_DERIV_WELL_FUNC(const double&, const char*);
double INTERP_FUNC(const double&, const char*);
double DERIV_INTERP_FUNC(const double&, const char*);
double SECOND_DERIV_INTERP_FUNC(const double&, const char*);
double INTERP_RATIO_FUNC(const double&, const char*, const char*);
double COMPL_INTERP_RATIO_FUNC(const double&, const char*, const char*);
double TRIPLE_WELL_FUNC(const double&, const double&, const double&);
double DERIV_TRIPLE_WELL_FUNC(const double&, const double&, const double&);
}

#endif
