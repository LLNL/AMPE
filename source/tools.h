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
#ifndef included_Tools
#define included_Tools

#include <map>
#include <string>

void listLocalToGlobal(std::map<int, double> &lcl_map,
                       std::map<int, double> &gbl_map);

void sumReduction(double *x, const int count);
double sumReduction(const double x);
int sumReduction(const int x);
void allGatherv(const int *x_in, int size_in, int *x_out, int size_out);
void allGatherv(const double *x_in, int size_in, double *x_out, int size_out);
void allGatherSetup(int size_in, const int size_out, int *&rcounts,
                    int *&disps);
void printDeprecated(const std::string &s_old, const std::string &s_new);

#endif
