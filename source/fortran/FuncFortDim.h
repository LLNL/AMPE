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
#ifndef FuncFortDim_H
#define FuncFortDim_H

// Link between C/C++ and Fortran files
//       name in             name in
//      C/C++ code            Fortran code
//      ----------            ------------
#if (NDIM == 2)
#define SINGLE_INDEX_STIFFNESS single_index_stiffness2d_
#define STORAGE_INDEX_STIFFNESS storage_index_stiffness2d_
#endif
#if (NDIM == 3)
#define SINGLE_INDEX_STIFFNESS single_index_stiffness3d_
#define STORAGE_INDEX_STIFFNESS storage_index_stiffness3d_
#endif

// Function argument list interfaces
extern "C" {
int SINGLE_INDEX_STIFFNESS(const int&, const int&);
int STORAGE_INDEX_STIFFNESS(const int&, const int&);
}

#endif
