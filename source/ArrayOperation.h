// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_ArrayOperation
#define included_ArrayOperation

template <class TYPE>
struct MultiplyOperation {
   void operator()(TYPE& vdst, const TYPE& vsrc) const;
};

template <class TYPE>
struct AddOperation {
   void operator()(TYPE& vdst, const TYPE& vsrc) const;
};

#include "ArrayOperation.cc"

#endif
