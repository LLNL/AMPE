// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef ArrayOperation_C
#define ArrayOperation_C

#include "ArrayOperation.h"

template <class TYPE>
void MultiplyOperation<TYPE>::operator()(TYPE& vdst, const TYPE& vsrc) const
{
   vdst = vdst * vsrc;
}

template <class TYPE>
void AddOperation<TYPE>::operator()(TYPE& vdst, const TYPE& vsrc) const
{
   vdst = vdst + vsrc;
}

#endif
