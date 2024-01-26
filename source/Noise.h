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
#ifndef included_Noise
#define included_Noise

#include "SAMRAI/hier/Patch.h"

using namespace SAMRAI;


class Noise
{
 public:
   virtual double gen() = 0;

   void setField(std::shared_ptr<hier::Patch> patch, const int data_id,
                 const int phi_id);
};

#endif
