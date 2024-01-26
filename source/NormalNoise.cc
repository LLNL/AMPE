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
#include "NormalNoise.h"

NormalNoise* NormalNoise::s_pinstance = 0;
std::unique_ptr<boost::variate_generator<NormalNoise::rng_type,
                                         NormalNoise::distribution_type> >
    NormalNoise::s_gen = 0;

NormalNoise::NormalNoise(unsigned seed)
{
   rng_type rng(seed);
   distribution_type nd(0.0, 1.0);

   s_gen = std::unique_ptr<gen_type>(new gen_type(rng, nd));
}
