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
#include "UniformNoise.h"

UniformNoise* UniformNoise::s_pinstance = 0;
std::unique_ptr<boost::variate_generator<UniformNoise::rng_type,
                                         UniformNoise::distribution_type> >
    UniformNoise::s_gen = 0;

UniformNoise::UniformNoise(unsigned seed)
{
   rng_type rng(seed);
   distribution_type nd(-0.5, 0.5);

   s_gen = std::unique_ptr<gen_type>(new gen_type(rng, nd));
}
