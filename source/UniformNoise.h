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
#include "Noise.h"

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>

class UniformNoise : public Noise
{
 public:
   static UniformNoise* instance(unsigned seed = 0)
   {
      if (s_pinstance == 0) {
         s_pinstance = new UniformNoise(seed);
      }
      return s_pinstance;
   }

   double gen() { return (*s_gen)(); }

 private:
   UniformNoise(unsigned seed);

   static UniformNoise* s_pinstance;

   typedef boost::minstd_rand rng_type;
   typedef boost::uniform_real<> distribution_type;
   typedef boost::variate_generator<rng_type, distribution_type> gen_type;

   static std::unique_ptr<
       boost::variate_generator<rng_type, distribution_type> >
       s_gen;
};
