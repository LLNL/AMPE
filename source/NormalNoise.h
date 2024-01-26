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

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

class NormalNoise : public Noise
{
 public:
   static void setup(const int data_id);

   static NormalNoise* instance(unsigned seed = 0)
   {
      if (s_pinstance == 0) {
         s_pinstance = new NormalNoise(seed);
      }
      return s_pinstance;
   }

   double gen() { return (*s_gen)(); }

 private:
   NormalNoise(unsigned seed);

   static NormalNoise* s_pinstance;

   typedef boost::mt19937 rng_type;
   typedef boost::normal_distribution<double> distribution_type;
   typedef boost::variate_generator<rng_type, distribution_type> gen_type;

   static std::unique_ptr<
       boost::variate_generator<rng_type, distribution_type> >
       s_gen;
};
