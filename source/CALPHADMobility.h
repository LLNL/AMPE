// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#ifndef included_CALPHADMobility
#define included_CALPHADMobility

// Mobility data for one species

#include "PhysicalConstants.h"

#include "SAMRAI/tbox/Database.h"
using namespace SAMRAI;

#include <cmath>
#include <cassert>

#include <string>
#include <vector>

// to store coefficients of linear polynomial
typedef std::pair<double, double> double_pair;

class CALPHADMobility
{
 private:
   std::string d_name;
   short d_nspecies;

   // each "double_pair" contains the two coefficients
   // (a0,a1) for the expansion a(T)=a0+R*T*ln(a1)
   std::vector<double_pair> d_q;
   std::vector<double_pair> d_q_extra;  // for T*ln(T) and 1/T terms
   std::vector<std::vector<double_pair> > d_qq0;
   std::vector<std::vector<double_pair> > d_qq1;
   std::vector<std::vector<double_pair> > d_qq2;
   std::vector<std::vector<double_pair> > d_qq3;

   double getQ(const unsigned short isp0, const double temperature) const
   {
      assert(isp0 < (unsigned short)d_q.size());
      assert(d_q[isp0].second > 0.);

      double val = d_q[isp0].first +
                   gas_constant_R_JpKpmol * temperature * log(d_q[isp0].second);

      if (!d_q_extra.empty()) {
         val += d_q_extra[isp0].first * temperature * log(temperature) +
                d_q_extra[isp0].second / temperature;
      }
      return val;
   }

   double getQQ0(const unsigned short isp0, const unsigned short isp1,
                 const double temperature) const
   {
      assert(d_qq0[isp0][isp1].second > 0.);

      return d_qq0[isp0][isp1].first + gas_constant_R_JpKpmol * temperature *
                                           log(d_qq0[isp0][isp1].second);
   }

   double getQQ1(const unsigned short isp0, const unsigned short isp1,
                 const double temperature) const
   {
      assert(d_qq1[isp0][isp1].second > 0.);

      return d_qq1[isp0][isp1].first + gas_constant_R_JpKpmol * temperature *
                                           log(d_qq1[isp0][isp1].second);
   }

   double getQQ2(const unsigned short isp0, const unsigned short isp1,
                 const double temperature) const
   {
      assert(d_qq2[isp0][isp1].second > 0.);

      return d_qq2[isp0][isp1].first + gas_constant_R_JpKpmol * temperature *
                                           log(d_qq2[isp0][isp1].second);
   }

   double getQQ3(const unsigned short isp0, const unsigned short isp1,
                 const double temperature) const
   {
      assert(d_qq3[isp0][isp1].second > 0.);

      return d_qq3[isp0][isp1].first + gas_constant_R_JpKpmol * temperature *
                                           log(d_qq3[isp0][isp1].second);
   }

   double getDeltaG(const double c0, const double c1, const double temp) const;

   double getDeltaG(const double c0, const double c1, const double c2,
                    const double temp) const;

 public:
   CALPHADMobility(const std::string& name) { d_name = name; }

   void initialize(std::shared_ptr<tbox::Database> db);

   // returns mobility in [um^2*mol/s*J]
   double getAtomicMobility(const double c0, const double c1,
                            const double temp) const
   {
      assert(temp > 0.);
      assert(c0 >= 0.);
      assert(c1 >= 0.);
      assert(c0 <= 1.);
      assert(c1 <= 1.);

      const double dG = getDeltaG(c0, c1, temp);
      const double rtinv = 1. / (gas_constant_R_JpKpmol * temp);
      return exp(dG * rtinv) * rtinv;
   }

   double getAtomicMobility(const double c0, const double c1, const double c2,
                            const double temp) const
   {
      assert(temp > 0.);

      const double dG = getDeltaG(c0, c1, c2, temp);
      const double rtinv = 1. / (gas_constant_R_JpKpmol * temp);
      return exp(dG * rtinv) * rtinv;
   }

   double getAtomicMobilityBinary(const double c0, const double temp) const
   {
      assert(temp > 0.);

      const double c1 = 1. - c0;
      return getAtomicMobility(c0, c1, temp);
   }

   void printDiffusionVsTemperature(const double tempmin, const double tempmax,
                                    std::ostream& os) const;
};

double computeDiffusionMobilityBinaryPhase(
    const double c0, const double temp,
    std::vector<CALPHADMobility>& calphad_mobilities_phase);
void computeDiffusionMobilityTernaryPhase(
    const double c0, const double c1, const double temp,
    std::vector<CALPHADMobility>& calphad_mobilities_phase,
    std::vector<double>& mobility);

#endif
