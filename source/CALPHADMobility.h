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

const double gas_constant = GASCONSTANT_R_JPKPMOL;

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

      double val =
          d_q[isp0].first + gas_constant * temperature * log(d_q[isp0].second);

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

      return d_qq0[isp0][isp1].first +
             gas_constant * temperature * log(d_qq0[isp0][isp1].second);
   }

   double getQQ1(const unsigned short isp0, const unsigned short isp1,
                 const double temperature) const
   {
      assert(d_qq1[isp0][isp1].second > 0.);

      return d_qq1[isp0][isp1].first +
             gas_constant * temperature * log(d_qq1[isp0][isp1].second);
   }

   double getQQ2(const unsigned short isp0, const unsigned short isp1,
                 const double temperature) const
   {
      assert(d_qq2[isp0][isp1].second > 0.);

      return d_qq2[isp0][isp1].first +
             gas_constant * temperature * log(d_qq2[isp0][isp1].second);
   }

   double getQQ3(const unsigned short isp0, const unsigned short isp1,
                 const double temperature) const
   {
      assert(d_qq3[isp0][isp1].second > 0.);

      return d_qq3[isp0][isp1].first +
             gas_constant * temperature * log(d_qq3[isp0][isp1].second);
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
      const double rtinv = 1. / (gas_constant * temp);
      return exp(dG * rtinv) * rtinv;
   }

   double getAtomicMobility(const double c0, const double c1, const double c2,
                            const double temp) const
   {
      assert(temp > 0.);

      const double dG = getDeltaG(c0, c1, c2, temp);
      const double rtinv = 1. / (gas_constant * temp);
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
