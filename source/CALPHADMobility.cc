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
#include "CALPHADMobility.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#include <map>


const double m2toum2 = 1.e12;

void CALPHADMobility::initialize(std::shared_ptr<tbox::Database> db)
{
   std::map<char, short> index_map;
   index_map.insert(std::pair<char, short>('A', 0));
   index_map.insert(std::pair<char, short>('B', 1));
   index_map.insert(std::pair<char, short>('C', 2));

   std::vector<std::string> smq;
   smq.push_back("qA");
   smq.push_back("qB");
   smq.push_back("qC");
   smq.push_back("qD");

   // count first order terms to get number of species
   d_nspecies = 0;
   for (std::vector<std::string>::const_iterator it = smq.begin();
        it != smq.end(); it++) {
      if (db->keyExists(*it)) {
         d_nspecies++;
      }
   }

   d_q.resize(d_nspecies);
   for (std::vector<std::string>::const_iterator it = smq.begin();
        it != smq.end(); it++) {
      if (db->keyExists(*it)) {
         size_t n = db->getArraySize(*it);
         assert(n >= 2);
         assert(n <= 4);

         double tmp[4];
         db->getDoubleArray(*it, &tmp[0], n);
         double_pair tmpp(tmp[0], tmp[1]);
         d_q[index_map[*(it->c_str() + 1)]] = tmpp;

         if (n > 2) {
            if (d_q_extra.empty()) d_q_extra.resize(d_nspecies);
            double_pair tmpp_extra(tmp[2], tmp[3]);
            d_q_extra[index_map[*(it->c_str() + 1)]] = tmpp_extra;
         }
      }
   }

   std::vector<std::string> smqq;
   smqq.push_back("q0AB");
   if (d_nspecies > 2) {
      smqq.push_back("q0AC");
      smqq.push_back("q0BC");
   }
   //... more codes for more than 3 species

   d_qq0.resize(d_nspecies);
   for (short i = 0; i < d_nspecies; i++)
      d_qq0[i].resize(d_nspecies);

   for (std::vector<std::string>::const_iterator it = smqq.begin();
        it != smqq.end(); it++) {
      double tmp[2] = {0., 1.};
      if (db->keyExists(*it)) {
         db->getDoubleArray(*it, &tmp[0], 2);
      }
      double_pair tmpp(tmp[0], tmp[1]);
      d_qq0[index_map[*(it->c_str() + 2)]][index_map[*(it->c_str() + 3)]] =
          tmpp;
   }

   smqq.clear();
   smqq.push_back("q1AB");
   if (d_nspecies > 2) {
      smqq.push_back("q1AC");
      smqq.push_back("q1BC");
   }
   //... more codes for more than 3 species

   d_qq1.resize(d_nspecies);
   for (short i = 0; i < d_nspecies; i++)
      d_qq1[i].resize(d_nspecies);

   for (std::vector<std::string>::const_iterator it = smqq.begin();
        it != smqq.end(); it++) {
      double tmp[2] = {0., 1.};
      if (db->keyExists(*it)) {
         db->getDoubleArray(*it, &tmp[0], 2);
      }
      double_pair tmpp(tmp[0], tmp[1]);
      d_qq1[index_map[*(it->c_str() + 2)]][index_map[*(it->c_str() + 3)]] =
          tmpp;
   }

   smqq.clear();
   smqq.push_back("q2AB");
   if (d_nspecies > 2) {
      smqq.push_back("q2AC");
      smqq.push_back("q2BC");
   }
   //... more codes for more than 3 species

   d_qq2.resize(d_nspecies);
   for (short i = 0; i < d_nspecies; i++)
      d_qq2[i].resize(d_nspecies);

   for (std::vector<std::string>::const_iterator it = smqq.begin();
        it != smqq.end(); it++) {
      double tmp[2] = {0., 1.};
      if (db->keyExists(*it)) {
         db->getDoubleArray(*it, &tmp[0], 2);
      }
      double_pair tmpp(tmp[0], tmp[1]);
      d_qq2[index_map[*(it->c_str() + 2)]][index_map[*(it->c_str() + 3)]] =
          tmpp;
   }

   smqq.clear();
   smqq.push_back("q3AB");
   if (d_nspecies > 2) {
      smqq.push_back("q3AC");
      smqq.push_back("q3BC");
   }
   //... more codes for more than 3 species

   d_qq3.resize(d_nspecies);
   for (short i = 0; i < d_nspecies; i++)
      d_qq3[i].resize(d_nspecies);

   for (std::vector<std::string>::const_iterator it = smqq.begin();
        it != smqq.end(); it++) {
      double tmp[2] = {0., 1.};
      if (db->keyExists(*it)) {
         db->getDoubleArray(*it, &tmp[0], 2);
      }
      double_pair tmpp(tmp[0], tmp[1]);
      d_qq3[index_map[*(it->c_str() + 2)]][index_map[*(it->c_str() + 3)]] =
          tmpp;
   }
}

//=======================================================================

double CALPHADMobility::getDeltaG(const double c0, const double c1,
                                  const double temp) const
{
   const double dc = c0 - c1;

   return c0 * getQ(0, temp) + c1 * getQ(1, temp) +
          c0 * c1 *
              (getQQ0(0, 1, temp) +
               dc * (getQQ1(0, 1, temp) +
                     dc * (getQQ2(0, 1, temp) + dc * getQQ3(0, 1, temp))));
}

double CALPHADMobility::getDeltaG(const double c0, const double c1,
                                  const double c2, const double temp) const
{
   return c0 * getQ(0, temp) + c1 * getQ(1, temp) + c2 * getQ(2, temp) +
          c0 * c1 * (getQQ0(0, 1, temp) + (c0 - c1) * getQQ1(0, 1, temp)) +
          c0 * c2 * (getQQ0(0, 2, temp) + (c0 - c2) * getQQ1(0, 2, temp)) +
          c1 * c2 * (getQQ0(1, 2, temp) + (c1 - c2) * getQQ1(1, 2, temp));
}

//=======================================================================

void CALPHADMobility::printDiffusionVsTemperature(const double tempmin,
                                                  const double tempmax,
                                                  std::ostream& os) const
{
   const int npts = 100;
   const double dtemp = (tempmax - tempmin) / (double)(npts - 1);

   for (int i = 0; i < npts; i++) {
      const double temperature = tempmin + i * dtemp;
      const double m = getAtomicMobility(1., 0., temperature);
      os << 10000. / temperature << "\t" << m * temperature * gas_constant
         << std::endl;
   }
}

//=======================================================================
// Meolans et al., Comput. Coupling of Phase Diagrams and Thermo. 32, p.268
// (2008), eq.(71)

double computeDiffusionMobilityBinaryPhase(
    const double c0, const double temp,
    std::vector<CALPHADMobility>& calphad_mobilities_phase)
{
   assert(calphad_mobilities_phase.size() == 2);
   assert(c0 >= 0.);
   assert(c0 <= 1.);

   const double c1 = 1. - c0;
   const double m0 =
       calphad_mobilities_phase[0].getAtomicMobility(c0, c1, temp);
   const double m1 =
       calphad_mobilities_phase[1].getAtomicMobility(c0, c1, temp);
   const double mm = c0 * m1 + c1 * m0;

   assert(m0 >= 0.);
   assert(m1 >= 0.);

   return c0 * c1 * mm * m2toum2;
}

void computeDiffusionMobilityTernaryPhase(
    const double c0, const double c1, const double temp,
    std::vector<CALPHADMobility>& calphad_mobilities_phase,
    std::vector<double>& mobility)
{
   assert(calphad_mobilities_phase.size() == 3);

   const double c2 = 1. - c0 - c1;
   const double m0 =
       calphad_mobilities_phase[0].getAtomicMobility(c0, c1, c2, temp);
   const double m1 =
       calphad_mobilities_phase[1].getAtomicMobility(c0, c1, c2, temp);
   const double m2 =
       calphad_mobilities_phase[2].getAtomicMobility(c0, c1, c2, temp);

   mobility[0] = c0 *
                 ((1. - c0) * (1. - c0) * m0 + c0 * c1 * m1 + c0 * c2 * m2) *
                 m2toum2;
   mobility[1] = mobility[2] =
       c0 * c1 * (-(1. - c0) * m0 - (1. - c1) * m1 + c2 * m2) * m2toum2;
   mobility[3] = c1 *
                 (c0 * c1 * m0 + (1. - c1) * (1. - c1) * m1 + c1 * c2 * m2) *
                 m2toum2;
}
