// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "RigidBodyForces.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

using namespace SAMRAI;

RigidBodyForces::RigidBodyForces(const QuatModelParameters& model_parameters,
                                 const int phase_id, const int weight_id)
    : d_model_parameters(model_parameters),
      d_phase_id(phase_id),
      d_weight_id(weight_id)
{
   // forces between solid particles only, not liquid phase (last order
   // parameter)
   const int nforces = d_model_parameters.norderpA();
   d_forces.resize(nforces);
   for (auto& f : d_forces)
      f.resize(nforces);
   // Get timers
   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_eval_pair_forces_timer =
       tman->getTimer("AMPE::RigidBodyForces::evaluatePairForces()");
}

void RigidBodyForces::evaluatePairForces(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   assert(d_phase_id >= 0);
   assert(!d_forces.empty());

   t_eval_pair_forces_timer->start();

#if (NDIM == 2)
   std::array<double, NDIM> zero{0., 0.};
#else
   std::array<double, NDIM> zero{0., 0., 0.};
#endif
   for (auto& fi : d_forces)
      for (auto& fij : fi)
         fij = zero;

   std::vector<double> forces(d_forces.size() * d_forces.size() * NDIM, 0.);
   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         evaluatePairForces(patch, forces);
      }
   }

   const double siffness = d_model_parameters.rbStiffness();
   int count = 0;
   for (auto& fi : d_forces)
      for (auto& fij : fi)
         for (auto& f : fij) {
            f = forces[count] * siffness;
            count++;
         }

   std::vector<double> tmp;
   for (auto& fi : d_forces)
      for (auto& fij : fi)
         for (auto& f : fij)
            tmp.push_back(f);

   sumReduction(&tmp[0], tmp.size());

   double* ptmp = &tmp[0];
   for (auto& fi : d_forces)
      for (auto& fij : fi)
         for (auto& f : fij) {
            f = *ptmp;
            ptmp++;
         }

   t_eval_pair_forces_timer->stop();
}

void RigidBodyForces::printPairForces(std::ostream& os)
{
   const double threshold = 1.e-8;

   os << "## Grain pair forces:" << std::endl;
   int counti = 0;
   int countj = 0;
   for (auto& fi : d_forces) {
      os << "RB Force on particle " << counti << " :";
      for (auto& fij : fi) {
         // print only forces larger than threshold
         bool flag = false;
         for (auto& f : fij)
            if (std::abs(f) > threshold) flag = true;
         if (flag) {
            os << " (" << countj << ") ";
            for (auto& f : fij)
               os << " " << f;
         }
         countj++;
      }
      os << std::endl;
      counti++;
      countj = 0;
   }
   os << std::endl;
}

void RigidBodyForces::getPairForce(const short i, const short j,
                                   std::array<double, NDIM>& f)
{
   f = d_forces[i][j];
}
