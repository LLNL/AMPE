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
#ifndef CompositionDiffusionStrategy_H
#define CompositionDiffusionStrategy_H

#include "SAMRAI/hier/PatchHierarchy.h"

using namespace SAMRAI;

enum class DiffusionInterpolationType { LINEAR, PBG, BIASED, UNDEFINED };

class CompositionDiffusionStrategy
{
 public:
   CompositionDiffusionStrategy(DiffusionInterpolationType interp_func_type)
       : d_interp_func_type(interp_func_type){};

   /*
    * compute actual diffusion by weighting diffusion in each phase
    * using phase variable
    */
   virtual void setDiffusion(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id) = 0;

   char interpChar() const
   {
      switch (d_interp_func_type) {
         case DiffusionInterpolationType::LINEAR: return 'l';
         case DiffusionInterpolationType::PBG: return 'p';
         case DiffusionInterpolationType::BIASED: return 'b';
         default:
            TBOX_ERROR(
                "Invalid interp_func_type for CompositionDiffusionStrategy");
            return 'u';
      }
   }

 private:
   DiffusionInterpolationType d_interp_func_type;
};

#endif
