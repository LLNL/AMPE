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
// Neumann-Kopp rule
#ifndef included_NKRHeatCapacityStrategy
#define included_NKRHeatCapacityStrategy

#include "HeatCapacityStrategy.h"

// HeatCapacityStrategy based on Neumann-Kopp rule
// Cp(A-B)=c*Cp(A)+(1-c)*Cp(B)+
class NKRHeatCapacityStrategy : public HeatCapacityStrategy
{
 public:
   NKRHeatCapacityStrategy(const std::vector<std::map<short, double> >& cp,
                           const int cp_id, const int composition_id,
                           const int temperature_id)
       : d_cp(cp),
         d_cp_id(cp_id),
         d_composition_id(composition_id),
         d_temperature_id(temperature_id)
   {
   }
   void setCurrentValue(std::shared_ptr<hier::PatchHierarchy> patch_hierarchy);

 private:
   // cp coefficients for T polynomial expansion for each species
   std::vector<std::map<short, double> > d_cp;

   int d_cp_id;
   int d_composition_id;
   int d_temperature_id;
};

#endif
