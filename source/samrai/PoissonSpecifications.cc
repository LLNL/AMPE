/*************************************************************************
 *
 * This file is adapted from the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE at https://github.com/LLNL/SAMRAI.
 *
 * Copyright:     (c) 1997-2021 Lawrence Livermore National Security, LLC
 * Description:   Specifications for the scalar Poisson equation
 *
 ************************************************************************/

#include "PoissonSpecifications.h"

void PoissonSpecifications::printClassData(std::ostream &stream) const
{
   stream << "PoissonSpecifications " << d_object_name << "\n"
          << "   D is ";
   if (d_D_id != -1) {
      stream << "variable with patch id " << d_D_id << "\n";
   } else {
      stream << "constant with value " << d_D_constant << "\n";
   }
   stream << "   M is ";
   if (d_M_id != -1) {
      stream << "variable with patch id " << d_M_id << "\n";
   } else {
      stream << "constant with value " << d_M_constant << "\n";
   }
   stream << "   C is ";
   if (d_C_zero) {
      stream << "zero\n";
   } else if (d_C_id != -1) {
      stream << "variable with patch id " << d_C_id << "\n";
   } else {
      stream << "constant with value " << d_C_constant << "\n";
   }
   return;
}
