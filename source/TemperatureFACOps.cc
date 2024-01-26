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
#include "TemperatureFACOps.h"

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "FuncFort.h"

#include <cassert>

//======================================================================

TemperatureFACOps::TemperatureFACOps(const std::string &object_name,
                                     std::shared_ptr<tbox::Database> database)
    : EllipticFACOps(tbox::Dimension(NDIM), object_name, database)
{
   return;
}

//======================================================================
// Eq.     M div (D grad u) + C u = f

void TemperatureFACOps::setOperatorCoefficients(const double m, const double c,
                                                const double d)
{
   assert(m > 0.);
   assert(fabs(d) > 0.);
   assert(d < 1.e15);

   // tbox::pout<<"TemperatureFACOps::setOperatorCoefficients: m="<<m<<",
   // c="<<c<<", d="<<d<<endl;
   setMConstant(m);
   setCConstant(c);
   setDConstant(d);
}
