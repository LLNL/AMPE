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
#include "ParabolicTools.h"

void readParabolicData(std::shared_ptr<SAMRAI::tbox::Database> input_db,
                       const std::string datasetname, double coeff[][2])
{
   std::shared_ptr<SAMRAI::tbox::Database> db =
       input_db->getDatabase(datasetname);
   coeff[0][0] = db->getDouble("a0");
   coeff[0][1] = db->getDouble("a1");
   coeff[1][0] = db->getDouble("b0");
   coeff[1][1] = db->getDouble("b1");
   coeff[2][0] = db->getDouble("c0");
   coeff[2][1] = db->getDouble("c1");
}
