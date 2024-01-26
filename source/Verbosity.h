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
#ifndef included_Verbosity
#define included_Verbosity

class Verbosity
{
 public:
   Verbosity(){};
   virtual ~Verbosity(){};
   virtual bool isSilent(void) = 0;
   virtual bool notSilent(void) = 0;
   virtual void setBasicLevel(int l) = 0;
   virtual void setAmrLevel(int l) = 0;
   virtual int basicLevel(void) = 0;
   virtual int amrLevel(void) = 0;
};

class QuatVerbosity : public Verbosity
{
 public:
   QuatVerbosity(void) : i_basic_level(1), i_amr_level(0){};
   ~QuatVerbosity(){};
   bool isSilent(void) { return basicLevel() == 0; }
   bool notSilent(void) { return !isSilent(); }
   void setBasicLevel(int l) { i_basic_level = l; }
   void setAmrLevel(int l)
   {
      i_amr_level = l;
      if (isSilent()) i_amr_level = 0;
   }
   int basicLevel(void) { return i_basic_level; }
   int amrLevel(void) { return i_amr_level; }

 private:
   int i_basic_level;
   int i_amr_level;
};

#endif
